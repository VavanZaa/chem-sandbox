"""Check accepted electron-fragment energy costs against GRI gas enthalpies.

This is a feed-only diagnostic. It never loads an HCN assay and does not tune
or change any source term. Standard reaction enthalpy is a necessary energy
floor, not a plasma heat balance or a rate validation.
"""

from __future__ import annotations

import json
from typing import Any

import cantera as ct
import numpy as np
from scipy.constants import Avogadro, elementary_charge

from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.schema import PrebioticApparatusRequest
from scripts.probe_barth_feed_pathway import make_probe_payload


def audit_source_thermochemistry(
    result: dict[str, Any], temperature_K: float,
) -> dict[str, Any]:
    """Audit the necessary enthalpy floor and N2 parent-event capacity."""

    phase = ct.Solution("gri30.yaml")
    phase.TP = temperature_K, ct.one_atm
    enthalpy_J_mol = np.asarray(phase.standard_enthalpies_RT) * ct.gas_constant * temperature_K / 1000.0
    by_name = dict(zip(phase.species_names, map(float, enthalpy_J_mol), strict=True))
    faraday_J_mol_eV = Avogadro * elementary_charge
    n2_dissociation_requirement_eV = (
        2 * by_name["N"] - by_name["N2"]
    ) / faraday_J_mol_eV
    rows = []
    for step in result["plasma_model"]["heavy_species_cycle_summaries"]:
        for source in step["source_extents"]:
            reaction_enthalpy_J_mol = (
                sum(by_name[name] for name in source["products"])
                - by_name[source["target"]]
            )
            requirement_eV = max(reaction_enthalpy_J_mol, 0.0) / faraday_J_mol_eV
            source_cost_eV = float(source["energy_cost_eV_per_event"])
            rows.append({
                "channel_id": source["channel_id"],
                "target": source["target"],
                "products": source["products"],
                "standard_reaction_enthalpy_eV": reaction_enthalpy_J_mol / faraday_J_mol_eV,
                "nonnegative_standard_enthalpy_requirement_eV": requirement_eV,
                "source_energy_cost_eV": source_cost_eV,
                "outcome_only_parent_manifold": source["energy_loss_accounted_by_parent_manifold"],
                "energy_deficit_eV": max(requirement_eV - source_cost_eV, 0.0),
                "accepted_extent_mol": float(source.get("actual_integrated_extent_mol", source["accepted_extent_mol"])),
            })
    deficits = [
        row for row in rows
        if row["energy_deficit_eV"] > 0.02
        and not row["outcome_only_parent_manifold"]
    ]
    n2_outcome_parent_steps = []
    for snapshot, heavy_step in zip(
        result["plasma_model"]["cycle_summaries"],
        result["plasma_model"]["heavy_species_cycle_summaries"],
        strict=True,
    ):
        parent_raw = sum(
            float(source["uncoupled_active_window_opportunity_mol"])
            for source in snapshot["source_terms"]
            if source["target"] == "N2"
            and source["role"] == "excitation"
            and not source.get("outcome_only")
            and float(source.get("threshold_eV") or 0.0) >= n2_dissociation_requirement_eV
        )
        scale = float(heavy_step["collision_energy_allocation"]["global_collision_energy_scale"])
        outcome_accepted = sum(
            float(source.get("actual_integrated_extent_mol", source["accepted_extent_mol"]))
            for source in heavy_step["source_extents"]
            if source["target"] == "N2"
            and source["products"] == ["N", "N"]
            and source["energy_loss_accounted_by_parent_manifold"]
        )
        n2_outcome_parent_steps.append({
            "substep": snapshot["substep"],
            "high_energy_parent_excitation_energy_scaled_opportunity_mol": parent_raw * scale,
            "n2_dissociation_outcome_accepted_mol": outcome_accepted,
            "parent_event_capacity_covers_outcome": outcome_accepted <= parent_raw * scale * (1 + 1e-8),
        })
    parent_deficits = [
        row for row in n2_outcome_parent_steps
        if not row["parent_event_capacity_covers_outcome"]
    ]
    return {
        "status": (
            "accepted_source_thermochemical_floor_passed"
            if not deficits and not parent_deficits
            else "accepted_source_thermochemical_floor_failed"
        ),
        "assay_loaded": False,
        "gas_temperature_K": temperature_K,
        "standard_state_basis": "GRI-Mech ideal-gas standard enthalpy at the scenario temperature",
        "accepted_channel_instances": len(rows),
        "distinct_accepted_channels": len({row["channel_id"] for row in rows}),
        "deficits": deficits,
        "n2_outcome_parent_manifold_coverage": n2_outcome_parent_steps,
        "n2_parent_coverage_deficits": parent_deficits,
        "n2_parent_minimum_capacity_to_outcome_ratio": min(
            (
                row["high_energy_parent_excitation_energy_scaled_opportunity_mol"]
                / row["n2_dissociation_outcome_accepted_mol"]
                for row in n2_outcome_parent_steps
                if row["n2_dissociation_outcome_accepted_mol"] > 0
            ),
            default=None,
        ),
        "channels": rows,
        "material_ledger_closed": result["material_ledger"]["closed"],
        "complete_thermal_energy_balance_closed": result["energy_ledger"]["complete_thermal_energy_balance_closed"],
        "science_boundary": (
            "Reaction enthalpy is a necessary energy floor only. This does not "
            "prove that a specific N2 parent excitation dissociates. Parent-event "
            "capacity is necessary but not sufficient. It does not validate "
            "electron energy transfer, excited products, reactor heat, "
            "spark exposure, HCN identity discovery, or experimental yield."
        ),
    }


def audit_barth_feed_source_thermochemistry() -> dict[str, Any]:
    payload = make_probe_payload(
        duration_s=60, cycles=1, substeps=2,
        gas_water_kla_per_s=1e-4, gas_water_solver="coupled_bdf",
        primary_source_forcing="inventory_limited",
        formaldehyde_hydration_model="winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    )
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    return audit_source_thermochemistry(
        result, float(payload["apparatus"]["discharge_temperature"]["value"]),
    )


def main() -> None:
    report = audit_barth_feed_source_thermochemistry()
    print(json.dumps(report, indent=2))
    if report["status"] == "accepted_source_thermochemical_floor_failed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
