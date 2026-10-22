"""Short CO2-fed closed-flask pathway diagnostic, not a Barth HCN comparison.

This uses only a feed/geometry card without an HCN measurement. Spark pulse
width, active volume, coupling, temperature and pressure are explicitly
hypothetical. A result from this script cannot validate an absolute yield.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
from typing import Any
from collections import Counter

from scipy.constants import R

from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.boltzmann import IST_CO2_O2_DATASET_ID
from app.prebiotic.schema import PrebioticApparatusRequest
from scripts.run_hcn_formation_gate import PROFILE


FEED_CARD = (
    Path(__file__).resolve().parents[1] / "data" / "benchmarks"
    / "barth_2024_feed_only_220725001_v1.json"
)


def make_probe_payload(
    *, duration_s: float = 60.0, spark_on: bool = True, substeps: int = 2,
    cycles: int = 1, gas_water_kla_per_s: float | None = None,
    gas_water_solver: str = "post_step_analytic",
    primary_source_forcing: str = "time_averaged",
    formaldehyde_hydration_model: str = "off",
    remove_carbon_for_control: bool = False,
    remove_nitrogen_for_control: bool = False,
) -> dict[str, Any]:
    if not 0 < duration_s <= 3600:
        raise ValueError("This diagnostic is limited to a positive duration of at most one hour")
    if not 1 <= substeps <= 48:
        raise ValueError("This diagnostic requires 1 to 48 integration substeps")
    if not 1 <= cycles <= 4:
        raise ValueError("This diagnostic requires 1 to 4 identical cycles")
    if gas_water_kla_per_s is not None and (
        not math.isfinite(gas_water_kla_per_s) or gas_water_kla_per_s < 0
    ):
        raise ValueError("The hypothetical gas/water kLa must be finite and nonnegative")
    if gas_water_solver not in {"post_step_analytic", "coupled_bdf"}:
        raise ValueError("Unknown gas/water contact solver")
    if primary_source_forcing not in {"time_averaged", "inventory_limited"}:
        raise ValueError("Unknown distributed primary source forcing")
    if gas_water_solver == "coupled_bdf" and (gas_water_kla_per_s is None or not spark_on):
        raise ValueError("Coupled contact requires an active spark and declared kLa")
    if formaldehyde_hydration_model not in {
        "off", "winkelman_2002_source_bounded",
        "winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    }:
        raise ValueError("Unknown formaldehyde hydration model")
    if formaldehyde_hydration_model != "off" and gas_water_solver != "coupled_bdf":
        raise ValueError("Formaldehyde hydration requires coupled contact")
    if remove_carbon_for_control and remove_nitrogen_for_control:
        raise ValueError("Choose one altered-feed negative control at a time")
    feed = json.loads(FEED_CARD.read_text(encoding="utf-8"))
    base = deepcopy(json.loads(PROFILE.read_text(encoding="utf-8")))
    assumptions = feed["diagnostic_only_assumptions"]
    measured_fractions = feed["measured_initial_gas_mole_fractions"]
    n2_fraction = 1.0 - sum(float(value) for value in measured_fractions.values())
    if not 0 < n2_fraction < 1:
        raise ValueError("Measured feed fractions leave no positive N2 balance")
    gas_fractions = {**measured_fractions, "N2": n2_fraction}
    if remove_carbon_for_control:
        gas_fractions["N2"] += gas_fractions.pop("CO2")
    if remove_nitrogen_for_control:
        gas_fractions["H2"] += gas_fractions.pop("N2")
    temperature_K = float(assumptions["gas_and_water_temperature_K"])
    headspace_L = float(assumptions["headspace_volume_L"])
    total_gas_mol = (
        float(assumptions["initial_pressure_Pa"]) * headspace_L * 1e-3
        / (R * temperature_K)
    )
    base["id"] = (
        "barth_feed_carbon_removed_negative_control"
        if remove_carbon_for_control else "barth_feed_only_pathway_probe_not_yield_validation"
    )
    if remove_nitrogen_for_control:
        base["id"] = "barth_feed_nitrogen_removed_negative_control"
    base["species"].extend((
        {"id": "CO2", "name": "Carbon dioxide", "composition": {"C": 1, "O": 2},
         "charge": 0, "smiles": "O=C=O", "phase_behavior": "permanent_gas"},
        {"id": "O2", "name": "Oxygen", "composition": {"O": 2},
         "charge": 0, "smiles": "O=O", "phase_behavior": "permanent_gas"},
    ))
    base["initial_inventory"] = [
        {"species_id": "H2O", "compartment": "boiler_liquid", "amount": {
            "value": feed["reported_water_volume_L"] * assumptions["water_density_mol_per_L"],
            "unit": "mol",
        }},
        *(
            {"species_id": species, "compartment": "gas_loop", "amount": {
                "value": total_gas_mol * fraction, "unit": "mol",
            }}
            for species, fraction in gas_fractions.items()
        ),
    ]
    base["apparatus"].update({
        "cycles": cycles,
        "operation_mode": "closed_stirred_flask",
        "integration_substeps_per_cycle": substeps,
        "primary_source_forcing": primary_source_forcing,
        "aqueous_acid_base_model": "hcn_nh3_water_equilibrium",
        "boiler_liquid_volume": {"value": feed["reported_water_volume_L"], "unit": "L"},
        "gas_loop_volume": {"value": headspace_L, "unit": "L"},
        "discharge_zone_volume": {
            "value": assumptions["discharge_zone_volume_L"], "unit": "L",
            "source": "hypothetical active zone for pathway diagnostic; not Barth measurement",
            "method": "declared geometry probe", "basis": "hypothetical",
        } if spark_on else None,
        "condensate_trap_volume": {"value": feed["reported_water_volume_L"], "unit": "L"},
        "gas_circulation_rate": {"value": 0, "unit": "L/s"},
        "pressure": {"value": assumptions["initial_pressure_Pa"], "unit": "Pa"},
        "boiler_temperature": {"value": temperature_K, "unit": "K"},
        "discharge_temperature": {"value": temperature_K, "unit": "K"},
        "condenser_temperature": {"value": temperature_K, "unit": "K"},
        "boil_duration": {"value": 0, "unit": "s"},
        "discharge_duration": {"value": duration_s / cycles, "unit": "s"},
        "condense_duration": {"value": 0, "unit": "s"},
        "condenser_phase_transfer_model": "step_equilibrium",
    })
    if gas_water_kla_per_s is not None:
        base["apparatus"].update({
            "closed_flask_gas_water_transfer_model": "finite_rate_contact",
            "closed_flask_finite_contact_solver": gas_water_solver,
            "closed_flask_liquid_side_kla": {
                "value": gas_water_kla_per_s, "unit": "1/s",
                "source": "hypothetical contact-rate sensitivity; not a Barth measurement",
                "method": "declared two-compartment Henry relaxation",
                "basis": "hypothetical",
            },
        })
    if formaldehyde_hydration_model != "off":
        base["apparatus"]["closed_flask_formaldehyde_hydration"] = formaldehyde_hydration_model
    if spark_on:
        base["spark"].update({
            "voltage": {"value": 49000, "unit": "V"},
            "current": {"value": 0.001, "unit": "A"},
            "pulse_width": {"value": assumptions["pulse_width_s"], "unit": "s"},
            "frequency": {"value": 500000, "unit": "Hz"},
            "coupling_efficiency": assumptions["electrical_coupling_efficiency"],
            "source": "nominal BD-50E voltage/current/frequency; pulse width and coupling hypothetical",
            "basis": "hypothetical",
        })
        base["spark"]["boltzmann"].update({
            "collision_set": IST_CO2_O2_DATASET_ID,
            "reduced_electric_field_Td": assumptions["reduced_electric_field_Td"],
        })
    else:
        base["spark"] = {"mode": "off"}
    return base


def run_probe(
    *, duration_s: float = 60.0, spark_on: bool = True, substeps: int = 2,
    cycles: int = 1, gas_water_kla_per_s: float | None = None,
    gas_water_solver: str = "post_step_analytic",
    primary_source_forcing: str = "time_averaged",
    formaldehyde_hydration_model: str = "off",
    remove_carbon_for_control: bool = False,
    remove_nitrogen_for_control: bool = False,
) -> dict[str, Any]:
    payload = make_probe_payload(
        duration_s=duration_s, spark_on=spark_on, substeps=substeps, cycles=cycles,
        gas_water_kla_per_s=gas_water_kla_per_s,
        gas_water_solver=gas_water_solver,
        primary_source_forcing=primary_source_forcing,
        formaldehyde_hydration_model=formaldehyde_hydration_model,
        remove_carbon_for_control=remove_carbon_for_control,
        remove_nitrogen_for_control=remove_nitrogen_for_control,
    )
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    branch = result["plasma_model"].get("propagated_charged_branch") or {}
    heavy_steps = result["plasma_model"]["heavy_species_cycle_summaries"]
    o1d_fates = [
        item["o1d_fate_diagnostic"]
        for item in result["plasma_model"]["cycle_summaries"]
        if item.get("o1d_fate_diagnostic") is not None
    ]
    accepted_neutral_sources = [
        item for step in heavy_steps for item in step["source_extents"]
    ]
    excluded = [
        channel
        for step in heavy_steps
        for channel in step["collision_energy_allocation"]["excluded_fragment_channels"]
    ]
    charged_compartments = branch.get("final_compartments_mol") or {}
    charged_cyanide_by_compartment = {
        place: float(amounts.get("HCN", 0.0)) + float(amounts.get("CN-", 0.0))
        for place, amounts in charged_compartments.items()
    }
    charged_roundoff_clip_mol = float(
        branch.get("unprojected_ion_roundoff_clip_total_mol") or 0.0
    )
    charged_controls_pass = (
        not spark_on or (
            (branch.get("material_ledger") or {}).get("closed") is True
            and branch.get("failed_substep_count") == 0
            and charged_roundoff_clip_mol <= 1e-15
        )
    )
    return {
        "status": (
            "carbon_removed_negative_control_not_yield_validation"
            if remove_carbon_for_control and result["material_ledger"]["closed"] and charged_controls_pass
            else "nitrogen_removed_negative_control_not_yield_validation"
            if remove_nitrogen_for_control and result["material_ledger"]["closed"] and charged_controls_pass
            else "feed_matched_short_pathway_diagnostic_not_yield_validation"
            if result["material_ledger"]["closed"] and charged_controls_pass
            else "pathway_probe_failed_material_or_charged_controls"
        ),
        "duration_s": duration_s,
        "spark_on": spark_on,
        "carbon_removed_negative_control": remove_carbon_for_control,
        "nitrogen_removed_negative_control": remove_nitrogen_for_control,
        "integration_substeps_per_cycle": substeps,
        "cycles": cycles,
        "total_integration_steps": cycles * substeps,
        "gas_water_contact_model": result["operating_conditions"]["closed_flask_gas_water_transfer_model"],
        "gas_water_contact_solver": result["operating_conditions"]["closed_flask_finite_contact_solver"],
        "primary_source_forcing": result["operating_conditions"]["primary_source_forcing"],
        "formaldehyde_hydration_model": result["operating_conditions"]["closed_flask_formaldehyde_hydration"],
        "formaldehyde_hydration_minimum_ideal_pH": result["operating_conditions"]["formaldehyde_hydration_minimum_ideal_pH"],
        "formaldehyde_hydration_pH_extrapolation_used": result["scientific_firewall"]["formaldehyde_hydration_pH_extrapolation_used"],
        "hypothetical_liquid_side_kla_per_s": gas_water_kla_per_s,
        "initial_feed_mol": {
            item["species_id"]: item["amount"]["value"]
            for item in payload["initial_inventory"] if item["compartment"] == "gas_loop"
        },
        "final_gas_CO_mol": result["final_compartments_mol"]["gas_loop"].get("CO", 0.0),
        "final_gas_NO_mol": result["final_compartments_mol"]["gas_loop"].get("NO", 0.0),
        "final_gas_HCN_mol": result["final_compartments_mol"]["gas_loop"].get("HCN", 0.0),
        "final_aqueous_HCN_plus_CN_minus_mol": result["hcn_formation"]["final_hcn_plus_cyanide_aqueous_mol"],
        "charged_final_HCN_plus_CN_minus_total_mol": (
            sum(charged_cyanide_by_compartment.values())
            if charged_compartments else None
        ),
        "charged_final_aqueous_HCN_plus_CN_minus_mol": (
            sum(charged_cyanide_by_compartment.get(place, 0.0) for place in (
                "boiler_liquid", "condensate_trap", "collected_condensate",
            )) if charged_compartments else None
        ),
        "HCN_numerical_reporting_floor_mol": result["hcn_formation"]["numerical_reporting_floor_mol"],
        "o1d_first_collision_audit": {
            "evaluated_substeps": len(o1d_fates),
            "status_counts": dict(Counter(item["status"] for item in o1d_fates)),
            "first_resolved_ground_O_fraction": next(
                (item.get("physical_ground_O_fraction") for item in o1d_fates
                 if item["status"] == "conditional_fate_resolved"), None,
            ),
            "maximum_unreviewed_bath_mole_fraction": max(
                (float(item.get("unreviewed_bath_mole_fraction") or 0.0)
                 for item in o1d_fates), default=0.0,
            ),
            "maximum_unreviewed_bath_amount_mol_by_species": {
                species: max(
                    float(item.get("unreviewed_bath_amounts_mol", {}).get(species, 0.0))
                    for item in o1d_fates
                )
                for species in sorted({
                    species for item in o1d_fates
                    for species in item.get("unreviewed_bath_amounts_mol", {})
                })
            },
            "products_injected_into_solver": False,
        },
        "excluded_CO2_excited_fragment_channels": [
            item for item in excluded if item.get("product_label") == "CO(X) + O(1D)"
            and item.get("reason") == "unresolved_excited_state"
        ],
        "accepted_neutral_fragment_extent_mol_by_target": {
            target: sum(
                float(item["accepted_extent_mol"])
                for item in accepted_neutral_sources if item["target"] == target
            )
            for target in ("N2", "CO2", "O2", "H2", "H2O")
        },
        "primary_negative_species_projection_total_absolute_mol": sum(
            float((step.get("solver") or {}).get("conservative_projection_total_absolute_mol") or 0.0)
            for step in heavy_steps
        ),
        "primary_negative_species_projection_substeps": sum(
            bool((step.get("solver") or {}).get("negative_species_clipped_mol"))
            for step in heavy_steps
        ),
        "primary_positivity_retry_substeps": sum(
            bool((step.get("solver") or {}).get("positivity_retry_used"))
            for step in heavy_steps
        ),
        "primary_HCN_projection_adjustment_total_absolute_mol": sum(
            abs(float((step.get("solver") or {}).get(
                "conservative_projection_adjustment_mol_by_species", {}
            ).get("HCN", 0.0)))
            for step in heavy_steps
        ),
        "primary_worst_local_element_residual_mol": max(
            (
                abs(value)
                for step in heavy_steps
                for value in step.get("local_primary_conserved_residual_mol", {}).values()
            ),
            default=0.0,
        ),
        "all_excluded_fragment_reasons": sorted({item.get("reason") for item in excluded}),
        "material_ledger_closed": result["material_ledger"]["closed"],
        "charged_branch_material_ledger_closed": (branch.get("material_ledger") or {}).get("closed"),
        "charged_branch_failed_substeps": branch.get("failed_substep_count"),
        "charged_branch_unprojected_ion_roundoff_clip_mol": charged_roundoff_clip_mol,
        "charged_branch_failure_examples": branch.get("failed_substep_examples"),
        "energy_basis": result["energy_ledger"]["deposited_energy_basis"],
        "thermal_balance_closed": result["energy_ledger"]["complete_thermal_energy_balance_closed"],
        "effective_formaldehyde_pool_without_explicit_hydration": result[
            "scientific_firewall"
        ]["effective_formaldehyde_pool_without_explicit_hydration"],
        "numerical_convergence_validated": False,
        "science_boundary": (
            "Feed fractions are from one published row, but pressure, temperature, spark pulse, "
            "coupling and active zone are hypothetical. HCN measurements are never loaded. "
            "Substep convergence is not established. This cannot be compared with the "
            "full-duration aqueous assay."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--duration-s", type=float, default=60.0)
    parser.add_argument("--spark-off", action="store_true")
    parser.add_argument("--substeps", type=int, default=2)
    parser.add_argument("--cycles", type=int, default=1)
    parser.add_argument("--gas-water-kla-per-s", type=float)
    parser.add_argument("--gas-water-solver", choices=("post_step_analytic", "coupled_bdf"),
                        default="post_step_analytic")
    parser.add_argument("--primary-source-forcing", choices=("time_averaged", "inventory_limited"),
                        default="time_averaged")
    parser.add_argument("--formaldehyde-hydration-model", choices=(
        "off", "winkelman_2002_source_bounded",
        "winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    ), default="off")
    parser.add_argument("--remove-carbon-for-control", action="store_true")
    parser.add_argument("--remove-nitrogen-for-control", action="store_true")
    arguments = parser.parse_args()
    print(json.dumps(run_probe(
        duration_s=arguments.duration_s,
        spark_on=not arguments.spark_off,
        substeps=arguments.substeps,
        cycles=arguments.cycles,
        gas_water_kla_per_s=arguments.gas_water_kla_per_s,
        gas_water_solver=arguments.gas_water_solver,
        primary_source_forcing=arguments.primary_source_forcing,
        formaldehyde_hydration_model=arguments.formaldehyde_hydration_model,
        remove_carbon_for_control=arguments.remove_carbon_for_control,
        remove_nitrogen_for_control=arguments.remove_nitrogen_for_control,
    ), indent=2))


if __name__ == "__main__":
    main()
