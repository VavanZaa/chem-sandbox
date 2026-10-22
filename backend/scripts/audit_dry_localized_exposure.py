"""Conditional two-zone dry-flow plasma sensitivity, never an HCN calibration.

One limiting interpretation of an electrical space-time fraction is a fixed
active gas-volume fraction during each chemistry substep, with complete remix
afterward. Real nanosecond streamers, their duty cycle, repeated gas encounters,
and local density are unmeasured here. This must not be called a reproduction.
"""
from __future__ import annotations

import json
from typing import Any

import numpy as np
from scipy.constants import R

from app.prebiotic.boltzmann import COMMUNITY_CH4_DATASET_ID
from app.prebiotic.coupled_plasma import solve_mixture_plasma_snapshot
from app.prebiotic.heavy_gas import run_gri_heavy_species_cycle
from app.prebiotic.schema import BoltzmannDischargeProgram
from scripts.audit_dry_electrical_exposure import audit_dry_electrical_exposure
from scripts.run_dry_flow_hcn_probe import (
    AMBIENT_FLOW_REFERENCE_K,
    AVERAGE_DISCHARGE_POWER_W,
    DOWNSTREAM_PROBE_TEMPERATURE_K,
    PRESSURE_PA,
    REACTOR_VOLUME_L,
    RESIDENCE_S,
    _gri_inventory,
)


LOCAL_FIELD_TD = 435.0


def run_dry_localized_exposure_probe(
    *,
    active_volume_fraction: float,
    local_electron_density_m3: float,
    methane_fraction: float = 0.01,
    substeps: int = 16,
    discharge_on: bool = True,
) -> dict[str, Any]:
    if not 0 < active_volume_fraction <= 1:
        raise ValueError("Active gas-volume fraction must lie in (0, 1]")
    if local_electron_density_m3 <= 0:
        raise ValueError("Local electron density must be positive")
    if not 0 <= methane_fraction <= 1:
        raise ValueError("Methane fraction must lie in [0, 1]")
    if not 1 <= substeps <= 128:
        raise ValueError("Substeps must lie in [1, 128]")
    temperature_K = DOWNSTREAM_PROBE_TEMPERATURE_K
    packet_residence_s = (
        RESIDENCE_S * AMBIENT_FLOW_REFERENCE_K / temperature_K
    )
    duration_s = packet_residence_s / substeps
    deposited_energy_J = (
        AVERAGE_DISCHARGE_POWER_W * packet_residence_s if discharge_on else 0.0
    )
    energy_per_step_J = deposited_energy_J / substeps
    total_feed_mol = PRESSURE_PA * REACTOR_VOLUME_L * 1e-3 / (R * temperature_K)
    prepared, inlet = _gri_inventory(total_feed_mol, methane_fraction)
    outlet = inlet.copy()
    controls = (
        BoltzmannDischargeProgram.model_validate({
            "collision_set": COMMUNITY_CH4_DATASET_ID,
            "reduced_electric_field_Td": LOCAL_FIELD_TD,
            "electron_number_density": {
                "value": local_electron_density_m3,
                "unit": "1/m^3",
                "source": "conditional local interpretation of paper electrical estimate",
                "method": "declared density ceiling; may be power-capped",
                "basis": "hypothetical",
            },
            "acknowledge_two_term_limit": True,
        })
        if discharge_on else None
    )
    source_consumption_mol_by_target: dict[str, float] = {}
    excluded_fragment_channels: dict[str, str] = {}
    actual_electron_densities_m3: list[float] = []
    allocated_collision_J = 0.0
    fragment_energy_J = 0.0
    active_volume_L = REACTOR_VOLUME_L * active_volume_fraction
    inactive_volume_L = REACTOR_VOLUME_L - active_volume_L
    for _ in range(substeps):
        active = outlet * active_volume_fraction
        inactive = outlet - active
        plasma = (
            solve_mixture_plasma_snapshot(
                controls=controls,
                identities=prepared["identities"],
                gas_amounts_mol=active,
                gas_volume_L=active_volume_L,
                gas_temperature_K=temperature_K,
                discharge_duration_s=duration_s,
                duty_fraction=1.0,
                deposited_energy_J=energy_per_step_J,
            ) if controls is not None else {"status": "discharge_off", "source_terms": []}
        )
        if controls is not None:
            actual_electron_densities_m3.append(float(plasma["electron_number_density_m3"]))
        active, chemistry = run_gri_heavy_species_cycle(
            active, prepared, plasma, temperature_K, active_volume_L,
            duration_s, energy_per_step_J, source_forcing="inventory_limited",
            include_diagnostic_shadows=False,
        )
        allocated_collision_J += float(
            chemistry["collision_energy_allocation"]["allocated_collision_energy_J"]
        )
        fragment_energy_J += float(chemistry["source_energy_used_J"])
        for item in chemistry["source_extents"]:
            target = str(item["target"])
            source_consumption_mol_by_target[target] = (
                source_consumption_mol_by_target.get(target, 0.0)
                + float(item.get("actual_integrated_extent_mol", item["accepted_extent_mol"]))
            )
        for item in chemistry["collision_energy_allocation"]["excluded_fragment_channels"]:
            excluded_fragment_channels[str(item["channel_id"])] = str(item["reason"])
        if inactive_volume_L > 0:
            inactive, _ = run_gri_heavy_species_cycle(
                inactive, prepared, {"status": "discharge_off", "source_terms": []},
                temperature_K, inactive_volume_L, duration_s, 0.0,
                source_forcing="inventory_limited", include_diagnostic_shadows=False,
            )
        outlet = active + inactive
    residual = prepared["atomic"] @ (outlet - inlet)
    tolerance = np.maximum(np.abs(prepared["atomic"] @ inlet) * 1e-12, 1e-12)
    ch4_inlet = float(inlet[prepared["index"]["CH4"]])
    ch4_outlet = float(outlet[prepared["index"]["CH4"]])
    hcn_outlet = float(outlet[prepared["index"]["HCN"]])
    return {
        "status": "conditional_two_zone_sensitivity_unmeasured_encounter_count",
        "scenario": {
            "methane_feed_fraction": methane_fraction,
            "local_field_Td": LOCAL_FIELD_TD,
            "declared_local_electron_density_m3": local_electron_density_m3,
            "active_volume_fraction_each_substep": active_volume_fraction,
            "complete_remix_after_each_substep": True,
            "local_density_is_independently_measured": False,
            "spatial_fraction_is_independently_measured": False,
            "joint_space_time_fraction_is_equated_to_continuous_active_volume": True,
            "individual_streamer_pulses_resolved": False,
            "temperature_K_uniform_proxy": temperature_K,
            "packet_residence_s": packet_residence_s,
            "substeps": substeps,
            "mixing_encounter_count": substeps,
            "mixing_encounter_count_independently_measured": False,
            "varying_substeps_changes_the_physical_mixing_scenario": True,
            "discharge_on": discharge_on,
        },
        "outlet": {
            "ch4_mol": ch4_outlet,
            "ch4_conversion_fraction": (
                (ch4_inlet - ch4_outlet) / ch4_inlet if ch4_inlet > 0 else None
            ),
            "c2h6_mol": float(outlet[prepared["index"]["C2H6"]]),
            "hcn_mol": hcn_outlet,
            "hcn_below_1e-12_mol_reporting_floor": hcn_outlet < 1e-12,
        },
        "source_consumption_mol_by_target": source_consumption_mol_by_target,
        "effective_local_electron_density_m3_range": (
            [min(actual_electron_densities_m3), max(actual_electron_densities_m3)]
            if actual_electron_densities_m3 else None
        ),
        "n2_unresolved_product_state_source_excluded": (
            excluded_fragment_channels.get("N2_dissociation_outcome_000")
            == "unresolved_electronic_product_state"
        ) if discharge_on and methane_fraction < 1.0 else True,
        "ledger": {
            "maximum_material_residual_mol": float(np.max(np.abs(residual), initial=0.0)),
            "residual_mol_by_conserved_quantity": {
                key: float(value)
                for key, value in zip(prepared["conserved"], residual, strict=True)
            },
            "tolerance_mol_by_conserved_quantity": {
                key: float(value)
                for key, value in zip(prepared["conserved"], tolerance, strict=True)
            },
            "material_closed": bool(np.all(np.abs(residual) <= tolerance)),
            "charge_closed": bool(
                abs(residual[prepared["conserved"].index("charge")])
                <= tolerance[prepared["conserved"].index("charge")]
            ),
            "declared_discharge_J": deposited_energy_J,
            "allocated_collision_J": allocated_collision_J,
            "resolved_fragment_J": fragment_energy_J,
            "collision_energy_upper_bound_closed": bool(
                fragment_energy_J <= allocated_collision_J * (1 + 1e-9)
                and allocated_collision_J <= deposited_energy_J * (1 + 1e-9)
            ),
            "full_thermal_balance_closed": False,
        },
        "published_hcn_measurements_used_as_inputs": False,
        "quantitative_hcn_or_conversion_validation_allowed": False,
    }


def audit_dry_localization_bracket() -> dict[str, Any]:
    electrical = audit_dry_electrical_exposure()
    fractions = electrical["conditional_uniform_drift_power"]
    densities = electrical["paper_printed_values"]
    geometric_density = electrical["arithmetic_audit"][
        "electron_density_from_printed_current_and_area_m3"
    ]
    candidates = (
        (
            "published_density_if_local",
            fractions["power_compatible_space_time_fraction_upper_if_published_density_local"],
            densities["electron_density_m3"],
        ),
        (
            "current_area_recomputed_density_if_local",
            fractions["power_compatible_space_time_fraction_upper_if_recomputed_density_local"],
            geometric_density,
        ),
    )
    rows = []
    for label, upper_fraction, density in candidates:
        prediction = run_dry_localized_exposure_probe(
            active_volume_fraction=upper_fraction,
            local_electron_density_m3=density,
        )
        rows.append({
            "electrical_interpretation": label,
            "conditional_power_compatible_joint_fraction_upper": upper_fraction,
            "fraction_used_as_active_volume_sensitivity": upper_fraction,
            "prediction": prediction,
        })
    return {
        "status": "conditional_localization_bracket_unmeasured_exposure_and_encounters",
        "rows": rows,
        "hcn_assay_loaded": False,
        "product_fitted_streamer_fraction_used": False,
        "independent_occupancy_measured": False,
        "mixing_encounter_count_independently_measured": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_localization_bracket(), indent=2))
