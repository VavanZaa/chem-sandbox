"""Audit one declared apparatus hour without fitting to historical products.

The first modeled discharge is traced by integrating signed reaction progress,
not by extrapolating the final instantaneous rate. This is a mechanism diagnostic,
not a measured yield or a change to the blind-physics solver.
"""

from __future__ import annotations

import json
import argparse
import copy

import cantera as ct
import numpy as np
from scipy.integrate import solve_ivp

from app.prebiotic.apparatus import _prepare, simulate_prebiotic_apparatus
from app.prebiotic.boltzmann import COMMUNITY_CH4_DATASET_ID, REACTIVE_N2_DATASET_ID
from app.prebiotic.heavy_gas import _neutral_fragment_extents, augment_prepared_with_gri
from app.prebiotic.precursor_benchmark import stribling_miller_measurement_profile
from app.prebiotic.schema import PrebioticApparatusRequest


def pulse_sensitivity(
    request: PrebioticApparatusRequest, snapshot: dict, deposited_energy_J: float
) -> dict:
    """Counterfactual pulse spacing in one well-mixed isothermal gas box.

    This omits condenser capture and charged-species feedback. Event counts are
    numerical timing scenarios, not a claim about parcel-resolved spark exposure.
    """
    prepared = _prepare(request)
    augment_prepared_with_gri(prepared)
    gas = np.zeros(len(prepared["ids"]), dtype=float)
    zone_fraction = float(snapshot["discharge_zone_fraction"])
    for item in snapshot["source_terms"]:
        target = item["target"]
        if target in prepared["index"]:
            gas[prepared["index"][target]] = (
                float(item["target_amount_mol"]) / zone_fraction
            )
    injected, _, _, _, _ = _neutral_fragment_extents(
        gas, prepared, snapshot, deposited_energy_J
    )
    phase = ct.Solution("gri30.yaml")
    positions = np.asarray([prepared["index"][name] for name in phase.species_names])
    initial = gas[positions]
    fragment_delta = injected[positions] - initial
    volume_m3 = request.apparatus.gas_loop_volume.as_unit("L") * 1e-3
    temperature_K = request.apparatus.discharge_temperature.as_unit("K")
    molecular_weights = np.asarray(phase.molecular_weights)
    duration_s = float(snapshot["whole_loop_mean_residence_time_s"])

    def rate(_time_s: float, amounts: np.ndarray) -> np.ndarray:
        bounded = np.maximum(amounts, 0.0)
        mass_kg = float(np.dot(bounded, molecular_weights) / 1000.0)
        phase.TDX = temperature_K, mass_kg / volume_m3, bounded
        return np.asarray(phase.net_production_rates) * 1000.0 * volume_m3

    def run_pulses(fragment_change: np.ndarray, pulses: int) -> dict[str, float | int]:
        amounts = initial.copy()
        interval_s = duration_s / pulses
        evaluation_count = 0
        for _ in range(pulses):
            amounts += fragment_change / pulses
            integration = solve_ivp(
                rate,
                (0.0, interval_s),
                amounts,
                method="BDF",
                rtol=1e-8,
                atol=1e-18,
                max_step=max(interval_s / 20, 1e-9),
            )
            if not integration.success:
                raise RuntimeError(integration.message)
            amounts = np.maximum(integration.y[:, -1], 0.0)
            evaluation_count += integration.nfev
        return {
            "HCN_mol": float(amounts[phase.species_index("HCN")]),
            "H2CN_mol": float(amounts[phase.species_index("H2CN")]),
            "NO_mol": float(amounts[phase.species_index("NO")]),
            "N2_mol": float(amounts[phase.species_index("N2")]),
            "function_evaluations": evaluation_count,
        }
    cases = {str(pulses): run_pulses(fragment_delta, pulses) for pulses in (1, 6, 60)}
    ablations = {}
    for label in ("water_atomic_O_fragments_disabled", "all_water_fragments_disabled"):
        variant = copy.deepcopy(snapshot)
        for source in variant["source_terms"]:
            if source["target"] != "H2O":
                continue
            if label == "all_water_fragments_disabled" or any(
                product["formula"] == "O" for product in source.get("products", [])
            ):
                source["reactive"] = False
        variant_injected, _, _, _, _ = _neutral_fragment_extents(
            gas, prepared, variant, deposited_energy_J
        )
        ablations[label] = run_pulses(variant_injected[positions] - initial, 1)
    return {
        "status": "counterfactual_pulse_timing_sensitivity_not_validated",
        "temperature_K": temperature_K,
        "gas_volume_L": volume_m3 * 1000,
        "modeled_chemical_duration_s": duration_s,
        "deposited_energy_J": deposited_energy_J,
        "condenser_transport_included": False,
        "electron_collisions_recomputed_between_pulses": False,
        "event_count_is_numerical_not_measured_parcel_exposure": True,
        "cases": cases,
        "source_ablations": ablations,
        "ablations_are_diagnostic_only_not_physical_controls": True,
    }


def energy_sensitivity(profile: dict, baseline: dict) -> dict:
    """One-hour current bracket, not a self-consistent circuit calculation."""
    cases = {}
    for factor in (0.1, 1.0, 10.0):
        if factor == 1.0:
            result = baseline
        else:
            variant = copy.deepcopy(profile)
            variant["spark"]["current"]["value"] *= factor
            result = simulate_prebiotic_apparatus(
                PrebioticApparatusRequest.model_validate(variant)
            )
        energy_J = result["energy_ledger"]["deposited_energy_used_J"]
        collection = result["final_compartments_mol"]["collected_condensate"]
        first_plasma = result["plasma_model"]["cycle_summaries"][0]
        cases[str(factor)] = {
            "deposited_energy_J": energy_J,
            "collected_HCN_mol": collection.get("HCN", 0.0),
            "collected_CH2O_mol": collection.get("CH2O", 0.0),
            "one_hour_collected_HCN_nmol_per_J": collection.get("HCN", 0.0) * 1e9 / energy_J,
            "one_hour_collected_CH2O_nmol_per_J": collection.get("CH2O", 0.0) * 1e9 / energy_J,
            "material_ledger_closed": result["material_ledger"]["closed"],
            "charged_ledger_closed": result["plasma_model"]["propagated_charged_branch"]
                ["material_ledger"]["closed"],
            "first_step_effective_electron_density_m3": first_plasma[
                "electron_number_density_m3"
            ],
            "first_step_drift_to_deposited_energy_ratio_before_cap": first_plasma[
                "drift_power_diagnostic"
            ]["drift_to_deposited_energy_ratio"],
        }
    return {
        "status": "hypothetical_one_hour_energy_sensitivity_not_historical_calorimetry",
        "all_other_declared_inputs_held_fixed": True,
        "historical_yield_used_to_choose_factors": False,
        "not_a_zero_time_energy_yield": True,
        "electron_density_energy_capped_per_case": True,
        "electron_density_recomputed_from_circuit": False,
        "circuit_and_electrode_physics_recomputed": False,
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pulse-sweep", action="store_true")
    parser.add_argument("--energy-sweep", action="store_true")
    parser.add_argument(
        "--collision-set",
        choices=(REACTIVE_N2_DATASET_ID, COMMUNITY_CH4_DATASET_ID),
        default=REACTIVE_N2_DATASET_ID,
    )
    arguments = parser.parse_args()
    profile = stribling_miller_measurement_profile()
    profile["spark"]["boltzmann"]["collision_set"] = arguments.collision_set
    profile["apparatus"].update(
        cycles=1,
        integration_substeps_per_cycle=2,
        primary_source_forcing="time_averaged",
        condenser_phase_transfer_model="coupled_flow_transport",
        condensate_nonwater_return_fraction=0,
        condenser_contact_water_holdup={
            "value": 0.1,
            "unit": "L",
            "source": "hypothetical contact-water sensitivity; not an apparatus measurement",
            "method": "declared scenario parameter",
            "basis": "hypothetical",
        },
    )
    request = PrebioticApparatusRequest.model_validate(profile)
    result = simulate_prebiotic_apparatus(request)
    first = result["plasma_model"]["heavy_species_cycle_summaries"][0]
    audit = first["integrated_pathway_diagnostics"]
    if audit is None:
        raise RuntimeError("The first-step integrated pathway audit is unavailable")
    output = {
        "scenario": "declared one-hour Stribling/Miller-family sensitivity, not exact 1986 apparatus",
        "collision_set": arguments.collision_set,
        "historical_product_data_used_by_solver": False,
        "hypothetical_operating_inputs": result["evidence"]["hypothetical_sources"],
        "first_step_duration_s": first["duration_s"],
        "first_step_source_energy_J": first["collision_energy_allocation"],
        "first_step_uniform_drift_power_audit": result["plasma_model"]["cycle_summaries"][0]
            ["drift_power_diagnostic"],
        "first_step_source_extents_mol": {
            item["channel_id"]: item["accepted_extent_mol"]
            for item in first["source_extents"]
        },
        "integrated_pathways": audit,
        "material_ledger_closed": result["material_ledger"]["closed"],
        "first_hour_collected_mol": {
            species: result["final_compartments_mol"]["collected_condensate"].get(species, 0)
            for species in ("HCN", "CH2O")
        },
    }
    if arguments.pulse_sweep:
        output["pulse_timing_sensitivity"] = pulse_sensitivity(
            request,
            result["plasma_model"]["cycle_summaries"][0],
            result["energy_ledger"]["deposited_energy_used_J"]
            / len(result["plasma_model"]["cycle_summaries"]),
        )
    if arguments.energy_sweep:
        output["electrical_energy_sensitivity"] = energy_sensitivity(profile, result)
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
