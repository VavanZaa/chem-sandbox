from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

import cantera as ct
import numpy as np
import pytest
from fastapi.testclient import TestClient
from scipy.constants import electron_mass, elementary_charge

from app.main import app
from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.blind_aqueous import (
    discover_aqueous_candidates,
    enumerate_neutral_tree_isomers,
)
from app.prebiotic.boltzmann import (
    HYBRID_DATASET_ID,
    IST_CO2_O2_DATASET_ID,
    REACTIVE_N2_DATASET_ID,
    _tabata_n2_analytic_fit_m2,
    collision_dataset_status,
    solve_self_consistent_eedf,
)
from app.prebiotic.coupled_plasma import _uniform_drift_power_diagnostic
from app.prebiotic.heavy_gas import _conservative_clip_pavan_ions
from app.prebiotic.candidate_ranking import (
    build_paired_starting_geometries,
    CandidateRankingJobManager,
    evaluate_frozen_candidate_ranking,
    select_amino_carboxylic_candidates,
)
from app.prebiotic.capabilities import prebiotic_capabilities
from app.prebiotic.plasma import (
    calculate_electron_collision_rate,
    calculate_rate_coefficient,
    validate_electron_collision_dataset,
)
from app.prebiotic.phase_partition import (
    capture_gas_water_by_flow,
    equilibrate_gas_water,
    enforce_water_vapor_capacity,
)
from app.prebiotic.historical_benchmark import (
    evaluate_historical_benchmark,
    historical_material_time_profile,
)
from app.prebiotic.post_discovery import (
    evaluate_frozen_aqueous_discovery,
    freeze_discovery_artifact,
)
from app.prebiotic.precursor_benchmark import (
    evaluate_precursor_energy_benchmark,
    stribling_miller_measurement_profile,
)
from app.prebiotic.schema import (
    BoltzmannEEDFRequest,
    ElectronCollisionDatasetRequest,
    ElectronCollisionRateRequest,
    PrebioticApparatusRequest,
)


def evidence(value: float, unit: str, source: str = "independent apparatus fixture") -> dict:
    return {
        "value": value,
        "unit": unit,
        "source": source,
        "method": "declared validation fixture",
        "basis": "hypothetical",
    }


def test_tiny_negative_pavan_ion_projection_conserves_atoms_and_charge() -> None:
    packaged = next(
        Path(directory) / "example_data/methane-plasma-pavan-2023.yaml"
        for directory in ct.get_data_directories()
        if (Path(directory) / "example_data/methane-plasma-pavan-2023.yaml").is_file()
    )
    phase = ct.Solution(str(packaged), "gri30-plasma")
    raw = np.zeros(phase.n_species)
    raw[phase.species_index("N2")] = 0.02
    raw[phase.species_index("e")] = 1e-9
    raw[phase.species_index("N2+")] = -1e-12
    final, clipped, projected = _conservative_clip_pavan_ions(raw, phase)
    atoms = np.asarray([
        [phase.n_atoms(position, element) for position in range(phase.n_species)]
        for element in phase.element_names
    ])
    charge = np.asarray([species.charge for species in phase.species()])
    assert clipped == {"N2+": pytest.approx(1e-12)}
    assert projected == {"N2+": pytest.approx(1e-12)}
    assert final[phase.species_index("N2+")] == 0
    assert np.max(np.abs(atoms @ (final - raw))) < 1e-17
    assert abs(charge @ (final - raw)) < 1e-17
    raw[phase.species_index("N2+")] = -3e-11
    with pytest.raises(ValueError, match="projection bound"):
        _conservative_clip_pavan_ions(raw, phase)


def test_unrepresentable_pavan_ion_roundoff_is_bounded_and_visible() -> None:
    packaged = next(
        Path(directory) / "example_data/methane-plasma-pavan-2023.yaml"
        for directory in ct.get_data_directories()
        if (Path(directory) / "example_data/methane-plasma-pavan-2023.yaml").is_file()
    )
    phase = ct.Solution(str(packaged), "gri30-plasma")
    raw = np.zeros(phase.n_species)
    raw[phase.species_index("N2")] = 0.02
    raw[phase.species_index("e")] = 1e-9
    raw[phase.species_index("CH4")] = 1e-23
    raw[phase.species_index("CH4+")] = -4e-21
    final, clipped, projected = _conservative_clip_pavan_ions(raw, phase)
    atoms = np.asarray([
        [phase.n_atoms(position, element) for position in range(phase.n_species)]
        for element in phase.element_names
    ])
    charge = np.asarray([species.charge for species in phase.species()])
    assert clipped["CH4+"] == pytest.approx(4e-21)
    assert "CH4+" not in projected
    assert final[phase.species_index("CH4+")] == 0
    assert np.max(np.abs(atoms @ (final - raw))) <= 2e-20
    assert abs(charge @ (final - raw)) <= 4e-21
    raw[phase.species_index("CH4+")] = -1.7e-18
    final, clipped, projected = _conservative_clip_pavan_ions(raw, phase)
    assert clipped["CH4+"] == pytest.approx(1.7e-18)
    assert "CH4+" not in projected
    assert np.max(np.abs(atoms @ (final - raw))) <= 7e-18
    assert abs(charge @ (final - raw)) <= 1.7e-18
    raw[phase.species_index("CH4+")] = -1.1e-17
    with pytest.raises(ValueError, match="projection"):
        _conservative_clip_pavan_ions(raw, phase)


def test_trace_pavan_ion_uses_conservative_bath_projection_when_parent_missing() -> None:
    packaged = next(
        Path(directory) / "example_data/methane-plasma-pavan-2023.yaml"
        for directory in ct.get_data_directories()
        if (Path(directory) / "example_data/methane-plasma-pavan-2023.yaml").is_file()
    )
    phase = ct.Solution(str(packaged), "gri30-plasma")
    raw = np.zeros(phase.n_species)
    for species, amount in {
        "N2": 0.02, "CO2": 0.002, "H2": 0.002, "H2O": 1e-4,
        "O2": 1e-4, "e": 1e-9, "CH4": 6e-21, "CH4+": -1.865e-17,
    }.items():
        raw[phase.species_index(species)] = amount
    final, clipped, projected = _conservative_clip_pavan_ions(raw, phase)
    conserved = np.asarray([
        [phase.n_atoms(position, element) for position in range(phase.n_species)]
        for element in phase.element_names
    ] + [[species.charge for species in phase.species()]])
    assert clipped["CH4+"] == pytest.approx(1.865e-17)
    assert projected["CH4+"] == pytest.approx(1.865e-17)
    assert np.min(final) >= 0
    assert np.max(np.abs(conserved @ (final - raw))) < 5e-18
    raw[phase.species_index("CH4+")] = -1.1e-15
    with pytest.raises(ValueError, match="projection bound"):
        _conservative_clip_pavan_ions(raw, phase)


def apparatus_payload() -> dict:
    return {
        "id": "closed_prebiotic_apparatus",
        "species": [
            {
                "id": "H2O",
                "name": "Water",
                "composition": {"H": 2, "O": 1},
                "charge": 0,
                "smiles": "O",
                "phase_behavior": "condensable",
                "evaporation_rate": evidence(0.2, "1/s"),
                "condensation_rate": evidence(0.4, "1/s"),
            },
            {
                "id": "CH4",
                "name": "Methane",
                "composition": {"C": 1, "H": 4},
                "charge": 0,
                "smiles": "C",
                "phase_behavior": "permanent_gas",
            },
        ],
        "reactions": [],
        "initial_inventory": [
            {
                "species_id": "H2O",
                "compartment": "boiler_liquid",
                "amount": {"value": 1, "unit": "mol"},
            },
            {
                "species_id": "CH4",
                "compartment": "gas_loop",
                "amount": {"value": 0.1, "unit": "mol"},
            },
        ],
        "apparatus": {
            "cycles": 3,
            "boiler_liquid_volume": {"value": 1, "unit": "L"},
            "gas_loop_volume": {"value": 2, "unit": "L"},
            "condensate_trap_volume": {"value": 0.2, "unit": "L"},
            "gas_circulation_rate": {"value": 1, "unit": "L/s"},
            "pressure": {"value": 1, "unit": "atm"},
            "boiler_temperature": {"value": 373.15, "unit": "K"},
            "discharge_temperature": {"value": 400, "unit": "K"},
            "condenser_temperature": {"value": 280, "unit": "K"},
            "boil_duration": {"value": 4, "unit": "s"},
            "discharge_duration": {"value": 2, "unit": "s"},
            "condense_duration": {"value": 3, "unit": "s"},
            "condensate_return_fraction": 1,
        },
        "spark": {"mode": "off"},
        "interventions": [],
        "mechanism_origin": "user_declared",
        "knowledge_policy": "blind_physics_only",
        "acknowledge_reduced_order_limit": True,
    }


def spark() -> dict:
    return {
        "mode": "reduced_order",
        "voltage": {"value": 20_000, "unit": "V"},
        "current": {"value": 0.1, "unit": "A"},
        "pulse_width": {"value": 0.001, "unit": "s"},
        "frequency": {"value": 100, "unit": "Hz"},
        "coupling_efficiency": 0.5,
        "mean_electron_energy": {"value": 5, "unit": "eV"},
        "source": "hypothetical electrical validation fixture",
        "basis": "hypothetical",
    }


def boltzmann_spark() -> dict:
    return {
        "mode": "boltzmann_eedf",
        "voltage": {"value": 20_000, "unit": "V"},
        "current": {"value": 0.01, "unit": "A"},
        "pulse_width": {"value": 0.001, "unit": "s"},
        "frequency": {"value": 100, "unit": "Hz"},
        "coupling_efficiency": 0.3,
        "boltzmann": {
            "solver": "bolos_two_term",
            "collision_set": "lxcat_morgan_2026-09-22",
            "reduced_electric_field_Td": 200,
            "electron_number_density": evidence(1e15, "1/m^3"),
            "energy_cells": 120,
            "maximum_iterations": 200,
            "acknowledge_two_term_limit": True,
        },
        "source": "hypothetical coupled-plasma validation fixture",
        "basis": "hypothetical",
    }


def reaction_payload() -> dict:
    payload = apparatus_payload()
    payload["species"] = [
        {
            "id": "A",
            "name": "Carbon-bearing precursor A",
            "composition": {"C": 1},
            "charge": 0,
            "phase_behavior": "permanent_gas",
        },
        {
            "id": "B",
            "name": "Carbon-bearing product B",
            "composition": {"C": 1},
            "charge": 0,
            "phase_behavior": "nonvolatile",
            "condensation_rate": evidence(1, "1/s"),
        },
    ]
    payload["reactions"] = [
        {
            "id": "electron_channel",
            "name": "Abstract balanced electron-impact channel",
            "phase": "gas",
            "stoichiometry": {"A": -1, "B": 1},
            "rate_law": {
                "kind": "electron_impact_reduced",
                "rate_constant": evidence(1, "1/s"),
                "energy_cost": evidence(1000, "J/mol"),
            },
        }
    ]
    payload["initial_inventory"] = [
        {
            "species_id": "A",
            "compartment": "gas_loop",
            "amount": {"value": 0.1, "unit": "mol"},
        }
    ]
    return payload


def test_closed_recirculating_apparatus_conserves_material() -> None:
    request = PrebioticApparatusRequest.model_validate(apparatus_payload())
    result = simulate_prebiotic_apparatus(request)
    assert result["material_ledger"]["closed"] is True
    assert result["cycles_completed"] == 3
    assert result["cycles"][0]["evaporated_mol"]["H2O"] > 0
    assert result["cycles"][0]["condensed_mol"]["H2O"] > 0
    assert result["final_total_mol"]["H2O"] == pytest.approx(1)
    assert result["scientific_firewall"]["known_products_accessible_to_discovery"] is False
    assert result["aqueous_cyanide_speciation"]["prediction_inventory_mutated"] is False


def test_apparatus_inventory_drives_cycle_resolved_boltzmann_source_terms() -> None:
    payload = apparatus_payload()
    payload["apparatus"]["cycles"] = 1
    payload["apparatus"]["discharge_zone_volume"] = evidence(0.02, "L")
    payload["spark"] = boltzmann_spark()
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    snapshot = result["cycles"][0]["plasma_snapshot"]
    assert snapshot["status"] == "apparatus_mixture_eedf_solved"
    assert sum(snapshot["mole_fractions"].values()) == pytest.approx(1)
    assert set(snapshot["mole_fractions"]) == {"H2O", "CH4"}
    assert snapshot["mean_electron_energy_eV"] > 0
    drift = snapshot["drift_power_diagnostic"]
    assert drift["estimated_electron_drift_energy_J"] > 0
    assert drift["electron_density_if_all_deposit_went_to_drift_m3"] > 0
    assert drift["used_to_scale_source_terms"] is True
    assert snapshot["electron_number_density_m3"] == pytest.approx(
        drift["electron_density_if_all_deposit_went_to_drift_m3"]
    )
    assert snapshot["electron_number_density_m3"] < snapshot[
        "declared_electron_number_density_m3"
    ]
    assert snapshot["reactive_channel_count"] > 0
    assert snapshot["balanced_reactive_channel_count"] > 0
    assert any(item["role"] == "dissociation" for item in snapshot["source_terms"])
    assert result["plasma_model"]["self_consistent_eedf_solved"] is True
    assert result["plasma_model"]["apparatus_mixture_coupled"] is True
    assert snapshot["discharge_zone_fraction"] == pytest.approx(0.01)
    assert snapshot["whole_loop_mean_residence_time_s"] == pytest.approx(0.02)
    assert result["cycles"][0]["heavy_species_snapshot"]["duration_s"] == pytest.approx(0.02)
    assert result["operating_conditions"][
        "effective_heavy_chemistry_residence_time_s_per_cycle"
    ] == pytest.approx(0.02)
    assert result["plasma_model"]["heavy_species_inventory_coupled"] is True
    assert result["scientific_firewall"]["known_products_accessible_to_discovery"] is True
    assert result["scientific_firewall"]["reaction_recipe_used_for_discovery"] is True
    assert result["scientific_firewall"]["preexisting_gri_hcn_reactions_accessible_to_solver"] is True
    assert result["scientific_firewall"]["blind_discovery_eligible"] is False
    heavy = result["cycles"][0]["heavy_species_snapshot"]
    assert heavy["status"] == "exploratory_neutral_gas_network_integrated"
    assert heavy["formed_mol"]["CH2O"] > 0
    assert heavy["source_energy_used_J"] <= result["cycles"][0]["deposited_electrical_energy_J"]
    allocation = heavy["collision_energy_allocation"]
    assert allocation["all_collision_raw_threshold_demand_J"] >= allocation["resolved_fragment_raw_threshold_demand_J"]
    assert allocation["resolved_fragment_energy_used_J"] <= allocation["allocated_collision_energy_J"]
    expected_scale = min(
        1.0,
        result["cycles"][0]["deposited_energy_used_J"]
        / allocation["all_collision_raw_threshold_demand_J"],
    )
    assert allocation["global_collision_energy_scale"] == pytest.approx(expected_scale)
    comparison = heavy["mechanism_comparison"]
    assert comparison["status"] == "multi_mechanism_shadow_complete"
    assert comparison["shadow"]["mechanism_species"] == 71
    assert comparison["shadow"]["retained_heavy_reactions"] == 387
    assert comparison["shadow"]["excluded_plasma_rate_reactions"] == 53
    assert len(comparison["shadow"]["mechanism_sha256"]) == 64
    assert comparison["scientific_firewall"]["primary_inventory_mutated_by_shadow"] is False
    assert comparison["scientific_firewall"]["comparison_is_independent_confidence"] is False
    literature = comparison["literature_hcn_shadow"]
    assert literature["status"] == "kamarinopoulou_2024_neutral_shadow_complete"
    assert literature["published_species"] == 83
    assert literature["retained_heavy_reactions"] == 189
    assert literature["fixed_condition_electron_reactions_excluded"] == 62
    assert literature["integrated_duration_s"] == pytest.approx(0.018)
    assert literature["residence_time_domain_limited"] is True
    assert len(literature["mechanism_sha256"]) == 64
    assert literature["scientific_firewall"]["paper_experimental_yields_accessible_to_solver"] is False
    assert literature["scientific_firewall"]["primary_inventory_mutated_by_shadow"] is False
    assert literature["scientific_firewall"]["hnc_ground_state_thermo_independently_validated"] is False
    assert any("HNC thermo record conflicts" in item for item in literature["warnings"])
    assert result["material_ledger"]["closed"] is True


def test_uniform_drift_power_audit_scales_with_electron_density() -> None:
    parameters = {
        "reduced_mobility": 1e24,
        "reduced_field_Td": 200.0,
        "neutral_amount_mol": 0.1,
        "volume_m3": 0.01,
        "active_duration_s": 0.1,
        "deposited_energy_J": 1.0,
    }
    first = _uniform_drift_power_diagnostic(
        **parameters, electron_density_m3=1e15
    )
    second = _uniform_drift_power_diagnostic(
        **parameters, electron_density_m3=2e15
    )
    assert second["estimated_electron_drift_energy_J"] == pytest.approx(
        2 * first["estimated_electron_drift_energy_J"]
    )
    assert first["electron_density_if_all_deposit_went_to_drift_m3"] == pytest.approx(
        second["electron_density_if_all_deposit_went_to_drift_m3"]
    )
    assert first["electron_mobility_m2_V_s"] * first["neutral_number_density_m3"] == pytest.approx(
        parameters["reduced_mobility"]
    )


def test_charge_balanced_pavan_branch_propagates_its_own_closed_inventory() -> None:
    payload = apparatus_payload()
    payload["apparatus"]["cycles"] = 2
    payload["apparatus"]["discharge_zone_volume"] = evidence(0.02, "L")
    payload["spark"] = boltzmann_spark()
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    branch = result["plasma_model"]["propagated_charged_branch"]
    assert result["plasma_model"]["propagated_charged_branch_coupled"] is True
    assert branch["status"] == "propagated_charge_balanced_pavan_branch_complete"
    assert branch["cycles_completed"] == 2
    assert branch["material_ledger"]["closed"] is True
    assert max(abs(value) for value in branch["material_ledger"]["residual"].values()) < 1e-10
    first = result["cycles"][0]["heavy_species_snapshot"]["mechanism_comparison"]["charged_shadow"]
    second = result["cycles"][1]["heavy_species_snapshot"]["mechanism_comparison"]["charged_shadow"]
    assert first["carried_from_previous_cycle"] is False
    assert second["carried_from_previous_cycle"] is True
    assert first["branch_plasma_snapshot"]["status"] == "apparatus_mixture_eedf_solved"
    assert second["branch_plasma_snapshot"]["status"] == "apparatus_mixture_eedf_solved"
    assert first["scientific_firewall"]["primary_inventory_mutated_by_shadow"] is False
    assert branch["scientific_firewall"]["quantitative_yield_claim_allowed"] is False


def test_continuous_overlap_split_steps_preserve_declared_time_energy_and_closure() -> None:
    payload = apparatus_payload()
    payload["apparatus"]["cycles"] = 2
    payload["apparatus"]["operation_mode"] = "continuous_overlap"
    payload["apparatus"]["integration_substeps_per_cycle"] = 2
    payload["apparatus"]["discharge_zone_volume"] = evidence(0.02, "L")
    payload["spark"] = boltzmann_spark()
    payload["spark"]["calorimetric_deposited_energy_per_cycle"] = evidence(10.0, "J")
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert result["cycles_completed"] == 2
    assert len(result["cycles"]) == 2
    assert len(result["plasma_model"]["cycle_summaries"]) == 4
    assert len(result["plasma_model"]["heavy_species_cycle_summaries"]) == 4
    assert result["total_duration_s"] == pytest.approx(18.0)
    assert result["cycles"][0]["integration_substeps"] == 2
    assert result["cycles"][0]["end_time_s"] == pytest.approx(9.0)
    assert result["cycles"][1]["end_time_s"] == pytest.approx(18.0)
    assert result["energy_ledger"]["deposited_electrical_energy_J"] == pytest.approx(24.0)
    assert result["energy_ledger"]["deposited_energy_used_J"] == pytest.approx(20.0)
    assert result["cycles"][0]["deposited_energy_used_J"] == pytest.approx(10.0)
    assert result["material_ledger"]["closed"] is True
    assert result["plasma_model"]["propagated_charged_branch"]["material_ledger"]["closed"] is True
    assert result["plasma_model"]["propagated_charged_branch"]["cycles_completed"] == 2
    assert result["operating_conditions"]["operation_mode"] == "continuous_overlap"
    assert result["operating_conditions"]["integration_substeps_per_cycle"] == 2
    assert result["operating_conditions"]["operator_splitting_convergence_assessed_in_this_run"] is False
    assert result["operating_conditions"]["effective_heavy_chemistry_residence_time_s_per_cycle"] == pytest.approx(0.02)
    assert result["operating_conditions"]["effective_heavy_chemistry_residence_time_s_per_integration_step"] == pytest.approx(0.01)


def test_time_averaged_primary_source_preserves_material_and_energy() -> None:
    payload = apparatus_payload()
    payload["apparatus"].update(
        cycles=2,
        operation_mode="continuous_overlap",
        integration_substeps_per_cycle=2,
        primary_source_forcing="time_averaged",
        discharge_zone_volume=evidence(0.02, "L"),
    )
    payload["spark"] = boltzmann_spark()
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    heavy_steps = result["plasma_model"]["heavy_species_cycle_summaries"]
    assert len(heavy_steps) == 4
    assert all(step["source_forcing"] == "time_averaged" for step in heavy_steps)
    assert [step["diagnostic_shadows_evaluated"] for step in heavy_steps] == [
        False, True, False, True
    ]
    assert all(step["solver"]["function_evaluations"] > 0 for step in heavy_steps)
    assert result["material_ledger"]["closed"] is True
    assert result["plasma_model"]["propagated_charged_branch"]["material_ledger"]["closed"] is True
    assert result["energy_ledger"]["deposited_energy_used_J"] == pytest.approx(24.0)
    assert result["operating_conditions"]["primary_source_forcing"] == "time_averaged"
    repeated = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    for species in ("HCN", "CH2O", "CH3OH"):
        assert repeated["final_compartments_mol"]["boiler_liquid"][species] == pytest.approx(
            result["final_compartments_mol"]["boiler_liquid"][species], rel=1e-6, abs=1e-12
        )


def test_cold_trap_retains_all_nonwater_without_species_targeting() -> None:
    payload = apparatus_payload()
    payload["apparatus"].update(cycles=2, condensate_nonwater_return_fraction=0)
    payload["species"].append({
        "id": "CH3OH", "name": "Methanol", "composition": {"C": 1, "H": 4, "O": 1},
        "charge": 0, "smiles": "CO", "phase_behavior": "condensable",
    })
    payload["initial_inventory"].extend([
        {"species_id": "CH3OH", "compartment": "condensate_trap", "amount": {"value": 0.02, "unit": "mol"}},
        {"species_id": "H2O", "compartment": "condensate_trap", "amount": {"value": 0.03, "unit": "mol"}},
    ])
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    final = result["final_compartments_mol"]
    assert final["condensate_trap"]["CH3OH"] == 0
    assert final["collected_condensate"]["CH3OH"] == pytest.approx(0.02)
    assert final["boiler_liquid"]["CH3OH"] == 0
    assert final["gas_loop"]["CH3OH"] == 0
    assert final["condensate_trap"]["H2O"] == 0
    assert result["material_ledger"]["closed"] is True
    assert result["operating_conditions"]["condensate_nonwater_return_fraction"] == 0
    assert any("generic separation scenario" in warning for warning in result["warnings"])
    payload["apparatus"].pop("condensate_nonwater_return_fraction")
    returned = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert returned["final_compartments_mol"]["condensate_trap"]["CH3OH"] == 0
    assert returned["final_compartments_mol"]["collected_condensate"]["CH3OH"] == 0
    assert returned["final_compartments_mol"]["boiler_liquid"]["CH3OH"] == pytest.approx(0.02)
    assert returned["material_ledger"]["closed"] is True


def test_flow_limited_cold_trap_isolation_is_resolution_stable_without_spark() -> None:
    payload = apparatus_payload()
    payload["species"].append({
        "id": "HCN", "name": "Hydrogen cyanide", "composition": {"H": 1, "C": 1, "N": 1},
        "charge": 0, "smiles": "C#N", "phase_behavior": "volatile",
    })
    payload["initial_inventory"].append({
        "species_id": "HCN", "compartment": "gas_loop", "amount": {"value": 1e-5, "unit": "mol"},
    })
    payload["apparatus"].update(
        cycles=2,
        operation_mode="continuous_overlap",
        condensate_nonwater_return_fraction=0,
        condenser_phase_transfer_model="flow_limited_capture",
        condenser_contact_water_holdup=evidence(0.1, "L"),
    )
    collected = []
    for substeps in (2, 8):
        payload["apparatus"]["integration_substeps_per_cycle"] = substeps
        result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
        assert result["material_ledger"]["closed"] is True
        assert result["final_compartments_mol"]["gas_loop"]["HCN"] > 0
        assert result["final_compartments_mol"]["boiler_liquid"]["HCN"] == 0
        collected.append(result["final_compartments_mol"]["collected_condensate"]["HCN"])
    assert collected[0] > 0
    assert collected[1] == pytest.approx(collected[0], rel=1e-3)
    payload["apparatus"].pop("condenser_contact_water_holdup")
    with pytest.raises(ValueError, match="requires a declared contact-water holdup"):
        PrebioticApparatusRequest.model_validate(payload)


def test_coupled_flow_transport_preserves_both_ledgers_and_one_hour_resolution() -> None:
    profile = stribling_miller_measurement_profile()
    profile["apparatus"].update(
        cycles=1,
        primary_source_forcing="time_averaged",
        condensate_nonwater_return_fraction=0,
        condenser_phase_transfer_model="coupled_flow_transport",
        condenser_contact_water_holdup=evidence(0.1, "L"),
    )
    amounts = []
    for substeps in (2, 8):
        profile["apparatus"]["integration_substeps_per_cycle"] = substeps
        result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(profile))
        charged = result["plasma_model"]["propagated_charged_branch"]
        assert result["material_ledger"]["closed"] is True
        assert charged["material_ledger"]["closed"] is True
        assert charged["material_ledger"]["maximum_absolute_residual_mol"] < 1e-11
        assert charged["failed_substep_count"] == 0
        assert charged["scientific_firewall"]["known_products_accessible_to_solver"] is True
        assert result["plasma_model"]["heavy_species_cycle_summaries"][0][
            "scientific_firewall"
        ]["known_products_accessible_to_solver"] is True
        assert result["total_duration_s"] == pytest.approx(3600)
        assert result["energy_ledger"]["deposited_energy_used_J"] == pytest.approx(2160)
        assert result["plasma_model"]["heavy_species_cycle_summaries"][0][
            "coupled_henry_transport"
        ]["status"] == "mean_field_henry_exchange_coupled_to_gri_bdf"
        pathway_audit = result["plasma_model"]["heavy_species_cycle_summaries"][0][
            "integrated_pathway_diagnostics"
        ]
        assert pathway_audit["adaptive_interval_count"] > 20
        assert pathway_audit["not_a_validated_yield_or_mechanism"] is True
        for species in ("N2", "N", "NO", "CH3", "H2CN", "HCN", "CH2O"):
            tracked = pathway_audit["species"][species]
            assert tracked["integrated_chemical_net_mol"] == pytest.approx(
                tracked["observed_chemical_net_mol"], abs=1e-11
            )
        amounts.append((
            result["final_compartments_mol"]["collected_condensate"]["HCN"],
            result["final_compartments_mol"]["collected_condensate"]["CH2O"],
            charged["final_compartments_mol"]["collected_condensate"]["HCN"],
            charged["final_compartments_mol"]["collected_condensate"]["CH2O"],
        ))
    assert all(value > 0 for value in amounts[0])
    assert amounts[1] == pytest.approx(amounts[0], rel=1e-3)
    profile["apparatus"]["primary_source_forcing"] = "impulse"
    with pytest.raises(ValueError, match="requires time-averaged"):
        PrebioticApparatusRequest.model_validate(profile)


def test_discharge_zone_cannot_exceed_gas_loop() -> None:
    payload = apparatus_payload()
    payload["apparatus"]["discharge_zone_volume"] = evidence(2.1, "L")
    with pytest.raises(ValueError, match="cannot exceed the gas loop"):
        PrebioticApparatusRequest.model_validate(payload)


def test_hybrid_apparatus_includes_n2_but_does_not_invent_nitrogen_fragments() -> None:
    payload = apparatus_payload()
    payload["apparatus"]["cycles"] = 1
    payload["species"].extend(
        [
            {
                "id": "H2",
                "name": "Hydrogen",
                "composition": {"H": 2},
                "charge": 0,
                "smiles": "[H][H]",
                "phase_behavior": "permanent_gas",
            },
            {
                "id": "N2",
                "name": "Nitrogen",
                "composition": {"N": 2},
                "charge": 0,
                "smiles": "N#N",
                "phase_behavior": "permanent_gas",
            },
        ]
    )
    payload["initial_inventory"].extend(
        [
            {"species_id": "H2", "compartment": "gas_loop", "amount": {"value": 0.1, "unit": "mol"}},
            {"species_id": "N2", "compartment": "gas_loop", "amount": {"value": 0.1, "unit": "mol"}},
        ]
    )
    payload["spark"] = boltzmann_spark()
    payload["spark"]["boltzmann"]["collision_set"] = HYBRID_DATASET_ID
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    snapshot = result["cycles"][0]["plasma_snapshot"]
    assert "N2" in snapshot["mole_fractions"]
    assert snapshot["collision_set"]["dataset_id"] == HYBRID_DATASET_ID
    assert any(item["target"] == "N2" for item in snapshot["source_terms"])
    accepted = result["cycles"][0]["heavy_species_snapshot"]["source_extents"]
    assert all(item["target"] != "N2" for item in accepted)
    assert result["material_ledger"]["closed"] is True


def test_tabata_total_nitrogen_rate_does_not_invent_ground_state_fragments() -> None:
    payload = apparatus_payload()
    payload["apparatus"]["cycles"] = 1
    payload["species"].extend(
        [
            {
                "id": "H2",
                "name": "Hydrogen",
                "composition": {"H": 2},
                "charge": 0,
                "smiles": "[H][H]",
                "phase_behavior": "permanent_gas",
            },
            {
                "id": "N2",
                "name": "Nitrogen",
                "composition": {"N": 2},
                "charge": 0,
                "smiles": "N#N",
                "phase_behavior": "permanent_gas",
            },
        ]
    )
    payload["initial_inventory"].extend(
        [
            {"species_id": "H2", "compartment": "gas_loop", "amount": {"value": 0.1, "unit": "mol"}},
            {"species_id": "N2", "compartment": "gas_loop", "amount": {"value": 0.1, "unit": "mol"}},
        ]
    )
    payload["spark"] = boltzmann_spark()
    payload["spark"]["boltzmann"]["collision_set"] = REACTIVE_N2_DATASET_ID
    payload["spark"]["boltzmann"]["energy_cells"] = 200
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    snapshot = result["cycles"][0]["plasma_snapshot"]
    dissociation = next(
        item
        for item in snapshot["source_terms"]
        if item["target"] == "N2" and item["product_label"] == "N + N"
    )
    assert dissociation["role"] == "dissociation"
    assert dissociation["atoms_and_product_charge_balanced"] is True
    assert dissociation["conditional_source_rate_mol_s"] > 0
    assert dissociation["outcome_only"] is True
    assert dissociation["product_state_unresolved"] is True
    assert dissociation["threshold_energy_lower_bound_J"] > 0
    assert dissociation["independent_threshold_energy_lower_bound_J"] == 0
    assert snapshot["overlapping_outcome_threshold_energy_J"] == pytest.approx(
        dissociation["threshold_energy_lower_bound_J"]
    )
    assert snapshot["threshold_energy_lower_bound_J"] == pytest.approx(
        sum(source["independent_threshold_energy_lower_bound_J"]
            for source in snapshot["source_terms"])
    )
    accepted = result["cycles"][0]["heavy_species_snapshot"]["source_extents"]
    assert all(item["target"] != "N2" for item in accepted)
    allocation = result["cycles"][0]["heavy_species_snapshot"]["collision_energy_allocation"]
    assert any(
        item["channel_id"] == dissociation["id"]
        and item["reason"] == "unresolved_electronic_product_state"
        for item in allocation["excluded_fragment_channels"]
    )
    assert allocation["overlapping_outcome_raw_threshold_J"] > 0
    assert allocation["all_collision_raw_threshold_demand_J"] > allocation["resolved_fragment_raw_threshold_demand_J"]
    assert allocation["resolved_fragment_energy_used_J"] < allocation["allocated_collision_energy_J"]
    shadow = result["cycles"][0]["heavy_species_snapshot"]["mechanism_comparison"]["shadow"]
    assert shadow["status"] == "correlated_pavan_neutral_shadow_complete"
    assert shadow["mapped_excitation_source_extents"]
    assert {item["product"] for item in shadow["mapped_excitation_source_extents"]} >= {
        "N2(A)", "N2(B)", "N2(C)", "N2(v1)"
    }
    assert shadow["mapped_excitation_energy_used_J"] > 0
    assert shadow["total_mapped_source_energy_used_J"] <= allocation["allocated_collision_energy_J"]
    assert shadow["scientific_firewall"]["primary_inventory_mutated_by_shadow"] is False
    charged = result["cycles"][0]["heavy_species_snapshot"]["mechanism_comparison"]["charged_shadow"]
    assert charged["status"] == "correlated_pavan_charge_balanced_shadow_complete"
    assert charged["retained_reactions"] == 406
    assert charged["retained_two_temperature_reactions"] == 19
    assert charged["disabled_embedded_electron_collision_reactions"] == 34
    assert any(
        item["channel_id"] == dissociation["id"]
        and item["reason"] == "electronic product-state branching unresolved"
        for item in charged["unsupported_source_channels"]
    )
    assert charged["source_extents"]
    assert charged["realized_independent_threshold_energy_J"] <= (
        charged["allocated_collision_energy_J"] * (1 + 1e-7)
    )
    assert charged["threshold_energy_is_full_thermal_balance"] is False
    assert all(
        item["accepted_extent_mol"] >= 0
        and item["nominal_inventory_capped_extent_mol"] > 0
        and item["time_averaged_source_rate_mol_s"]
        == pytest.approx(item["accepted_extent_mol"] / charged["duration_s"])
        for item in charged["source_extents"]
    )
    assert abs(charged["source_charge_residual_mol_e"]) < 1e-14
    assert abs(charged["final_charge_residual_mol_e"]) < 1e-13
    assert charged["max_atom_residual_mol"] < 1e-12
    assert charged["scientific_firewall"]["embedded_electron_collision_rates_used"] is False
    assert charged["scientific_firewall"]["external_eedf_sources_double_counted"] is False
    assert charged["scientific_firewall"]["primary_inventory_mutated_by_shadow"] is False
    assert charged["final_instantaneous_pathway_diagnostics"]["HCN"]["formation"]
    assert charged["final_instantaneous_pathway_diagnostics"]["HCN"]["consumption"]
    assert result["material_ledger"]["closed"] is True


def test_reviewed_phase_partition_conserves_each_species() -> None:
    prepared = {
        "ids": ["H2O", "HCN"],
        "index": {"H2O": 0, "HCN": 1},
    }
    gas = np.asarray([0.0, 1e-5])
    liquid = np.asarray([10.0, 0.0])
    cold_gas, cold_liquid, report = equilibrate_gas_water(
        gas, liquid, prepared, 280.0, 101_325.0, 5.0
    )
    assert cold_gas + cold_liquid == pytest.approx(gas + liquid)
    assert cold_liquid[1] > 0
    assert report["dataset_id"] == "sander_nist_henry_prebiotic_subset_v3"
    assert report["partial_pressure_basis"] == "fixed_volume_ideal_gas_nRT_over_V"


def test_reviewed_carbon_dioxide_henry_transfer_reaches_water() -> None:
    prepared = {
        "ids": ["H2O", "CO2"],
        "index": {"H2O": 0, "CO2": 1},
    }
    gas = np.asarray([0.0, 1e-3])
    liquid = np.asarray([55.5, 0.0])
    outlet_gas, outlet_liquid, report = equilibrate_gas_water(
        gas, liquid, prepared, 298.15, 101_325.0, 5.0
    )
    assert outlet_gas + outlet_liquid == pytest.approx(gas + liquid)
    assert 0 < outlet_liquid[1] < gas[1]
    assert report["dataset_id"] == "sander_nist_henry_prebiotic_subset_v3"


def test_henry_partition_uses_current_fixed_volume_not_stale_pressure() -> None:
    prepared = {"ids": ["H2O", "HCN"], "index": {"H2O": 0, "HCN": 1}}
    gas = np.asarray([0.0, 1e-5])
    liquid = np.asarray([10.0, 0.0])
    _, low_liquid, low_report = equilibrate_gas_water(
        gas, liquid, prepared, 298.15, 10_000.0, 5.0
    )
    _, high_liquid, high_report = equilibrate_gas_water(
        gas, liquid, prepared, 298.15, 1_000_000.0, 5.0
    )
    assert low_liquid == pytest.approx(high_liquid)
    assert low_report["actual_ideal_gas_pressure_Pa"] == pytest.approx(
        high_report["actual_ideal_gas_pressure_Pa"]
    )
    _, small_volume_liquid, _ = equilibrate_gas_water(
        gas, liquid, prepared, 298.15, 10_000.0, 2.5
    )
    assert small_volume_liquid[1] > low_liquid[1]


def test_flow_limited_henry_capture_preserves_material_and_elapsed_flow() -> None:
    prepared = {"ids": ["H2O", "N2", "HCN"], "index": {"H2O": 0, "N2": 1, "HCN": 2}}
    original_gas = np.asarray([0.0, 0.1, 1e-5])
    original_liquid = np.zeros(3)
    def run(steps: int) -> tuple[np.ndarray, np.ndarray]:
        gas, liquid = original_gas.copy(), original_liquid.copy()
        for _ in range(steps):
            gas, liquid, report = capture_gas_water_by_flow(
                gas, liquid, prepared, 280.0, 101_325.0, 5.0, 0.2,
                10.0 / steps, 0.1, 10.0,
            )
            assert report["status"] == "flow_limited_henry_capture_applied"
        return gas, liquid
    coarse_gas, coarse_liquid = run(1)
    fine_gas, fine_liquid = run(16)
    assert coarse_gas + coarse_liquid == pytest.approx(original_gas + original_liquid)
    assert fine_gas + fine_liquid == pytest.approx(original_gas + original_liquid)
    assert fine_liquid[2] == pytest.approx(coarse_liquid[2], rel=1e-4)
    still_gas, still_liquid, report = capture_gas_water_by_flow(
        original_gas, original_liquid, prepared, 280.0, 101_325.0, 5.0, 0.0,
        10.0, 0.1, 10.0,
    )
    assert report["status"] == "no_flow_or_contact_water"
    assert still_gas == pytest.approx(original_gas)
    assert still_liquid == pytest.approx(original_liquid)


def test_flow_capture_ignores_stale_pressure_for_fixed_volume() -> None:
    prepared = {"ids": ["H2O", "N2", "HCN"], "index": {"H2O": 0, "N2": 1, "HCN": 2}}
    gas = np.asarray([0.0, 0.1, 1e-5])
    liquid = np.zeros(3)
    results = [
        capture_gas_water_by_flow(
            gas, liquid, prepared, 280.0, pressure, 5.0, 0.2, 10.0, 0.1, 10.0
        )
        for pressure in (10_000.0, 1_000_000.0)
    ]
    assert results[0][0] == pytest.approx(results[1][0])
    assert results[0][1] == pytest.approx(results[1][1])
    assert results[0][2]["partial_pressure_basis"] == "fixed_volume_ideal_gas_nRT_over_V"


def test_water_vapor_capacity_returns_supersaturated_water_to_liquid() -> None:
    prepared = {"ids": ["H2O", "CH4"], "index": {"H2O": 0, "CH4": 1}}
    gas = np.asarray([10.0, 0.1])
    liquid = np.asarray([20.0, 0.0])
    capped_gas, capped_liquid, report = enforce_water_vapor_capacity(
        gas, liquid, prepared, 373.15, 5.0
    )
    assert capped_gas + capped_liquid == pytest.approx(gas + liquid)
    assert capped_gas[0] == pytest.approx(report["maximum_water_vapor_mol"])
    assert report["water_partial_pressure_Pa"] == pytest.approx(
        report["saturation_pressure_Pa"]
    )
    assert report["recondensed_water_mol"] > 9.0


def test_blind_aqueous_enumerator_generates_before_motif_classification() -> None:
    candidates = enumerate_neutral_tree_isomers({"C": 2, "H": 5, "N": 1, "O": 2})
    assert candidates
    assert all(item["formula"] == "C2H5NO2" for item in candidates)
    assert any("alpha_amino_carboxylic_acid" in item["motifs"] for item in candidates)

    prepared = {
        "ids": ["H2O", "HCN", "CH2O"],
        "index": {"H2O": 0, "HCN": 1, "CH2O": 2},
        "identities": [
            {"id": "H2O", "composition": {"H": 2, "O": 1}},
            {"id": "HCN", "composition": {"H": 1, "C": 1, "N": 1}},
            {"id": "CH2O", "composition": {"H": 2, "C": 1, "O": 1}},
        ],
    }
    report = discover_aqueous_candidates(np.asarray([1.0, 1e-6, 1e-6]), prepared)
    assert report["status"] == "product_agnostic_candidates_enumerated"
    assert report["scientific_firewall"]["known_products_accessible_to_generator"] is False
    assert report["scientific_firewall"]["reaction_templates_accessible_to_generator"] is False
    assert report["post_discovery_motif_counts"]["alpha_amino_carboxylic_acid"] > 0


def test_blind_enumerator_preserves_both_alanine_enantiomers() -> None:
    candidates = enumerate_neutral_tree_isomers({"C": 3, "H": 7, "N": 1, "O": 2})
    alpha_amino = {
        item["smiles"]: {center["descriptor"] for center in item["stereochemistry"]}
        for item in candidates
        if "alpha_amino_carboxylic_acid" in item["motifs"]
    }
    assert "C[C@H](N)C(=O)O" in alpha_amino
    assert "C[C@@H](N)C(=O)O" in alpha_amino
    assert {next(iter(value)) for value in alpha_amino.values() if value} >= {"R", "S"}


def test_post_discovery_energy_filter_is_generic_and_keeps_beta_alanine() -> None:
    records = []
    for composition in (
        {"C": 2, "H": 5, "N": 1, "O": 2},
        {"C": 3, "H": 7, "N": 1, "O": 2},
    ):
        candidates = enumerate_neutral_tree_isomers(composition)
        records.append({"formula": candidates[0]["formula"], "candidates": candidates})
    selected = select_amino_carboxylic_candidates({"formula_records": records})
    smiles = {item["smiles"] for item in selected}
    assert len(selected) == 7
    assert "NCC(=O)O" in smiles
    assert "NCCC(=O)O" in smiles
    assert "C[C@H](N)C(=O)O" in smiles
    assert "C[C@@H](N)C(=O)O" in smiles


def test_enantiomer_starting_geometries_are_exact_reflections() -> None:
    candidates = [
        {
            "formula": "C3H7NO2",
            "smiles": "C[C@H](N)C(=O)O",
            "stereochemistry": [{"atom_index": 1, "descriptor": "S"}],
        },
        {
            "formula": "C3H7NO2",
            "smiles": "C[C@@H](N)C(=O)O",
            "stereochemistry": [{"atom_index": 1, "descriptor": "R"}],
        },
    ]
    geometries, pairs = build_paired_starting_geometries(candidates, 20260923)
    left = geometries[candidates[0]["smiles"]].get_positions()
    right = geometries[candidates[1]["smiles"]].get_positions()
    assert len(pairs) == 1
    assert np.abs(left[:, 0]) == pytest.approx(np.abs(right[:, 0]))
    assert left[:, 1:] == pytest.approx(right[:, 1:])


def test_candidate_ranking_job_survives_manager_restart(tmp_path: Path) -> None:
    manager = CandidateRankingJobManager(tmp_path)
    created = manager.create()
    restored = CandidateRankingJobManager(tmp_path).get(created.id)
    assert restored is not None
    assert restored.id == created.id
    assert restored.status == "queued"


def test_post_discovery_evidence_cannot_change_frozen_discovery() -> None:
    prepared = {
        "ids": ["H2O", "NH3", "HCN", "CH2O"],
        "index": {"H2O": 0, "NH3": 1, "HCN": 2, "CH2O": 3},
        "identities": [
            {"id": "H2O", "composition": {"H": 2, "O": 1}},
            {"id": "NH3", "composition": {"H": 3, "N": 1}},
            {"id": "HCN", "composition": {"H": 1, "C": 1, "N": 1}},
            {"id": "CH2O", "composition": {"H": 2, "C": 1, "O": 1}},
        ],
    }
    inventory = np.asarray([10.0, 0.1, 1e-4, 1e-4])
    discovery = discover_aqueous_candidates(inventory, prepared)
    digest_before, bytes_before = freeze_discovery_artifact(discovery)
    report = evaluate_frozen_aqueous_discovery(
        discovery,
        inventory,
        prepared,
        10.0 * 0.01801528,
        7 * 24 * 3600,
        (280.0, 373.15),
    )
    digest_after, bytes_after = freeze_discovery_artifact(discovery)
    assert digest_before == digest_after == report["discovery_artifact_sha256"]
    assert bytes_before == bytes_after
    assert report["candidate_product"]["present_in_frozen_discovery"] is True
    assert report["scientific_firewall"]["evaluation_feedback_to_discovery"] is False
    assert report["scientific_firewall"]["historical_benchmark_yields_used"] is False
    assert report["scientific_firewall"]["quantitative_yield_claim_allowed"] is False


def test_post_discovery_cycle_replay_accumulates_without_feedback() -> None:
    prepared = {
        "ids": ["H2O", "NH3", "HCN", "CH2O"],
        "index": {"H2O": 0, "NH3": 1, "HCN": 2, "CH2O": 3},
        "identities": [
            {"id": "H2O", "composition": {"H": 2, "O": 1}},
            {"id": "NH3", "composition": {"H": 3, "N": 1}},
            {"id": "HCN", "composition": {"H": 1, "C": 1, "N": 1}},
            {"id": "CH2O", "composition": {"H": 2, "C": 1, "O": 1}},
        ],
    }
    inventory = np.asarray([10.0, 0.1, 1e-4, 1e-4])
    discovery = discover_aqueous_candidates(inventory, prepared)
    snapshots = []
    for cycle in range(1, 4):
        fraction = cycle / 3
        snapshots.append(
            {
                "cycle": cycle,
                "start_time_s": (cycle - 1) * 100_000,
                "end_time_s": cycle * 100_000,
                "compartments_mol": {
                    "boiler_liquid": {
                        "H2O": 10.0,
                        "NH3": 0.1,
                        "HCN": 1e-4 * fraction,
                        "CH2O": 1e-4 * fraction,
                    },
                    "condensate_trap": {},
                },
            }
        )
    report = evaluate_frozen_aqueous_discovery(
        discovery,
        inventory,
        prepared,
        10.0 * 0.01801528,
        300_000,
        (280.0, 373.15),
        snapshots,
    )
    nominal = report["hydrolysis_sensitivity"]["scenarios"]["nominal"]
    assert report["hydrolysis_sensitivity"]["cycles_replayed"] == 3
    assert len(nominal["cycle_trace"]) == 3
    assert nominal["candidate_product_mol"] > 0
    assert nominal["candidate_product_mol"] < nominal["cumulative_aminonitrile_formed_mol"]
    assert nominal["cycle_trace"][-1]["candidate_product_mol"] >= nominal["cycle_trace"][0]["candidate_product_mol"]


def test_historical_benchmark_is_held_out_and_does_not_claim_reproduction() -> None:
    candidates = {
        "formula_records": [
            {
                "candidates": [
                    {"smiles": "NCC(=O)O"},
                    {"smiles": "C[C@H](N)C(=O)O"},
                    {"smiles": "NCCC(=O)O"},
                ]
            }
        ]
    }
    report = evaluate_historical_benchmark(
        candidates,
        "a" * 64,
        {"H2O": 200 / 18.01528, "H2": 0.01, "CH4": 0.02, "NH3": 0.02},
        {"HCN": 1e-5, "CH2O": 1e-5, "CH3CHO": 1e-6},
        5.0,
        1600.0,
    )
    assert report["enumeration_coverage"] == 1.0
    assert report["enumeration_coverage_is_prediction_score"] is False
    assert report["historical_conditions_matched"] is False
    assert report["validated_reproduction"] is False
    assert report["scientific_firewall"]["benchmark_observations_available_to_discovery"] is False
    assert report["scientific_firewall"]["benchmark_yields_used_for_fitting"] is False


def test_frozen_class_ranking_is_scored_only_after_prediction_hash() -> None:
    ranking = {
        "status": "frozen_candidate_energy_ranking_complete",
        "discovery_artifact_sha256": "a" * 64,
        "candidate_filter": "neutral primary/secondary amine plus carboxylic acid",
        "engines_succeeded": ["fixture_a", "fixture_b"],
        "formula_local_consensus": [
            {"formula": "C2H5NO2", "smiles": "NCC(=O)O", "mean_formula_local_rank": 2.0, "rank_range": [2, 2]},
            {"formula": "C3H7NO2", "smiles": "C[C@H](N)C(=O)O", "mean_formula_local_rank": 2.5, "rank_range": [2, 3]},
            {"formula": "C3H7NO2", "smiles": "C[C@@H](N)C(=O)O", "mean_formula_local_rank": 2.5, "rank_range": [2, 3]},
            {"formula": "C3H7NO2", "smiles": "NCCC(=O)O", "mean_formula_local_rank": 4.0, "rank_range": [4, 4]},
        ],
    }
    report = evaluate_frozen_candidate_ranking(ranking, top_k=3)
    assert len(report["prediction_artifact_sha256"]) == 64
    assert report["class_filter_coverage"] == 1.0
    assert report["top_k_formula_local_recall"] == pytest.approx(2 / 3)
    assert report["validated_reproduction"] is False
    assert report["scientific_firewall"]["prediction_frozen_before_benchmark_read"] is True
    assert report["scientific_firewall"]["global_blind_product_prediction_claimed"] is False


def test_precursor_energy_benchmark_refuses_mismatched_feed_score() -> None:
    report = evaluate_precursor_energy_benchmark(
        {"CH4": 1.0, "H2": 0.5, "H2O": 10.0, "NH3": 1.0},
        {"HCN": 1e-8, "CH2O": 2e-8},
        10.0,
    )
    assert report["comparable"] is False
    assert report["species_comparison"]["HCN"]["simulated_nmol_per_J"] == pytest.approx(1.0)
    assert report["species_comparison"]["HCN"]["ratio_to_reference"] is None
    assert report["scientific_firewall"]["mismatched_feed_scored"] is False


def test_precursor_energy_benchmark_still_refuses_mismatched_measurement_basis() -> None:
    report = evaluate_precursor_energy_benchmark(
        {"CH4": 1.0, "H2": 1.0, "H2O": 10.0, "N2": 1.0},
        {"HCN": 1e-7, "CH2O": 2e-8},
        10.0,
    )
    assert report["composition_matched"] is True
    assert report["comparable"] is False
    assert report["feed_checks"]["zero_time_yield_basis_matched"] is False
    assert report["feed_checks"]["calorimetric_energy_basis_matched"] is False
    assert report["species_comparison"]["HCN"]["ratio_to_reference"] is None


def test_precursor_energy_benchmark_keeps_exact_fit_diagnostic_until_protocol_known() -> None:
    cycles = []
    energy_per_cycle_J = 10.0
    for cycle in range(1, 169):
        time_s = cycle * 3600.0
        cumulative_energy_J = cycle * energy_per_cycle_J
        hcn_yield = 10.0 - 1e-6 * time_s
        ch2o_yield = 2.0 - 1e-6 * time_s
        cycles.append(
            {
                "cycle": cycle,
                "end_time_s": time_s,
                "deposited_energy_used_J": energy_per_cycle_J,
                "compartments_mol": {
                    "boiler_liquid": {
                        "HCN": hcn_yield * cumulative_energy_J / 1e9,
                        "CH2O": ch2o_yield * cumulative_energy_J / 1e9,
                    },
                    "condensate_trap": {},
                },
            }
        )
    report = evaluate_precursor_energy_benchmark(
        {"CH4": 1.0, "H2": 1.0, "H2O": 10.0, "N2": 1.0},
        cycles[-1]["compartments_mol"]["boiler_liquid"],
        168 * energy_per_cycle_J,
        cycles=cycles,
        energy_basis="calorimetric_measurement",
        calorimetric_evidence={
            "value": energy_per_cycle_J,
            "unit": "J",
            "source": "synthetic calorimetric validation fixture",
            "method": "resistance-heater calorimetry",
            "basis": "experimental",
        },
    )
    assert report["comparable"] is False
    assert report["reference_protocol_complete"] is False
    assert report["reference_protocol_match_assessed"] is False
    assert report["status"] == "diagnostic_only_reference_protocol_incomplete"
    assert "exact_sample_times_and_regression_weighting" in report["missing_reference_protocol_fields"]
    assert report["feed_checks"]["zero_time_yield_basis_matched"] is True
    assert report["feed_checks"]["calorimetric_energy_basis_matched"] is True
    assert report["zero_time_extrapolation"]["HCN"]["sample_count"] == 168
    assert report["illustrative_sparse_sampling_sensitivity"]["zero_time_extrapolation"]["HCN"]["sample_count"] == 5
    assert report["illustrative_sparse_sampling_sensitivity"]["zero_time_extrapolation"]["HCN"]["intercept_nmol_per_J"] == pytest.approx(10.0)
    assert report["scientific_firewall"]["sparse_sampling_sensitivity_used_for_simulation_or_scoring"] is False
    assert report["species_comparison"]["HCN"]["simulated_nmol_per_J"] == pytest.approx(10.0)
    assert report["species_comparison"]["CH2O"]["simulated_nmol_per_J"] == pytest.approx(2.0)
    assert report["species_comparison"]["HCN"]["ratio_to_reference"] is None
    assert report["scientific_firewall"]["incomplete_reference_protocol_scored"] is False

    pressure_mismatched = evaluate_precursor_energy_benchmark(
        {"CH4": 1.0, "H2": 1.0, "H2O": 10.0, "N2": 1.0},
        cycles[-1]["compartments_mol"]["boiler_liquid"],
        168 * energy_per_cycle_J,
        cycles=cycles,
        energy_basis="calorimetric_measurement",
        calorimetric_evidence={
            "value": energy_per_cycle_J,
            "unit": "J",
            "source": "synthetic calorimetric validation fixture",
            "method": "resistance-heater calorimetry",
            "basis": "experimental",
        },
        declared_reference_pressure_Pa=101325.0,
        simulated_gas_pressure_range_Pa=[202000.0, 204000.0],
    )
    assert pressure_mismatched["pressure_self_consistent"] is False
    assert pressure_mismatched["comparable"] is False
    assert pressure_mismatched["species_comparison"]["HCN"]["ratio_to_reference"] is None
    assert any("inventory-derived gas pressure" in item for item in pressure_mismatched["remaining_blockers"])

    pressure_unassessed = evaluate_precursor_energy_benchmark(
        {"CH4": 1.0, "H2": 1.0, "H2O": 10.0, "N2": 1.0},
        cycles[-1]["compartments_mol"]["boiler_liquid"],
        168 * energy_per_cycle_J,
        cycles=cycles,
        energy_basis="calorimetric_measurement",
        calorimetric_evidence={"basis": "experimental"},
        declared_reference_pressure_Pa=101325.0,
    )
    assert pressure_unassessed["pressure_self_consistent"] is None
    assert pressure_unassessed["comparable"] is False


def test_calorimetric_energy_overrides_electrical_estimate_without_hiding_it() -> None:
    payload = apparatus_payload()
    payload["spark"] = spark()
    payload["spark"]["calorimetric_deposited_energy_per_cycle"] = evidence(
        5.0, "J", "resistance-heater calorimetry fixture"
    ) | {"basis": "experimental", "method": "resistance-heater temperature-rise comparison"}
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert result["energy_ledger"]["deposited_energy_basis"] == "calorimetric_measurement"
    assert result["energy_ledger"]["calorimetric_deposited_energy_J"] == pytest.approx(15.0)
    assert result["energy_ledger"]["deposited_energy_used_J"] == pytest.approx(15.0)
    assert result["energy_ledger"]["deposited_electrical_energy_J"] != pytest.approx(15.0)
    assert result["cycles"][0]["deposited_energy_used_J"] == pytest.approx(5.0)


def test_historical_profile_is_versioned_valid_and_matches_material_time_targets() -> None:
    profile = historical_material_time_profile()
    request = PrebioticApparatusRequest.model_validate(profile)
    assert request.apparatus.cycles == 168
    assert request.apparatus.operation_mode == "continuous_overlap"
    assert request.apparatus.integration_substeps_per_cycle == 2
    cycle_s = sum(
        quantity.as_unit("s")
        for quantity in (
            request.apparatus.boil_duration,
            request.apparatus.discharge_duration,
            request.apparatus.condense_duration,
        )
    )
    assert request.apparatus.cycles * cycle_s == pytest.approx(7 * 24 * 3600)
    response = TestClient(app).get("/api/prebiotic/benchmarks/miller-urey-1953/profile")
    assert response.status_code == 200
    assert response.json() == profile


def test_stribling_profile_spans_measurement_window_without_invented_calorimetry() -> None:
    profile = stribling_miller_measurement_profile()
    request = PrebioticApparatusRequest.model_validate(profile)
    assert request.apparatus.cycles == 168
    assert request.apparatus.operation_mode == "continuous_overlap"
    assert request.apparatus.integration_substeps_per_cycle == 2
    cycle_s = sum(
        quantity.as_unit("s")
        for quantity in (
            request.apparatus.boil_duration,
            request.apparatus.discharge_duration,
            request.apparatus.condense_duration,
        )
    )
    assert cycle_s == pytest.approx(3600)
    assert request.apparatus.cycles * cycle_s == pytest.approx(7 * 24 * 3600)
    inventory = {
        item.species_id: item.amount.as_unit("mol") for item in request.initial_inventory
    }
    assert inventory["H2"] / inventory["CH4"] == pytest.approx(1.0)
    assert {"CH4", "N2", "H2O", "H2"} <= set(inventory)
    assert request.spark.calorimetric_deposited_energy_per_cycle is None
    assert request.spark.boltzmann is not None
    assert request.spark.boltzmann.collision_set == REACTIVE_N2_DATASET_ID
    response = TestClient(app).get(
        "/api/prebiotic/benchmarks/stribling-miller-1986/profile"
    )
    assert response.status_code == 200
    assert response.json() == profile


def test_zero_flow_control_prevents_apparatus_transfer() -> None:
    payload = apparatus_payload()
    payload["apparatus"]["gas_circulation_rate"] = {"value": 0, "unit": "L/s"}
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert result["numerical_controls"]["boil_pass_fraction"] == 0
    assert result["cycles"][0]["evaporated_mol"]["H2O"] == 0
    assert result["final_compartments_mol"]["boiler_liquid"]["H2O"] == pytest.approx(1)


def test_condenser_off_control_does_not_accumulate_condensate() -> None:
    payload = apparatus_payload()
    payload["species"][0].pop("condensation_rate")
    payload["apparatus"]["condensate_return_fraction"] = 0
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert result["cycles"][0]["evaporated_mol"]["H2O"] > 0
    assert result["cycles"][0]["condensed_mol"]["H2O"] == 0
    assert result["final_compartments_mol"]["condensate_trap"]["H2O"] == 0
    assert result["material_ledger"]["closed"] is True


def test_smiles_must_match_declared_formula_and_charge() -> None:
    payload = apparatus_payload()
    payload["species"][0]["composition"] = {"H": 2}
    with pytest.raises(ValueError, match="SMILES disagrees"):
        simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))


def test_reduced_discharge_requires_spark_to_activate_channel() -> None:
    payload = reaction_payload()
    off = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert off["final_total_mol"]["B"] == 0
    payload["spark"] = spark()
    on = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert on["final_total_mol"]["B"] > 0
    assert on["energy_ledger"]["deposited_electrical_energy_J"] > 0
    assert on["energy_ledger"]["status"] == "within_declared_deposited_energy"
    assert on["scientific_firewall"]["quantitative_yield_claim_allowed"] is False


def test_discharge_reaction_is_stopped_by_deposited_energy_budget() -> None:
    payload = reaction_payload()
    payload["spark"] = spark()
    payload["reactions"][0]["rate_law"]["energy_cost"] = evidence(1e9, "J/mol")
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    maximum_extent = (
        result["energy_ledger"]["deposited_electrical_energy_J"] / 1e9
    )
    assert result["reaction_extents_mol"]["electron_channel"] <= maximum_extent * (1 + 1e-7)
    assert result["energy_ledger"]["reaction_was_energy_limited"] is True
    assert result["energy_ledger"]["status"] == "energy_limited_to_declared_deposit"


def test_unbalanced_discharge_channel_is_rejected() -> None:
    payload = reaction_payload()
    payload["species"][1]["composition"] = {"H": 1}
    with pytest.raises(ValueError, match="does not conserve"):
        simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))


def test_mid_run_addition_and_sampling_are_in_material_ledger() -> None:
    payload = apparatus_payload()
    payload["interventions"] = [
        {
            "cycle": 2,
            "additions": [
                {
                    "species_id": "CH4",
                    "compartment": "gas_loop",
                    "amount": {"value": 0.02, "unit": "mol"},
                }
            ],
            "sample_compartment": "gas_loop",
            "sample_fraction": 0.1,
            "note": "Add feed after the first complete circulation cycle",
        }
    ]
    result = simulate_prebiotic_apparatus(
        PrebioticApparatusRequest.model_validate(payload)
    )
    assert result["material_ledger"]["added_mol"]["CH4"] == pytest.approx(0.02)
    assert result["material_ledger"]["sampled_mol"]["CH4"] > 0
    assert result["material_ledger"]["closed"] is True
    assert result["cycles"][1]["intervention_notes"]


def test_achiral_parallel_network_does_not_create_numerical_enantiomer_bias() -> None:
    payload = reaction_payload()
    payload["species"] = [
        {
            "id": "P",
            "name": "Achiral precursor",
            "composition": {"C": 3, "H": 7, "N": 1, "O": 2},
            "charge": 0,
            "phase_behavior": "permanent_gas",
        },
        {
            "id": "L",
            "name": "Left stereoisomer",
            "composition": {"C": 3, "H": 7, "N": 1, "O": 2},
            "charge": 0,
            "smiles": "N[C@@H](C)C(=O)O",
            "phase_behavior": "nonvolatile",
            "enantiomer_of": "D",
        },
        {
            "id": "D",
            "name": "Right stereoisomer",
            "composition": {"C": 3, "H": 7, "N": 1, "O": 2},
            "charge": 0,
            "smiles": "N[C@H](C)C(=O)O",
            "phase_behavior": "nonvolatile",
            "enantiomer_of": "L",
        },
    ]
    payload["reactions"] = [
        {
            "id": side,
            "name": f"Symmetric {side} channel",
            "phase": "gas",
            "stoichiometry": {"P": -1, side: 1},
            "rate_law": {
                "kind": "electron_impact_reduced",
                "rate_constant": evidence(0.5, "1/s"),
                "energy_cost": evidence(500, "J/mol"),
            },
        }
        for side in ("L", "D")
    ]
    payload["initial_inventory"] = [
        {
            "species_id": "P",
            "compartment": "gas_loop",
            "amount": {"value": 0.1, "unit": "mol"},
        }
    ]
    payload["spark"] = spark()
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert result["final_total_mol"]["L"] == pytest.approx(result["final_total_mol"]["D"])
    assert result["stereochemistry_diagnostics"][0]["racemic_within_numerical_tolerance"] is True


def test_blind_discovery_origin_requires_frozen_manifest() -> None:
    payload = apparatus_payload()
    payload["mechanism_origin"] = "blind_discovery_artifact"
    with pytest.raises(ValueError, match="manifest"):
        PrebioticApparatusRequest.model_validate(payload)
    payload["discovery_manifest_sha256"] = "a" * 64
    request = PrebioticApparatusRequest.model_validate(payload)
    result = simulate_prebiotic_apparatus(request)
    assert result["scientific_firewall"]["blind_discovery_eligible"] is True


def test_frozen_manifest_does_not_make_gri_plasma_product_blind() -> None:
    payload = apparatus_payload()
    payload["apparatus"]["cycles"] = 1
    payload["spark"] = boltzmann_spark()
    payload["mechanism_origin"] = "blind_discovery_artifact"
    payload["discovery_manifest_sha256"] = "a" * 64
    result = simulate_prebiotic_apparatus(
        PrebioticApparatusRequest.model_validate(payload)
    )
    assert result["scientific_firewall"]["blind_discovery_eligible"] is False
    assert result["scientific_firewall"]["preexisting_gri_hcn_reactions_accessible_to_solver"] is True
    assert any("not de novo product-identity discovery" in warning for warning in result["warnings"])


def test_prebiotic_api_and_cantera_plasma_capability() -> None:
    client = TestClient(app)
    capability = client.get("/api/prebiotic/capabilities")
    assert capability.status_code == 200
    data = capability.json()
    assert data["apparatus"]["available"] is True
    assert data["cantera_cross_section_plasma"]["electron_collision_rate"] is True
    assert data["cantera_cross_section_plasma"]["reference_example_loaded"] is True
    assert data["cantera_cross_section_plasma"]["pavan_methane_plasma_loaded"] is True
    assert data["cantera_cross_section_plasma"]["pavan_methane_plasma_species"] == 71
    assert data["cantera_cross_section_plasma"]["pavan_methane_plasma_reactions"] == 440
    assert data["cantera_cross_section_plasma"]["pavan_charge_balanced_shadow_connected"] is True
    assert data["cantera_cross_section_plasma"]["pavan_retained_two_temperature_reactions"] == 19
    assert data["cantera_cross_section_plasma"]["pavan_disabled_embedded_electron_collision_reactions"] == 34
    assert data["cantera_cross_section_plasma"]["kamarinopoulou_hcn_shadow_connected"] is True
    assert data["cantera_cross_section_plasma"]["kamarinopoulou_hcn_species"] == 83
    assert data["cantera_cross_section_plasma"]["kamarinopoulou_hcn_neutral_reactions"] == 189
    assert data["cantera_cross_section_plasma"]["kamarinopoulou_hcn_fixed_electron_reactions_excluded"] == 62
    assert data["boltzmann_eedf"]["heavy_species_shadow_is_independent_confidence"] is False
    assert data["boltzmann_eedf"]["propagated_charged_branch"] is True
    assert data["boltzmann_eedf"]["propagated_charged_branch_self_consistent_eedf"] is True
    assert data["boltzmann_eedf"]["propagated_charged_branch_gas_water_partition"] is True
    assert data["boltzmann_eedf"]["continuous_overlap_operator_splitting"] is True
    assert data["boltzmann_eedf"]["automatic_timestep_convergence_claim"] is False
    response = client.post("/api/prebiotic/apparatus/simulate", json=apparatus_payload())
    assert response.status_code == 200, response.text
    assert response.json()["material_ledger"]["closed"] is True


def test_recorded_prebiotic_results_are_read_only_and_labelled_provisional() -> None:
    response = TestClient(app).get("/api/prebiotic/validation/recorded")
    assert response.status_code == 200
    record = response.json()
    assert record["status"] == "recorded_scenarios_not_validated_reproduction"
    assert any(run["cycles_completed"] == 168 for run in record["runs"])
    assert len(record["runs"]) == 6
    assert {run["condensate_nonwater_return_fraction"] for run in record["runs"]} == {0.0, 1.0}
    assert len(record["convergence_studies"]) == 6
    power_run = next(
        item for item in record["runs"]
        if item.get("model_revision") == "electron_drift_power_density_upper_bound"
    )
    assert power_run["primary"]["hcn_zero_time_intercept_nmol_per_J"] < 0.1
    assert power_run["first_step_uniform_drift_power_audit"]["drift_to_deposited_energy_ratio"] > 600
    assert power_run["full_week_resolution_check"][
        "maximum_absolute_relative_change_in_four_zero_time_intercepts"
    ] < 1e-4
    power_convergence = next(
        item for item in record["convergence_studies"]
        if item.get("model_revision") == "electron_drift_power_density_upper_bound"
    )
    assert max(abs(value) for value in power_convergence["relative_change_2_to_32"].values()) < 2e-4
    flow_study = next(
        item for item in record["convergence_studies"]
        if item.get("condenser_phase_transfer_model") == "flow_limited_capture"
    )
    assert flow_study["convergence_status"] == "not_converged"
    assert flow_study["relative_change_16_to_32"]["primary_water_hcn"] > 0.5
    coupled_study = next(
        item for item in record["convergence_studies"]
        if item.get("condenser_phase_transfer_model") == "coupled_flow_transport"
    )
    assert coupled_study["convergence_status"] == "24h_amounts_converged_at_tested_resolutions_only"
    assert max(abs(value) for value in coupled_study["relative_change_2_to_32"].values()) < 2e-4
    coupled_week = next(
        item for item in record["runs"]
        if item.get("condenser_phase_transfer_model") == "coupled_flow_transport"
    )
    assert coupled_week["full_week_resolution_check"]["maximum_absolute_relative_change_in_four_zero_time_intercepts"] < 1e-4
    assert coupled_week["quantitative_yield_claim_allowed"] is False
    assert len(record["sensitivity_studies"]) == 1
    contact = record["sensitivity_studies"][0]
    assert contact["reference_yields_used_for_solver_fitting"] is False
    assert contact["primary_hcn_high_to_low_relative_span"] > 0.5
    assert all(run["quantitative_yield_claim_allowed"] is False for run in record["runs"])
    assert record["scientific_firewall"]["reference_yields_used_for_discovery"] is False
    assert record["scientific_firewall"]["reference_yields_used_for_solver_fitting"] is False


def test_kamarinopoulou_reference_import_is_hash_verified_and_balanced() -> None:
    root = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "reference"
        / "kamarinopoulou_2024"
    )
    manifest = json.loads((root / "import_manifest.json").read_text(encoding="utf-8"))
    pdf_hash = hashlib.sha256((root / "sciadv.adl4246_sm.pdf").read_bytes()).hexdigest()
    yaml_path = root / "kamarinopoulou_2024_neutral.yaml"
    yaml_hash = hashlib.sha256(yaml_path.read_bytes()).hexdigest()
    assert pdf_hash == manifest["source_pdf_sha256"]
    assert yaml_hash == manifest["generated_files"][yaml_path.name]
    assert manifest["published_reaction_records"] == 251
    assert manifest["neutral_heavy_reaction_records"] == 189
    assert manifest["fixed_condition_electron_reaction_records"] == 62
    assert manifest["scientific_firewall"]["paper_experimental_yields_imported_into_mechanism"] is False

    phase = ct.Solution(str(yaml_path))
    assert phase.n_species == 83
    assert phase.n_reactions == 189
    for reaction in phase.reactions():
        for element in phase.element_names:
            reactant_atoms = sum(
                phase.n_atoms(species, element) * coefficient
                for species, coefficient in reaction.reactants.items()
            )
            product_atoms = sum(
                phase.n_atoms(species, element) * coefficient
                for species, coefficient in reaction.products.items()
            )
            assert reactant_atoms == pytest.approx(product_atoms), reaction.equation


def collision_payload() -> dict:
    return {
        "id": "abstract_ionization_channel",
        "target_species": "X",
        "product_species": "X+",
        "kind": "ionization",
        "energy_levels_eV": [0, 5, 10, 15, 20],
        "cross_sections_m2": [0, 0, 1e-22, 4e-21, 2e-21],
        "source": "independent cross-section validation fixture",
        "method": "tabulated electron collision measurements",
        "basis": "experimental",
        "mechanism_origin": "user_declared",
    }


def test_cross_section_dataset_is_validated_by_cantera_and_hashed() -> None:
    request = ElectronCollisionDatasetRequest.model_validate(collision_payload())
    result = validate_electron_collision_dataset(request)
    assert result["cantera_rate_type"] == "ElectronCollisionPlasmaRate"
    assert result["points"] == 5
    assert result["first_nonzero_energy_eV"] == 10
    assert len(result["sha256"]) == 64
    assert result["scientific_firewall"]["usable_for_rate_without_eedf"] is False


def test_cross_section_dataset_rejects_bad_grid_and_unknown_answer_fields() -> None:
    payload = collision_payload()
    payload["energy_levels_eV"] = [0, 5, 5, 20, 30]
    with pytest.raises(ValueError, match="strictly increasing"):
        ElectronCollisionDatasetRequest.model_validate(payload)
    payload = collision_payload()
    payload["expected_product"] = "glycine"
    response = TestClient(app).post(
        "/api/prebiotic/plasma/cross-sections/validate", json=payload
    )
    assert response.status_code == 422
    assert "extra_forbidden" in response.text


def test_cross_section_api() -> None:
    response = TestClient(app).post(
        "/api/prebiotic/plasma/cross-sections/validate", json=collision_payload()
    )
    assert response.status_code == 200, response.text
    assert response.json()["cantera_version"] == "3.2.0"


def collision_rate_payload() -> dict:
    return {
        "dataset": {
            "id": "constant_cross_section_fixture",
            "target_species": "A",
            "product_species": "A+",
            "kind": "ionization",
            "energy_levels_eV": [0, 100],
            "cross_sections_m2": [1e-20, 1e-20],
            "source": "analytic constant-cross-section validation fixture",
            "method": "constant curve for an independent closed-form test",
            "basis": "computed",
            "mechanism_origin": "user_declared",
        },
        "eedf": {
            "mode": "maxwellian",
            "mean_energy": {"value": 5, "unit": "eV"},
            "source": "analytic Maxwell-Boltzmann validation fixture",
            "method": "normalized three-dimensional energy distribution",
            "basis": "computed",
        },
        "electron_number_density": evidence(1e15, "1/m^3"),
    }


def test_maxwellian_constant_cross_section_matches_closed_form_mean_speed() -> None:
    request = ElectronCollisionRateRequest.model_validate(collision_rate_payload())
    result = calculate_electron_collision_rate(request)
    temperature_eV = 2 * 5 / 3
    mean_speed = math.sqrt(8 * elementary_charge * temperature_eV / (math.pi * electron_mass))
    expected = 1e-20 * mean_speed
    assert result["rate_coefficient_m3_s"] == pytest.approx(expected, rel=2e-3)
    assert result["pseudo_first_order_rate_s_inv"] == pytest.approx(
        result["rate_coefficient_m3_s"] * 1e15
    )
    assert result["scientific_firewall"]["quantitative_yield_claim_allowed"] is False


def test_maxwellian_rate_agrees_with_cantera_packaged_collision_mechanism() -> None:
    mechanism = next(
        Path(directory) / "example_data" / "oxygen-plasma-itikawa.yaml"
        for directory in ct.get_data_directories()
        if (Path(directory) / "example_data" / "oxygen-plasma-itikawa.yaml").is_file()
    )
    phase = ct.Solution(
        str(mechanism), "isotropic-electron-energy-plasma", transport_model=None
    )
    phase.TPX = 300, ct.one_atm, "O2:1.0, e:5e-3, O2+:5e-3"
    phase.isotropic_shape_factor = 1.0
    phase.mean_electron_energy = 5.0
    reaction_index, reaction = next(
        (index, candidate)
        for index, candidate in enumerate(phase.reactions())
        if isinstance(candidate.rate, ct.ElectronCollisionPlasmaRate)
    )
    dataset = ElectronCollisionDatasetRequest.model_validate(
        {
            "id": "cantera_oxygen_reference",
            "target_species": "O2",
            "product_species": "O+O-",
            "kind": "attachment",
            "energy_levels_eV": reaction.rate.energy_levels.tolist(),
            "cross_sections_m2": reaction.rate.cross_sections.tolist(),
            "source": "Cantera packaged oxygen-plasma-itikawa.yaml",
            "method": "packaged electron-collision reference mechanism",
            "basis": "literature",
        }
    )
    result = calculate_rate_coefficient(
        dataset,
        ElectronCollisionRateRequest.model_validate(
            {**collision_rate_payload(), "dataset": dataset.model_dump(mode="json")}
        ).eedf,
    )
    # L/mol/s and m3/kmol/s are numerically identical.
    assert result["rate_coefficient_L_mol_s"] == pytest.approx(
        phase.forward_rate_constants[reaction_index], rel=0.02
    )


def test_tabulated_eedf_is_normalized_and_reports_derived_mean() -> None:
    payload = collision_rate_payload()
    payload["eedf"] = {
        "mode": "tabulated",
        "energy_levels_eV": [0, 5, 10],
        "probability_density_eV_inv": [0, 4, 0],
        "source": "triangular EEDF validation fixture",
        "method": "tabulated probability density",
        "basis": "computed",
    }
    result = calculate_electron_collision_rate(
        ElectronCollisionRateRequest.model_validate(payload)
    )
    assert result["eedf"]["mean_energy_eV"] == pytest.approx(5)
    assert result["eedf"]["cross_section_support_probability"] == pytest.approx(1)


def electron_collision_apparatus_payload() -> dict:
    payload = apparatus_payload()
    payload["species"] = [
        {
            "id": "A",
            "name": "Neutral carbon test species",
            "composition": {"C": 1},
            "charge": 0,
            "phase_behavior": "permanent_gas",
        },
        {
            "id": "A+",
            "name": "Carbon test cation",
            "composition": {"C": 1},
            "charge": 1,
            "phase_behavior": "permanent_gas",
        },
        {
            "id": "e-",
            "name": "Electron",
            "composition": {},
            "charge": -1,
            "phase_behavior": "electron",
        },
    ]
    payload["reactions"] = [
        {
            "id": "A_ionization",
            "name": "Abstract electron-impact ionization",
            "phase": "gas",
            "stoichiometry": {"A": -1, "A+": 1, "e-": 1},
            "rate_law": {
                "kind": "electron_collision",
                "target_species_id": "A",
                "electron_species_id": "e-",
                "collision_dataset": {
                    "id": "A_ionization_xs",
                    "target_species": "A",
                    "product_species": "A+",
                    "kind": "ionization",
                    "energy_levels_eV": [0, 5, 10, 20, 100],
                    "cross_sections_m2": [0, 0, 1e-24, 2e-24, 1e-24],
                    "source": "abstract apparatus integration fixture",
                    "method": "tabulated cross-section fixture",
                    "basis": "computed",
                    "mechanism_origin": "user_declared",
                },
                "energy_cost": evidence(1, "J/mol"),
            },
        }
    ]
    payload["initial_inventory"] = [
        {
            "species_id": "A",
            "compartment": "gas_loop",
            "amount": {"value": 0.1, "unit": "mol"},
        },
        {
            "species_id": "e-",
            "compartment": "gas_loop",
            "amount": {"value": 1e-12, "unit": "mol"},
        },
    ]
    return payload


def test_cross_section_eedf_channel_is_coupled_to_apparatus_and_conserves_charge() -> None:
    payload = electron_collision_apparatus_payload()
    off = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert off["final_total_mol"]["A+"] == 0
    payload["spark"] = {
        "mode": "cross_section_eedf",
        "voltage": {"value": 20_000, "unit": "V"},
        "current": {"value": 0.001, "unit": "A"},
        "pulse_width": {"value": 0.001, "unit": "s"},
        "frequency": {"value": 100, "unit": "Hz"},
        "coupling_efficiency": 0.5,
        "eedf": collision_rate_payload()["eedf"],
        "source": "abstract apparatus electrical fixture",
        "basis": "hypothetical",
    }
    on = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert on["final_total_mol"]["A+"] > 0
    assert on["material_ledger"]["closed"] is True
    assert on["plasma_model"]["cross_section_channels"][0]["rate_coefficient_m3_s"] > 0
    assert on["plasma_model"]["self_consistent_eedf_solved"] is False


def test_collision_rate_api_rejects_unproven_answer_fields() -> None:
    payload = collision_rate_payload()
    response = TestClient(app).post("/api/prebiotic/plasma/rate-coefficient", json=payload)
    assert response.status_code == 200, response.text
    payload["expected_yield"] = 0.8
    response = TestClient(app).post("/api/prebiotic/plasma/rate-coefficient", json=payload)
    assert response.status_code == 422
    assert "extra_forbidden" in response.text


def boltzmann_payload(field_Td: float = 200.0) -> dict:
    return {
        "solver": "bolos_two_term",
        "collision_set": "lxcat_morgan_2026-09-22",
        "mole_fractions": {"H2": 0.4, "H2O": 0.3, "CH4": 0.2, "NH3": 0.1},
        "gas_temperature_K": 400.0,
        "reduced_electric_field_Td": field_Td,
        "energy_max_eV": 60.0,
        "energy_cells": 400,
        "initial_electron_temperature_eV": 2.0,
        "maximum_iterations": 200,
        "relative_tolerance": 1e-6,
        "acknowledge_two_term_limit": True,
    }


def test_morgan_collision_set_hash_and_all_target_coverage() -> None:
    status = collision_dataset_status()
    assert status["available"] is True
    assert status["sha256"] == "86e256c3ec690cb3ff759529f099959df777f4d0a331f769925e5265c9e8efd6"
    assert status["coverage"]["process_count"] == 46
    assert status["coverage"]["all_required_targets_complete"] is True
    assert set(status["coverage"]["targets"]) == {"H2", "H2O", "CH4", "NH3"}
    assert status["coverage"]["any_unresolved_dissociation_labels"] is True
    assert status["coverage"]["all_required_targets_have_explicit_neutral_fragment_source"] is False
    assert status["coverage"]["targets"]["CH4"]["unresolved_dissociation_channel_count"] == 4
    assert status["coverage"]["targets"]["CH4"]["explicit_neutral_fragment_labels"] == []
    assert {
        target
        for target, item in status["coverage"]["targets"].items()
        if not item["has_explicit_neutral_fragment_source"]
    } == {"H2", "CH4", "NH3"}
    assert status["coverage"]["targets"]["CH4"]["completeness_scope"] == "electron-energy-distribution input only"


def test_hybrid_collision_set_adds_hash_verified_complete_n2_coverage() -> None:
    status = collision_dataset_status(HYBRID_DATASET_ID)
    assert status["available"] is True
    assert status["sha256"] == "c099abe09c174bd0f789bd840dca4861799b085d49b38445201eb24f0188d84e"
    assert status["coverage"]["process_count"] == 71
    assert status["coverage"]["all_required_targets_complete"] is True
    assert set(status["coverage"]["targets"]) == {"H2", "H2O", "CH4", "NH3", "N2"}
    assert status["coverage"]["targets"]["N2"]["process_count"] == 25


def test_reactive_n2_collision_set_adds_hash_verified_tabata_outcome() -> None:
    status = collision_dataset_status(REACTIVE_N2_DATASET_ID)
    assert status["available"] is True
    assert status["sha256"] == "ef8213da28c9ea0db3b8969a4d95ec5dddd4fc1d4042cb03c7f93df6518e9ae1"
    assert status["coverage"]["process_count"] == 72
    assert status["coverage"]["targets"]["N2"]["process_count"] == 26
    assert status["coverage"]["targets"]["N2"]["neutral_dissociation_present"] is True
    assert status["coverage"]["targets"]["N2"]["has_explicit_neutral_fragment_source"] is True
    assert status["coverage"]["targets"]["N2"]["has_resolved_ground_state_fragment_source"] is False


def test_reactive_n2_eedf_reports_reviewed_neutral_dissociation_rate() -> None:
    payload = boltzmann_payload()
    payload["collision_set"] = REACTIVE_N2_DATASET_ID
    payload["mole_fractions"] = {"H2": 0.25, "H2O": 0.25, "CH4": 0.25, "N2": 0.25}
    result = solve_self_consistent_eedf(BoltzmannEEDFRequest.model_validate(payload))
    channel = next(
        item
        for item in result["rate_coefficients"]
        if item["target"] == "N2" and item["kind"] == "dissociation"
    )
    assert result["solver"]["converged"] is True
    assert result["collision_set"]["selected_process_count"] == (
        result["collision_set"]["eedf_process_count"] + 1
    )
    assert result["collision_set"]["post_solve_outcome_process_count"] == 1
    assert channel["product_label"] == "N + N"
    assert channel["product_state_unresolved"] is True
    assert channel["first_nonzero_energy_eV"] == 12.0
    assert channel["rate_coefficient_m3_s"] > 0
    assert "excluded from EEDF" in channel["coupling"]
    support = channel["product_state_energy_support"]
    assert support["measurement_energy_window_eV"] == [18.5, 148.5]
    assert support["measurement_window_extends_above_computed_grid"] is True
    assert support["state_branching_fraction_inferred"] is False
    rate_fractions = [
        support["fraction_of_computed_rate_below_measured_window"],
        support["fraction_of_computed_rate_within_measured_window"],
        support["fraction_of_computed_rate_above_measured_window"],
    ]
    assert all(0 <= value <= 1 for value in rate_fractions)
    assert sum(rate_fractions) == pytest.approx(1.0, abs=1e-10)
    analytic = channel["analytic_fit_diagnostic"]
    assert analytic["rate_coefficient_m3_s"] > 0
    assert analytic["used_for_eedf_or_chemical_source"] is False
    assert analytic["independent_experimental_measurement"] is False
    assert analytic["product_state_energy_support"]["state_branching_fraction_inferred"] is False


def test_qst_n2_analytic_fit_retains_official_energy_domain_and_units() -> None:
    values = _tabata_n2_analytic_fit_m2(np.array([11.99, 12.0, 14.0, 200.01]))
    assert values[0] == 0
    assert values[1] == pytest.approx(7.8560328e-23, rel=1e-6)
    assert values[2] > values[1]
    assert values[3] == 0


def test_barth_initial_n2_product_state_energy_support_is_grid_stable() -> None:
    mixture = {
        "H2": 0.00943, "CO2": 0.0571, "O2": 0.000356,
        "H2O": 0.016, "N2": 0.917114,
    }
    base = {
        "collision_set": IST_CO2_O2_DATASET_ID,
        "mole_fractions": mixture,
        "gas_temperature_K": 298.15,
        "reduced_electric_field_Td": 200,
        "acknowledge_two_term_limit": True,
    }
    channels = []
    for maximum, cells in ((60, 400), (150, 1000)):
        result = solve_self_consistent_eedf(BoltzmannEEDFRequest.model_validate({
            **base, "energy_max_eV": maximum, "energy_cells": cells,
        }))
        assert result["solver"]["converged"] is True
        channels.append(next(
            item for item in result["rate_coefficients"]
            if item["target"] == "N2" and item.get("product_state_unresolved")
        ))
    coarse, extended = channels
    analytic = coarse["analytic_fit_diagnostic"]
    assert 1.08 < analytic["tabular_linear_over_analytic_fit_ratio"] < 1.10
    assert 0.50 < analytic["product_state_energy_support"][
        "fraction_of_computed_rate_below_measured_window"
    ] < 0.53
    assert extended["rate_coefficient_m3_s"] == pytest.approx(
        coarse["rate_coefficient_m3_s"], rel=0.01,
    )
    coarse_support = coarse["product_state_energy_support"]
    extended_support = extended["product_state_energy_support"]
    assert coarse_support["measurement_window_extends_above_computed_grid"] is True
    assert extended_support["measurement_window_extends_above_computed_grid"] is False
    assert extended_support["fraction_of_computed_rate_below_measured_window"] == pytest.approx(
        coarse_support["fraction_of_computed_rate_below_measured_window"], abs=0.001,
    )
    assert 0.5 < extended_support["fraction_of_computed_rate_below_measured_window"] < 0.6
    assert extended_support["state_branching_fraction_inferred"] is False


def test_reactive_n2_dataset_rejects_underresolved_energy_grid() -> None:
    payload = boltzmann_payload()
    payload["collision_set"] = REACTIVE_N2_DATASET_ID
    payload["energy_cells"] = 120
    with pytest.raises(ValueError, match="at least 200 energy cells"):
        BoltzmannEEDFRequest.model_validate(payload)


def test_reactive_n2_dissociation_rate_is_grid_converged() -> None:
    payload = boltzmann_payload()
    payload["collision_set"] = REACTIVE_N2_DATASET_ID
    payload["mole_fractions"] = {"H2": 0.25, "H2O": 0.25, "CH4": 0.25, "N2": 0.25}
    payload["energy_cells"] = 200
    coarse = solve_self_consistent_eedf(BoltzmannEEDFRequest.model_validate(payload))
    payload["energy_cells"] = 400
    fine = solve_self_consistent_eedf(BoltzmannEEDFRequest.model_validate(payload))
    coarse_rate = next(
        item["rate_coefficient_m3_s"]
        for item in coarse["rate_coefficients"]
        if item["target"] == "N2" and item["kind"] == "dissociation"
    )
    fine_rate = next(
        item["rate_coefficient_m3_s"]
        for item in fine["rate_coefficients"]
        if item["target"] == "N2" and item["kind"] == "dissociation"
    )
    assert coarse_rate == pytest.approx(fine_rate, rel=0.01)


def test_hybrid_eedf_solves_n2_precursor_feed_without_product_labels() -> None:
    payload = boltzmann_payload()
    payload["collision_set"] = HYBRID_DATASET_ID
    payload["mole_fractions"] = {"H2": 0.25, "H2O": 0.25, "CH4": 0.25, "N2": 0.25}
    result = solve_self_consistent_eedf(BoltzmannEEDFRequest.model_validate(payload))
    assert result["solver"]["converged"] is True
    assert result["collision_set"]["dataset_id"] == HYBRID_DATASET_ID
    assert sum(item["target"] == "N2" for item in result["rate_coefficients"]) == 24
    assert result["scientific_firewall"]["known_products_accessible_to_solver"] is False


def test_morgan_only_schema_rejects_n2_but_hybrid_accepts_it() -> None:
    payload = boltzmann_payload()
    payload["mole_fractions"] = {"N2": 1.0}
    with pytest.raises(ValueError, match="selected collision set"):
        BoltzmannEEDFRequest.model_validate(payload)
    payload["collision_set"] = HYBRID_DATASET_ID
    assert BoltzmannEEDFRequest.model_validate(payload).mole_fractions == {"N2": 1.0}


def test_self_consistent_eedf_uses_full_mixture_and_converges() -> None:
    request = BoltzmannEEDFRequest.model_validate(boltzmann_payload())
    result = solve_self_consistent_eedf(request)
    assert result["status"] == "self_consistent_eedf_solved"
    assert result["solver"]["engine"] == "BOLOS"
    assert result["solver"]["converged"] is True
    assert result["eedf"]["normalization"] == pytest.approx(1.0, abs=2e-6)
    assert result["eedf"]["mean_energy_eV"] > 0
    assert len(result["rate_coefficients"]) == 42
    assert all(item["rate_coefficient_m3_s"] >= 0 for item in result["rate_coefficients"])
    assert result["scientific_firewall"]["known_products_accessible_to_solver"] is False
    assert result["scientific_firewall"]["quantitative_yield_claim_allowed"] is False
    assert result["independent_reference"]["available_for_this_mixture"] is False


def test_eedf_changes_with_reduced_electric_field() -> None:
    low = solve_self_consistent_eedf(
        BoltzmannEEDFRequest.model_validate(boltzmann_payload(50.0))
    )
    high = solve_self_consistent_eedf(
        BoltzmannEEDFRequest.model_validate(boltzmann_payload(500.0))
    )
    assert high["eedf"]["mean_energy_eV"] > low["eedf"]["mean_energy_eV"]


def test_cantera_provides_independent_ch4_h2o_reference_only() -> None:
    payload = boltzmann_payload()
    payload["mole_fractions"] = {"CH4": 0.5, "H2O": 0.5}
    result = solve_self_consistent_eedf(BoltzmannEEDFRequest.model_validate(payload))
    reference = result["independent_reference"]
    assert reference["available_for_this_mixture"] is True
    assert reference["mean_energy_eV"] > 0
    assert reference["comparison_scope"].startswith("independent trend check")


def test_boltzmann_eedf_api_and_scientific_firewall() -> None:
    client = TestClient(app)
    response = client.post("/api/prebiotic/plasma/eedf/solve", json=boltzmann_payload())
    assert response.status_code == 200, response.text
    assert response.json()["collision_set"]["coverage"]["all_required_targets_complete"] is True
    payload = boltzmann_payload()
    payload["expected_product"] = "glycine"
    response = client.post("/api/prebiotic/plasma/eedf/solve", json=payload)
    assert response.status_code == 422
    assert "extra_forbidden" in response.text
