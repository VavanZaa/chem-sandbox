"""Isolated, assay-blind checks for source-bounded aqueous formaldehyde chemistry."""

import math
from copy import deepcopy

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from app.prebiotic.formaldehyde_hydration import (
    METHANEDIOL_AQ,
    augment_prepared_with_methanediol,
    formaldehyde_contact_hydration_rates_mol_s,
    formaldehyde_hydration_parameters,
    formaldehyde_hydration_enthalpy_J_mol,
    formaldehyde_hydration_rate_mol_s,
    free_formaldehyde_henry_from_effective,
    source_bounded_formaldehyde_contact,
    validate_diagnostic_formaldehyde_range,
)
from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.schema import PrebioticApparatusRequest
from scripts.probe_barth_feed_pathway import make_probe_payload, run_probe


def test_winkelman_rate_and_equilibrium_at_room_temperature() -> None:
    rates = formaldehyde_hydration_parameters(298.15, 5.5)
    assert rates["within_measured_domain"] is True
    assert rates["explicit_extrapolation"] is False
    assert rates["hydration_rate_per_s"] == pytest.approx(10.7885, rel=2e-4)
    assert rates["equilibrium_hydrate_to_free_ratio"] == pytest.approx(1270.655, rel=2e-4)
    assert rates["dehydration_rate_per_s"] == pytest.approx(
        rates["hydration_rate_per_s"] / rates["equilibrium_hydrate_to_free_ratio"]
    )
    free_henry = free_formaldehyde_henry_from_effective(
        3000.0, rates["equilibrium_hydrate_to_free_ratio"]
    )
    assert free_henry == pytest.approx(2.359, rel=2e-3)
    assert free_henry * (1.0 + rates["equilibrium_hydrate_to_free_ratio"]) == pytest.approx(3000.0)
    assert formaldehyde_hydration_enthalpy_J_mol() == pytest.approx(-31336.0, rel=2e-4)


@pytest.mark.parametrize("temperature,pH", [(292.99, 5.5), (333.01, 5.5), (298.15, 4.99), (298.15, 7.01)])
def test_winkelman_rate_fails_closed_outside_evidence_domain(temperature: float, pH: float) -> None:
    with pytest.raises(ValueError, match="measured temperature/pH domain"):
        formaldehyde_hydration_parameters(temperature, pH)
    flagged = formaldehyde_hydration_parameters(temperature, pH, allow_extrapolation=True)
    assert flagged["within_measured_domain"] is False
    assert flagged["explicit_extrapolation"] is True


def test_isolated_hydration_equilibrium_and_exact_reaction_stoichiometry() -> None:
    rates = formaldehyde_hydration_parameters(298.15, 5.5)
    k_h = rates["hydration_rate_per_s"]
    k_d = rates["dehydration_rate_per_s"]
    initial = np.array([2.0e-8, 1.0e-8, 2.7])  # free CH2O, methanediol, liquid H2O

    def rhs(_time: float, state: np.ndarray) -> np.ndarray:
        rate = formaldehyde_hydration_rate_mol_s(*state, k_h, k_d)
        return np.array([-rate, rate, -rate])

    result = solve_ivp(rhs, (0.0, 300.0), initial, method="BDF", rtol=1e-11, atol=1e-22)
    assert result.success
    free, hydrate, water = result.y[:, -1]
    assert hydrate / free == pytest.approx(rates["equilibrium_hydrate_to_free_ratio"], rel=1e-6)
    assert free + hydrate == pytest.approx(initial[0] + initial[1], abs=1e-18)
    assert water + hydrate == pytest.approx(initial[2] + initial[1], abs=1e-13)
    assert free + 2 * hydrate + water == pytest.approx(
        initial[0] + 2 * initial[1] + initial[2], abs=1e-13
    )


def test_reverse_dehydration_returns_water_and_free_formaldehyde() -> None:
    rates = formaldehyde_hydration_parameters(298.15, 5.5)
    reverse = formaldehyde_hydration_rate_mol_s(
        0.0, 1e-6, 2.7,
        rates["hydration_rate_per_s"], rates["dehydration_rate_per_s"],
    )
    assert reverse < 0
    assert math.isfinite(reverse)
    with pytest.raises(ValueError, match="requires liquid water"):
        formaldehyde_hydration_rate_mol_s(1e-6, 0.0, 0.0, 1.0, 0.001)


def test_source_bounded_contact_uses_free_not_effective_henry() -> None:
    contact = source_bounded_formaldehyde_contact(
        temperature_K=298.15, pH=5.5, gas_volume_L=0.95,
        liquid_water_mol=2.7, liquid_side_kla_per_s=1e-4,
    )
    assert contact["effective_total_pool_henry_mol_kg_bar"] == pytest.approx(3000.0)
    assert contact["derived_free_henry_mol_kg_bar"] == pytest.approx(2.359, rel=2e-3)
    assert contact["free_liquid_to_gas_ratio"] == pytest.approx(2.99435, rel=2e-4)
    assert contact["cross_source_derived_not_measured"] is True
    assert contact["kla_is_independently_supplied_not_fitted"] is True
    with pytest.raises(ValueError, match="measured temperature/pH domain"):
        source_bounded_formaldehyde_contact(
            temperature_K=298.15, pH=4.8, gas_volume_L=0.95,
            liquid_water_mol=2.7, liquid_side_kla_per_s=1e-4,
        )


def test_diagnostic_pH_extrapolation_is_bounded_and_separately_flagged() -> None:
    validate_diagnostic_formaldehyde_range(298.15, 4.8)
    rates = formaldehyde_hydration_parameters(298.15, 4.8, allow_extrapolation=True)
    assert rates["within_measured_domain"] is False
    assert rates["explicit_extrapolation"] is True
    for temperature, pH in ((298.15, 3.99), (298.15, 7.01), (292.99, 4.8)):
        with pytest.raises(ValueError, match="bounded"):
            validate_diagnostic_formaldehyde_range(temperature, pH)


def test_distinct_methanediol_identity_has_reaction_stoichiometry_and_no_gas_transfer() -> None:
    prepared = {
        "ids": ["CH2O", "H2O"],
        "index": {"CH2O": 0, "H2O": 1},
        "identities": [
            {"composition": {"C": 1, "H": 2, "O": 1}, "charge": 0},
            {"composition": {"H": 2, "O": 1}, "charge": 0},
        ],
        "transfer": [{}, {}],
        "stoichiometry": np.zeros((2, 0)),
    }
    augment_prepared_with_methanediol(prepared)
    assert prepared["index"][METHANEDIOL_AQ] == 2
    assert prepared["identities"][2]["composition"] == {"C": 1, "H": 4, "O": 2}
    assert prepared["identities"][2]["charge"] == 0
    assert prepared["transfer"][2]["rates_s_inv"] == {
        "evaporation": 0.0, "condensation": 0.0, "wall_loss": 0.0,
    }
    assert prepared["stoichiometry"].shape == (3, 0)
    assert np.array_equal(prepared["atomic"] @ np.array([-1.0, -1.0, 1.0]), np.zeros(4))
    augment_prepared_with_methanediol(prepared)
    assert len(prepared["ids"]) == 3
    invalid = deepcopy(prepared)
    invalid["transfer"][2]["rates_s_inv"]["evaporation"] = 1.0
    with pytest.raises(ValueError, match="reserved"):
        augment_prepared_with_methanediol(invalid)


@pytest.mark.parametrize("initial_carbon", [(1e-6, 0.0, 0.0), (0.0, 0.0, 1e-6)])
def test_isolated_contact_hydration_equilibrium_desorption_and_atom_balance(
    initial_carbon: tuple[float, float, float],
) -> None:
    rates = formaldehyde_hydration_parameters(298.15, 5.5)
    K_h = rates["equilibrium_hydrate_to_free_ratio"]
    H_free = free_formaldehyde_henry_from_effective(3000.0, K_h)
    water_start = 2.7
    pressure_per_gas_mol_bar = 8.314462618 * 298.15 / (0.95e-3) / 1e5
    K_free = H_free * pressure_per_gas_mol_bar * water_start * 0.01801528
    initial = np.array([*initial_carbon, water_start])

    def rhs(_time: float, state: np.ndarray) -> np.ndarray:
        changes = formaldehyde_contact_hydration_rates_mol_s(
            *(max(float(value), 0.0) for value in state),
            free_liquid_to_gas_ratio=K_free,
            liquid_side_kla_per_s=1e-4,  # declared sensitivity, not a measured flask value
            hydration_rate_per_s=rates["hydration_rate_per_s"],
            dehydration_rate_per_s=rates["dehydration_rate_per_s"],
        )
        return np.asarray(changes)

    result = solve_ivp(rhs, (0.0, 100000.0), initial, method="BDF", rtol=1e-10, atol=1e-21)
    assert result.success
    gas, free, hydrate, water = result.y[:, -1]
    assert min(gas, free, hydrate, water) >= -1e-18
    assert free / gas == pytest.approx(K_free, rel=2e-3)
    assert hydrate / free == pytest.approx(K_h, rel=2e-3)
    carbon = np.array([1, 1, 1, 0], dtype=float)
    hydrogen = np.array([2, 2, 4, 2], dtype=float)
    oxygen = np.array([1, 1, 2, 1], dtype=float)
    for atoms in (carbon, hydrogen, oxygen):
        assert float(atoms @ result.y[:, -1]) == pytest.approx(float(atoms @ initial), abs=2e-12)


def _coupled_hydration_payload(*, duration_s: float = 60, substeps: int = 2,
                               remove_carbon: bool = False,
                               model: str = "winkelman_2002_source_bounded") -> dict:
    payload = make_probe_payload(
        duration_s=duration_s, substeps=substeps,
        gas_water_kla_per_s=1e-4, gas_water_solver="coupled_bdf",
        remove_carbon_for_control=remove_carbon,
    )
    payload["apparatus"]["closed_flask_formaldehyde_hydration"] = model
    return payload


def test_coupled_hydration_rejects_missing_pH_or_simultaneous_contact() -> None:
    payload = _coupled_hydration_payload()
    payload["apparatus"]["closed_flask_finite_contact_solver"] = "post_step_analytic"
    with pytest.raises(ValueError, match="requires coupled BDF contact"):
        PrebioticApparatusRequest.model_validate(payload)
    payload["apparatus"]["closed_flask_finite_contact_solver"] = "coupled_bdf"
    payload["apparatus"]["aqueous_acid_base_model"] = "off"
    with pytest.raises(ValueError, match="aqueous acid-base pH"):
        PrebioticApparatusRequest.model_validate(payload)


def test_coupled_apparatus_hydrates_seeded_aqueous_formaldehyde_with_both_ledgers_closed() -> None:
    # A component stress test, not an unseeded HCN-formation or yield run.
    primary_amounts = []
    charged_amounts = []
    for substeps in (2, 4, 8):
        payload = _coupled_hydration_payload(substeps=substeps)
        payload["species"].append({
            "id": "CH2O", "name": "Formaldehyde", "composition": {"C": 1, "H": 2, "O": 1},
            "charge": 0, "smiles": "C=O", "phase_behavior": "condensable",
        })
        payload["initial_inventory"].append({
            "species_id": "CH2O", "compartment": "boiler_liquid",
            "amount": {"value": 1e-6, "unit": "mol"},
        })
        result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
        charged = result["plasma_model"]["propagated_charged_branch"]
        assert result["material_ledger"]["closed"] is True
        assert charged["material_ledger"]["closed"] is True
        assert charged["failed_substep_count"] == 0
        primary_amounts.append(result["final_total_mol"][METHANEDIOL_AQ])
        charged_amounts.append(sum(
            amounts.get(METHANEDIOL_AQ, 0.0)
            for amounts in charged["final_compartments_mol"].values()
        ))
    assert all(9e-7 < amount < 1e-6 for amount in primary_amounts + charged_amounts)
    assert (max(primary_amounts) - min(primary_amounts)) / primary_amounts[-1] < 1e-4
    assert (max(charged_amounts) - min(charged_amounts)) / charged_amounts[-1] < 1e-4
    assert result["scientific_firewall"]["effective_formaldehyde_pool_without_explicit_hydration"] is False
    assert result["scientific_firewall"]["quantitative_yield_claim_allowed"] is False
    heat = result["energy_ledger"]["formaldehyde_hydration_partial_heat"]
    assert heat["not_included_in_temperature_feedback"] is True
    assert heat["net_hydration_extent_mol"] == pytest.approx(result["final_total_mol"][METHANEDIOL_AQ])
    assert heat["heat_release_if_isothermal_J"] == pytest.approx(
        -heat["reaction_enthalpy_J_mol"] * heat["net_hydration_extent_mol"]
    )
    assert result["energy_ledger"]["complete_thermal_energy_balance_closed"] is False


def test_coupled_hydration_carbon_removed_control_stays_carbon_free() -> None:
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(
        _coupled_hydration_payload(remove_carbon=True)
    ))
    charged = result["plasma_model"]["propagated_charged_branch"]
    assert result["material_ledger"]["closed"] is True
    assert charged["material_ledger"]["closed"] is True
    assert charged["failed_substep_count"] == 0
    assert result["final_total_mol"][METHANEDIOL_AQ] == 0
    assert result["final_total_mol"]["HCN"] == 0


def test_coupled_hydration_rejects_long_run_when_pH_leaves_measured_domain() -> None:
    payload = _coupled_hydration_payload(duration_s=1500, substeps=20)
    with pytest.raises(ValueError, match="formaldehyde hydration left measured domain.*pH="):
        simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    diagnostic = _coupled_hydration_payload(
        duration_s=1500, substeps=20,
        model="winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    )
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(diagnostic))
    charged = result["plasma_model"]["propagated_charged_branch"]
    assert result["material_ledger"]["closed"] is True
    assert charged["material_ledger"]["closed"] is True
    assert charged["failed_substep_count"] == 0
    assert 4.0 < result["operating_conditions"]["formaldehyde_hydration_minimum_ideal_pH"] < 5.0
    assert result["scientific_firewall"]["formaldehyde_hydration_pH_extrapolation_used"] is True
    assert result["scientific_firewall"]["quantitative_yield_claim_allowed"] is False


@pytest.mark.parametrize("source_forcing", ("time_averaged", "inventory_limited"))
def test_one_hour_nitrogen_removed_control_has_no_cyanide_without_clipping_oxygen(
    source_forcing: str,
) -> None:
    report = run_probe(
        duration_s=3600, substeps=48, gas_water_kla_per_s=1e-4,
        gas_water_solver="coupled_bdf",
        primary_source_forcing=source_forcing,
        formaldehyde_hydration_model="winkelman_2002_pH4_to_5_extrapolation_diagnostic",
        remove_nitrogen_for_control=True,
    )
    assert report["status"] == "nitrogen_removed_negative_control_not_yield_validation"
    assert report["final_aqueous_HCN_plus_CN_minus_mol"] == 0
    assert report["charged_final_aqueous_HCN_plus_CN_minus_mol"] == 0
    assert report["material_ledger_closed"] is True
    assert report["charged_branch_material_ledger_closed"] is True
    assert report["charged_branch_failed_substeps"] == 0


def test_inventory_limited_wet_source_tracks_actual_extent_within_energy_budget() -> None:
    payload = _coupled_hydration_payload()
    payload["apparatus"]["primary_source_forcing"] = "inventory_limited"
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert result["material_ledger"]["closed"] is True
    charged = result["plasma_model"]["propagated_charged_branch"]
    assert charged["material_ledger"]["closed"] is True
    assert charged["failed_substep_count"] == 0
    steps = result["plasma_model"]["heavy_species_cycle_summaries"]
    assert all(step["source_forcing"] == "inventory_limited" for step in steps)
    assert any(step["source_extents"] for step in steps)
    assert any(
        item["actual_integrated_extent_mol"] > 0
        for step in steps for item in step["source_extents"]
    )
    for step in steps:
        allocation = step["collision_energy_allocation"]
        assert allocation["resolved_fragment_energy_used_J"] <= (
            allocation["allocated_collision_energy_J"] * (1 + 1e-7)
        )
        assert all(
            item["actual_integrated_extent_mol"] >= 0
            for item in step["source_extents"]
        )
