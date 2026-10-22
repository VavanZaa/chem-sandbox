import json
import math
from copy import deepcopy
from pathlib import Path

import pytest
import numpy as np
from scipy.integrate import solve_ivp

from app.prebiotic.apparatus import (
    _aggregate_integration_steps,
    _apply_coupled_henry_transfer,
    simulate_prebiotic_apparatus,
)
from app.prebiotic.coupled_plasma import _channel_identity
from app.prebiotic.heavy_gas import (
    _conservative_project_nonnegative_neutrals,
    _pavan_source_product_name,
)
from app.prebiotic.hcn_formation import (
    evaluate_aqueous_hcn_exposure,
    evaluate_hcn_formation,
)
from app.prebiotic.precursor_benchmark import _analytical_amount_mol
from app.prebiotic.phase_partition import (
    closed_flask_coupled_henry_coefficients,
    transfer_gas_water_finite_rate,
)
from app.prebiotic.schema import PrebioticApparatusRequest
from scripts.audit_hcn_exposure_sensitivity import make_sensitivity_payload
from scripts.audit_barth_feed_convergence import summarize_convergence
from scripts.audit_barth_transfer_sensitivity import summarize_transfer_sensitivity
from scripts.probe_barth_feed_pathway import FEED_CARD, make_probe_payload, run_probe
from scripts.run_hcn_formation_gate import make_formation_scenario, run_gate, summarize_resolution_sweep


def test_exposure_sensitivity_changes_only_hypothetical_input_not_feed() -> None:
    base = make_formation_scenario(cycles=1, substeps=2)
    original = deepcopy(base)
    zone = make_sensitivity_payload(
        base, name="half_zone", axis="discharge_zone_volume", factor=0.5,
    )
    density = make_sensitivity_payload(
        base, name="double_density", axis="electron_number_density", factor=2.0,
    )
    coupling = make_sensitivity_payload(
        base, name="half_coupling", axis="electrical_coupling_efficiency", factor=0.5,
    )
    assert base == original
    assert all(item["initial_inventory"] == original["initial_inventory"] for item in (zone, density, coupling))
    assert zone["apparatus"]["discharge_zone_volume"]["value"] == pytest.approx(0.025)
    assert density["spark"]["boltzmann"]["electron_number_density"]["value"] == pytest.approx(2e15)
    assert coupling["spark"]["coupling_efficiency"] == pytest.approx(0.15)
    assert all(item["spark"]["basis"] == "hypothetical" for item in (zone, density, coupling))
    with pytest.raises(ValueError, match="exceeds physical fraction"):
        make_sensitivity_payload(
            base, name="invalid", axis="electrical_coupling_efficiency", factor=4.0,
        )


def test_barth_2024_aqueous_measurements_are_held_out_and_unit_explicit() -> None:
    path = Path(__file__).resolve().parents[1] / "data" / "benchmarks" / "barth_2024_aqueous_cyanide_v1.json"
    benchmark = json.loads(path.read_text(encoding="utf-8"))
    records = benchmark["records_with_reported_hcn"]
    assert benchmark["benchmark_id"] == "barth_2024_aqueous_cyanide_v1"
    assert "Held-out evaluation only" in benchmark["access_policy"]
    assert benchmark["field_units"]["HCN"] == "micromolar aqueous assay"
    assert len(records) == len({row["ID"] for row in records}) == 13
    assert all(row["HCN"] > 0 and row["duration"] > 0 for row in records)
    record = next(row for row in records if row["ID"] == "220725001")
    aqueous_mol = record["HCN"] * 1e-6 * benchmark["paper_apparatus"]["liquid_water_volume_L"]
    assert aqueous_mol == pytest.approx(38.85e-6)


def test_lxcat_carbon_oxygen_state_labels_keep_charge_and_excitation() -> None:
    ground_oxygen = _channel_identity({"O": 2}, "excitation", "O(3P)+O(3P)")
    assert ground_oxygen["atoms_and_product_charge_balanced"] is True
    assert all(not product["excited"] for product in ground_oxygen["products"])
    excited_oxygen = _channel_identity({"O": 2}, "excitation", "O(3P)+O(1D)")
    assert excited_oxygen["atoms_and_product_charge_balanced"] is True
    assert [product["excited"] for product in excited_oxygen["products"]] == [False, True]
    co2_dissociation = _channel_identity({"C": 1, "O": 2}, "excitation", "CO(X) + O(1D)")
    assert co2_dissociation["atoms_and_product_charge_balanced"] is True
    assert [product["excited"] for product in co2_dissociation["products"]] == [False, True]
    for target, label in [({"O": 2}, "O^- + O(3P)"), ({"C": 1, "O": 2}, "CO(X) + O(-,X)")]:
        attachment = _channel_identity(target, "attachment", label)
        assert attachment["atoms_and_product_charge_balanced"] is True
        assert sorted(product["charge"] for product in attachment["products"]) == [-1, 0]
    cation = _channel_identity({"C": 1, "O": 2}, "ionization", "CO2(+,X)")
    assert cation["atoms_and_product_charge_balanced"] is True
    assert cation["products"][0]["charge"] == 1
    assert cation["products"][0]["excited"] is False


def test_parenthesized_lxcat_cation_maps_to_ion_not_neutral_pavan_species() -> None:
    supported = {"CO2", "CO2+", "O2", "O2+", "O"}
    for target, label, ion in (
        ({"C": 1, "O": 2}, "CO2(+,X)", "CO2+"),
        ({"O": 2}, "O2(+,X)", "O2+"),
    ):
        source = {"product_label": label}
        product = _channel_identity(target, "ionization", label)["products"][0]
        assert _pavan_source_product_name(source, product, supported) == ion
    unsupported_atomic_ion = _channel_identity({"O": 1}, "ionization", "O(+,X)")["products"][0]
    assert _pavan_source_product_name(
        {"product_label": "O(+,X)"}, unsupported_atomic_ion, supported,
    ) is None


def test_trace_negative_gri_species_projection_conserves_atoms_without_new_products() -> None:
    species = ["CO2", "H2O", "H2", "O2", "CH2(S)"]
    elements = np.asarray([
        [1, 0, 0, 0, 1],
        [0, 2, 2, 0, 2],
        [2, 1, 0, 2, 0],
    ], dtype=float)
    raw = np.asarray([0.002, 2.7, 0.0003, 1e-5, -9.7e-13])
    projected, report = _conservative_project_nonnegative_neutrals(raw, elements, species)
    assert np.all(projected >= 0)
    assert projected[-1] == 0
    assert np.max(np.abs(elements @ (projected - raw))) < 5e-16
    assert report["negative_species_clipped_mol"] == {"CH2(S)": pytest.approx(-9.7e-13)}
    assert report["conservative_projection_total_absolute_mol"] > 0
    with pytest.raises(ValueError, match="negative inventory"):
        _conservative_project_nonnegative_neutrals(
            np.asarray([0.002, 2.7, 0.0003, 1e-5, -1.1e-12]), elements, species,
        )
    with pytest.raises(ValueError, match="projection failed atom conservation"):
        _conservative_project_nonnegative_neutrals(
            np.asarray([0.002, -9.7e-13]), np.asarray([[1, 0], [0, 2]]),
            ["CO2", "H2O"],
        )


def test_feed_only_numerical_screen_requires_both_branches_and_closed_controls() -> None:
    reports = [
        {
            "total_integration_steps": steps,
            "cycles": cycles,
            "integration_substeps_per_cycle": 48,
            "duration_s": 3600,
            "spark_on": True,
            "energy_basis": "declared_electrical_estimate",
            "initial_feed_mol": {"CO2": 0.002, "N2": 0.03},
            "status": "feed_matched_short_pathway_diagnostic_not_yield_validation",
            "material_ledger_closed": True,
            "charged_branch_material_ledger_closed": True,
            "charged_branch_failed_substeps": 0,
            "charged_branch_unprojected_ion_roundoff_clip_mol": 1e-18,
            "primary_negative_species_projection_total_absolute_mol": 1e-13,
            "final_aqueous_HCN_plus_CN_minus_mol": primary,
            "charged_final_aqueous_HCN_plus_CN_minus_mol": charged,
        }
        for steps, cycles, primary, charged in (
            (48, 1, 8.0, 8.2), (96, 2, 10.0, 10.1), (192, 4, 9.7, 9.9),
        )
    ]
    assert summarize_convergence(reports)["status"] == (
        "one_hour_feed_only_numerical_screen_passed_not_yield_validation"
    )
    failed = deepcopy(reports)
    failed[-1]["charged_branch_failed_substeps"] = 1
    assert summarize_convergence(failed)["status"] == "one_hour_feed_only_numerical_screen_failed"
    failed = deepcopy(reports)
    failed[-1]["charged_final_aqueous_HCN_plus_CN_minus_mol"] = 8.0
    assert summarize_convergence(failed)["status"] == "one_hour_feed_only_numerical_screen_failed"
    failed = deepcopy(reports)
    failed[-1]["initial_feed_mol"]["CO2"] *= 2
    assert summarize_convergence(failed)["status"] == "one_hour_feed_only_numerical_screen_failed"
    failed = deepcopy(reports)
    failed[-1]["gas_water_contact_solver"] = "coupled_bdf"
    assert summarize_convergence(failed)["status"] == "one_hour_feed_only_numerical_screen_failed"


def test_short_barth_feed_pathway_probe_holds_out_hcn_and_keeps_charged_closure() -> None:
    feed = json.loads(FEED_CARD.read_text(encoding="utf-8"))
    assert "HCN" not in feed["measured_initial_gas_mole_fractions"]
    assert "measured_hcn" not in feed
    payload = make_probe_payload(duration_s=60)
    assert all(item["species_id"] != "CH4" for item in payload["initial_inventory"])
    assert all(item["species_id"] != "HCN" for item in payload["initial_inventory"])
    report = run_probe(duration_s=60)
    assert report["status"] == "feed_matched_short_pathway_diagnostic_not_yield_validation"
    assert report["material_ledger_closed"] is True
    assert report["charged_branch_material_ledger_closed"] is True
    assert report["charged_branch_failed_substeps"] == 0
    assert report["charged_branch_unprojected_ion_roundoff_clip_mol"] <= 1e-15
    assert report["primary_worst_local_element_residual_mol"] < 1e-12
    assert report["primary_negative_species_projection_total_absolute_mol"] >= 0
    assert report["energy_basis"] == "declared_electrical_estimate"
    assert report["thermal_balance_closed"] is False
    assert report["effective_formaldehyde_pool_without_explicit_hydration"] is True
    assert report["numerical_convergence_validated"] is False
    assert report["excluded_CO2_excited_fragment_channels"]
    fate = report["o1d_first_collision_audit"]
    assert fate["evaluated_substeps"] == 2
    assert fate["status_counts"]["conditional_fate_resolved"] == 1
    assert fate["status_counts"]["unavailable"] == 1
    assert 0 < fate["first_resolved_ground_O_fraction"] < 1
    assert fate["maximum_unreviewed_bath_amount_mol_by_species"]["H2O2"] > 0
    assert fate["products_injected_into_solver"] is False
    assert report["final_aqueous_HCN_plus_CN_minus_mol"] < report["HCN_numerical_reporting_floor_mol"]
    dark_control = run_probe(duration_s=60, spark_on=False)
    assert dark_control["status"] == "feed_matched_short_pathway_diagnostic_not_yield_validation"
    assert dark_control["material_ledger_closed"] is True
    assert dark_control["final_gas_HCN_mol"] == 0
    assert dark_control["final_aqueous_HCN_plus_CN_minus_mol"] == 0


def _closed_flask_phase_test_request() -> dict:
    profile = Path(__file__).resolve().parents[1] / "data" / "benchmarks" / "stribling_miller_measurement_profile.json"
    payload = deepcopy(json.loads(profile.read_text(encoding="utf-8")))
    payload["id"] = "closed_flask_phase_test_not_a_formation_prediction"
    payload["species"].append({
        "id": "CO2", "name": "Carbon dioxide", "composition": {"C": 1, "O": 2},
        "charge": 0, "smiles": "O=C=O", "phase_behavior": "permanent_gas",
    })
    payload["initial_inventory"] = [
        {"species_id": "H2O", "compartment": "boiler_liquid", "amount": {"value": 2.775, "unit": "mol"}},
        {"species_id": "N2", "compartment": "gas_loop", "amount": {"value": 0.039, "unit": "mol"}},
        {"species_id": "CO2", "compartment": "gas_loop", "amount": {"value": 0.0001, "unit": "mol"}},
    ]
    payload["apparatus"].update({
        "cycles": 1, "operation_mode": "closed_stirred_flask",
        "integration_substeps_per_cycle": 2,
        "boiler_liquid_volume": {"value": 0.05, "unit": "L"},
        "gas_loop_volume": {"value": 1, "unit": "L"},
        "discharge_zone_volume": None,
        "condensate_trap_volume": {"value": 0.05, "unit": "L"},
        "gas_circulation_rate": {"value": 0, "unit": "L/s"},
        "boiler_temperature": {"value": 298.15, "unit": "K"},
        "discharge_temperature": {"value": 298.15, "unit": "K"},
        "condenser_temperature": {"value": 298.15, "unit": "K"},
        "boil_duration": {"value": 0, "unit": "s"},
        "discharge_duration": {"value": 3600, "unit": "s"},
        "condense_duration": {"value": 0, "unit": "s"},
        "condenser_phase_transfer_model": "step_equilibrium",
    })
    payload["spark"] = {"mode": "off"}
    return payload


def test_closed_stirred_flask_is_no_condenser_mass_conserving_phase_limit() -> None:
    result = simulate_prebiotic_apparatus(
        PrebioticApparatusRequest.model_validate(_closed_flask_phase_test_request())
    )
    final = result["final_compartments_mol"]
    assert result["material_ledger"]["closed"] is True
    assert result["total_duration_s"] == pytest.approx(3600)
    assert result["operating_conditions"]["operation_mode"] == "closed_stirred_flask"
    assert final["boiler_liquid"]["CO2"] > 0
    assert final["gas_loop"]["CO2"] < 0.0001
    assert all(amount == pytest.approx(0) for amount in final["condensate_trap"].values())
    assert all(amount == pytest.approx(0) for amount in final["collected_condensate"].values())
    assert result["hcn_formation"]["net_hcn_generated_mol"] == pytest.approx(0)
    assert result["hcn_formation"]["aqueous_hydrolysis_screen"]["prediction_inventory_mutated"] is False
    assert result["hcn_formation"]["aqueous_hydrolysis_screen"]["combined_screen_mol"]["nominal_mol"] == 0
    assert result["energy_ledger"]["complete_thermal_energy_balance_closed"] is False


def _finite_contact_request(kla_per_s: float) -> dict:
    payload = _closed_flask_phase_test_request()
    payload["apparatus"]["closed_flask_gas_water_transfer_model"] = "finite_rate_contact"
    payload["apparatus"]["closed_flask_liquid_side_kla"] = {
        "value": kla_per_s, "unit": "1/s", "source": "test-only mass-transfer fixture",
        "method": "two-compartment analytic relaxation check", "basis": "hypothetical",
    }
    return payload


def test_closed_flask_finite_contact_has_zero_and_equilibrium_limits() -> None:
    zero = simulate_prebiotic_apparatus(
        PrebioticApparatusRequest.model_validate(_finite_contact_request(0))
    )
    finite = simulate_prebiotic_apparatus(
        PrebioticApparatusRequest.model_validate(_finite_contact_request(1e-4))
    )
    fast = simulate_prebiotic_apparatus(
        PrebioticApparatusRequest.model_validate(_finite_contact_request(1))
    )
    equilibrium = simulate_prebiotic_apparatus(
        PrebioticApparatusRequest.model_validate(_closed_flask_phase_test_request())
    )
    initial_co2 = 1e-4
    zero_liquid = zero["final_compartments_mol"]["boiler_liquid"]["CO2"]
    finite_liquid = finite["final_compartments_mol"]["boiler_liquid"]["CO2"]
    fast_liquid = fast["final_compartments_mol"]["boiler_liquid"]["CO2"]
    target_liquid = equilibrium["final_compartments_mol"]["boiler_liquid"]["CO2"]
    ratio = target_liquid / (initial_co2 - target_liquid)
    expected = target_liquid * (1 - math.exp(-1e-4 * (1 + ratio) * 3600))
    assert zero_liquid == pytest.approx(0)
    assert finite_liquid == pytest.approx(expected, rel=1e-9)
    assert fast_liquid == pytest.approx(target_liquid, rel=1e-9)
    for result in (zero, finite, fast):
        final = result["final_compartments_mol"]
        assert final["gas_loop"]["CO2"] + final["boiler_liquid"]["CO2"] == pytest.approx(initial_co2)
        assert result["material_ledger"]["closed"] is True
        assert "apparatus:closed_flask_liquid_side_kla" in result["evidence"]["hypothetical_sources"]
        assert result["operating_conditions"]["closed_flask_gas_water_transfer_model"] == "finite_rate_contact"


def test_closed_flask_finite_contact_requires_declared_kla() -> None:
    payload = _finite_contact_request(1e-4)
    del payload["apparatus"]["closed_flask_liquid_side_kla"]
    with pytest.raises(ValueError, match="requires declared liquid-side kLa"):
        PrebioticApparatusRequest.model_validate(payload)


def test_feed_only_transfer_audit_keeps_hcn_assay_out_of_solver() -> None:
    equilibrium = make_probe_payload(duration_s=3600, substeps=48)
    finite = make_probe_payload(duration_s=3600, substeps=48, gas_water_kla_per_s=1e-4)
    assert finite["initial_inventory"] == equilibrium["initial_inventory"]
    assert finite["spark"] == equilibrium["spark"]
    assert "HCN" not in {item["species_id"] for item in finite["initial_inventory"]}
    assert finite["apparatus"]["closed_flask_liquid_side_kla"]["basis"] == "hypothetical"
    assert finite["apparatus"]["closed_flask_liquid_side_kla"]["value"] == 1e-4
    base_report = {
        "status": "feed_matched_short_pathway_diagnostic_not_yield_validation",
        "duration_s": 3600, "spark_on": True,
        "initial_feed_mol": {"N2": 0.03, "CO2": 0.002},
        "energy_basis": "declared_electrical_estimate",
        "total_integration_steps": 48,
        "cycles": 1, "integration_substeps_per_cycle": 48,
        "material_ledger_closed": True,
        "charged_branch_material_ledger_closed": True,
        "charged_branch_failed_substeps": 0,
        "gas_water_contact_model": "finite_rate_contact",
        "gas_water_contact_solver": "post_step_analytic",
        "hypothetical_liquid_side_kla_per_s": 1e-4,
        "final_aqueous_HCN_plus_CN_minus_mol": 1e-11,
        "final_gas_HCN_mol": 2e-12,
        "HCN_numerical_reporting_floor_mol": 1e-12,
    }
    equilibrium_report = deepcopy(base_report)
    equilibrium_report["gas_water_contact_model"] = "step_equilibrium"
    equilibrium_report["hypothetical_liquid_side_kla_per_s"] = None
    reports = [base_report, equilibrium_report]
    summary = summarize_transfer_sensitivity(reports)
    assert summary["material_and_charged_controls_pass"] is True
    assert summary["experimental_hcn_amount_loaded"] is False
    mismatch = deepcopy(reports)
    mismatch[-1]["initial_feed_mol"]["CO2"] *= 2
    assert summarize_transfer_sensitivity(mismatch)["material_and_charged_controls_pass"] is False
    failure = summarize_transfer_sensitivity(reports, [{
        "hypothetical_liquid_side_kla_per_s": 1e-3,
        "error": "negative inventory", "numerical_result_accepted": False,
    }])
    assert failure["material_and_charged_controls_pass"] is False
    assert failure["failed_scenarios"][0]["numerical_result_accepted"] is False


def test_finite_contact_allows_hcn_desorption_without_creating_matter() -> None:
    prepared = {"index": {"H2O": 0, "HCN": 1}}
    gas, liquid, audit = transfer_gas_water_finite_rate(
        np.array([0.0, 0.0]), np.array([2.775, 1e-6]), prepared,
        298.15, 101325.0, 1.0, 3600.0, 1e-4,
    )
    assert gas[1] > 0
    assert 0 < liquid[1] < 1e-6
    assert gas[1] + liquid[1] == pytest.approx(1e-6)
    assert audit["transferred_to_water_mol"]["HCN"] < 0


def test_coupled_closed_contact_rhs_matches_declared_kla_balance() -> None:
    coefficients = closed_flask_coupled_henry_coefficients(
        {"index": {"H2O": 0, "HCN": 1}}, 1.0, 298.15, 2.775, 1e-4,
    )
    entry = coefficients["HCN"]
    gas_mol, liquid_mol = 1e-6, 2e-6
    ratio = entry["boiler_liquid_to_gas_equilibrium_ratio"]
    coupled_rhs = entry["boiler_relaxation_rate_per_chemistry_s"] * (
        gas_mol - (gas_mol + liquid_mol) / (1 + ratio)
    )
    assert coupled_rhs == pytest.approx(1e-4 * (ratio * gas_mol - liquid_mol))
    assert entry["condenser_capture_rate_per_chemistry_s"] == 0
    assert closed_flask_coupled_henry_coefficients(
        {"index": {"H2O": 0, "HCN": 1}}, 1.0, 298.15, 2.775, 0,
    ) == {}


def test_coupled_closed_contact_ode_matches_analytic_transfer() -> None:
    prepared = {"index": {"H2O": 0, "HCN": 1}}
    kla = 1e-4
    coefficients = closed_flask_coupled_henry_coefficients(
        prepared, 1.0, 298.15, 2.775, kla,
    )["HCN"]
    ratio = coefficients["boiler_liquid_to_gas_equilibrium_ratio"]
    rate = coefficients["boiler_relaxation_rate_per_chemistry_s"]
    gas_initial, liquid_initial = 1e-6, 2e-6

    def rhs(_time: float, amounts: np.ndarray) -> list[float]:
        exchange = rate * (amounts[0] - sum(amounts) / (1 + ratio))
        return [-exchange, exchange]

    coupled = solve_ivp(rhs, (0, 120), [gas_initial, liquid_initial],
                        rtol=1e-11, atol=1e-20)
    assert coupled.success is True
    analytic_gas, analytic_liquid, _ = transfer_gas_water_finite_rate(
        np.array([0.0, gas_initial]), np.array([2.775, liquid_initial]),
        prepared, 298.15, 101325.0, 1.0, 120.0, kla,
    )
    assert coupled.y[0, -1] == pytest.approx(analytic_gas[1], rel=1e-7)
    assert coupled.y[1, -1] == pytest.approx(analytic_liquid[1], rel=1e-7)
    assert sum(coupled.y[:, -1]) == pytest.approx(gas_initial + liquid_initial)


def test_coupled_closed_contact_rejects_material_trace_clipping() -> None:
    compartments = {
        "boiler_liquid": np.zeros(1),
        "condensate_trap": np.zeros(1),
    }
    bad_report = {"boiler_delta_mol": {"HCN": -1e-12}, "captured_mol": {}}
    with pytest.raises(ValueError, match="negative inventory"):
        _apply_coupled_henry_transfer(
            compartments, {"HCN": 0}, bad_report,
            max_negative_roundoff_mol=1e-18,
        )
    tiny = {"boiler_delta_mol": {"HCN": -1e-19}, "captured_mol": {}}
    _apply_coupled_henry_transfer(
        {"boiler_liquid": np.zeros(1), "condensate_trap": np.zeros(1)},
        {"HCN": 0}, tiny, max_negative_roundoff_mol=1e-18,
    )
    assert tiny["boiler_negative_roundoff_clipped_mol"]["HCN"] == pytest.approx(1e-19)


def test_coupled_closed_contact_requires_boltzmann_source() -> None:
    payload = _finite_contact_request(1e-4)
    payload["apparatus"]["closed_flask_finite_contact_solver"] = "coupled_bdf"
    with pytest.raises(ValueError, match="requires Boltzmann EEDF mode"):
        PrebioticApparatusRequest.model_validate(payload)


def test_short_feed_only_coupled_contact_closes_both_material_ledgers() -> None:
    report = run_probe(
        duration_s=60, substeps=2, gas_water_kla_per_s=1e-4,
        gas_water_solver="coupled_bdf",
    )
    assert report["gas_water_contact_solver"] == "coupled_bdf"
    assert report["material_ledger_closed"] is True
    assert report["charged_branch_material_ledger_closed"] is True
    assert report["charged_branch_failed_substeps"] == 0
    assert report["final_aqueous_HCN_plus_CN_minus_mol"] >= 0


def test_coupled_contact_carbon_removed_control_forms_no_hcn() -> None:
    report = run_probe(
        duration_s=60, substeps=2, gas_water_kla_per_s=1e-4,
        gas_water_solver="coupled_bdf", remove_carbon_for_control=True,
    )
    assert report["status"] == "carbon_removed_negative_control_not_yield_validation"
    assert "CO2" not in report["initial_feed_mol"]
    assert report["material_ledger_closed"] is True
    assert report["charged_branch_material_ledger_closed"] is True
    assert report["charged_branch_failed_substeps"] == 0
    assert report["final_gas_HCN_mol"] == 0
    assert report["final_aqueous_HCN_plus_CN_minus_mol"] == 0
    assert report["charged_final_aqueous_HCN_plus_CN_minus_mol"] == 0


def test_closed_stirred_flask_rejects_fake_boiler_stage() -> None:
    payload = _closed_flask_phase_test_request()
    payload["apparatus"]["boil_duration"] = {"value": 1, "unit": "s"}
    with pytest.raises(ValueError, match="cannot boil"):
        PrebioticApparatusRequest.model_validate(payload)


def _closed_flask_declared_gas_reaction_request(*, electron_driven: bool) -> dict:
    payload = _closed_flask_phase_test_request()
    payload["species"].append({
        "id": "N2_marker", "name": "Test-only N2 state marker",
        "composition": {"N": 2}, "charge": 0, "phase_behavior": "permanent_gas",
    })
    payload["initial_inventory"] = [
        {"species_id": "N2", "compartment": "gas_loop", "amount": {"value": 0.01, "unit": "mol"}},
    ]
    payload["apparatus"]["discharge_duration"] = {"value": 10, "unit": "s"}
    payload["reactions"] = [{
        "id": "state_marker_conversion", "name": "Abstract balanced gas-kinetics fixture",
        "phase": "gas", "stoichiometry": {"N2": -1, "N2_marker": 1},
        "rate_law": {
            "kind": "electron_impact_reduced" if electron_driven else "mass_action",
            "rate_constant": {
                "value": 0.1, "unit": "1/s", "source": "test-only rate fixture",
                "method": "analytical first-order check", "basis": "hypothetical",
            },
        },
    }]
    if electron_driven:
        payload["apparatus"]["discharge_zone_volume"] = {
            "value": 0.1, "unit": "L", "source": "test-only active-volume fixture",
            "method": "analytical exposure check", "basis": "hypothetical",
        }
        payload["spark"] = {
            "mode": "reduced_order", "voltage": {"value": 1000, "unit": "V"},
            "current": {"value": 0.001, "unit": "A"},
            "pulse_width": {"value": 0.001, "unit": "s"},
            "frequency": {"value": 100, "unit": "Hz"},
            "coupling_efficiency": 0.5,
            "mean_electron_energy": {"value": 5, "unit": "eV"},
            "source": "test-only spark fixture", "basis": "hypothetical",
        }
    return payload


@pytest.mark.parametrize("electron_driven, expected_fraction", [(False, 1.0), (True, 0.1)])
def test_closed_flask_gas_kinetics_do_not_vanish_with_zero_circulation(
    electron_driven: bool, expected_fraction: float,
) -> None:
    payload = _closed_flask_declared_gas_reaction_request(electron_driven=electron_driven)
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert result["numerical_controls"]["discharge_pass_fraction"] == 0
    assert result["numerical_controls"]["declared_electron_channel_exposure_fraction"] == expected_fraction
    effective_rate = 0.1 * (0.01 if electron_driven else 1.0)
    expected_extent = 0.01 * (1 - math.exp(-effective_rate * 10))
    assert result["reaction_extents_mol"]["state_marker_conversion"] == pytest.approx(expected_extent, rel=1e-5)
    assert result["material_ledger"]["closed"] is True


def test_closed_flask_active_spark_requires_declared_zone_volume() -> None:
    payload = _closed_flask_declared_gas_reaction_request(electron_driven=True)
    payload["apparatus"]["discharge_zone_volume"] = None
    with pytest.raises(ValueError, match="requires a declared discharge-zone volume"):
        PrebioticApparatusRequest.model_validate(payload)


def test_closed_flask_boltzmann_refuses_unsupported_initial_carbon_feed() -> None:
    payload = _closed_flask_phase_test_request()
    profile = Path(__file__).resolve().parents[1] / "data" / "benchmarks" / "stribling_miller_measurement_profile.json"
    payload["spark"] = json.loads(profile.read_text(encoding="utf-8"))["spark"]
    payload["apparatus"]["discharge_zone_volume"] = {
        "value": 0.05, "unit": "L", "source": "test-only spark zone",
        "method": "collision-coverage failure check", "basis": "hypothetical",
    }
    with pytest.raises(ValueError, match="no collision-set coverage for initial neutral gas feed: CO2"):
        simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))


def test_closed_stirred_flask_spark_uses_supported_feed_without_condenser() -> None:
    payload = _closed_flask_phase_test_request()
    profile = Path(__file__).resolve().parents[1] / "data" / "benchmarks" / "stribling_miller_measurement_profile.json"
    payload["spark"] = json.loads(profile.read_text(encoding="utf-8"))["spark"]
    payload["species"] = [item for item in payload["species"] if item["id"] != "CO2"]
    payload["initial_inventory"] = [
        {"species_id": "H2O", "compartment": "boiler_liquid", "amount": {"value": 2.775, "unit": "mol"}},
        {"species_id": "H2O", "compartment": "gas_loop", "amount": {"value": 0.0005, "unit": "mol"}},
        {"species_id": "CH4", "compartment": "gas_loop", "amount": {"value": 0.01, "unit": "mol"}},
        {"species_id": "N2", "compartment": "gas_loop", "amount": {"value": 0.025, "unit": "mol"}},
        {"species_id": "H2", "compartment": "gas_loop", "amount": {"value": 0.005, "unit": "mol"}},
    ]
    payload["apparatus"]["discharge_zone_volume"] = {
        "value": 0.05, "unit": "L", "source": "test-only active-zone hypothesis",
        "method": "not measured", "basis": "hypothetical",
    }
    payload["apparatus"]["discharge_duration"] = {"value": 1200, "unit": "s"}
    payload["apparatus"]["aqueous_acid_base_model"] = "hcn_nh3_water_equilibrium"
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert result["material_ledger"]["closed"] is True
    assert result["numerical_controls"]["heavy_gas_chemistry_residence_fraction"] == 1
    assert result["cycles"][0]["heavy_species_snapshot"]["duration_s"] == pytest.approx(600)
    assert result["cycles"][0]["plasma_snapshot"]["discharge_zone_fraction"] == pytest.approx(0.05)
    assert result["plasma_model"]["self_consistent_eedf_solved"] is True
    assert result["plasma_model"]["propagated_charged_branch"]["material_ledger"]["closed"] is True
    assert result["hcn_formation"]["initial_hcn_mol"] == 0
    assert all(amount == pytest.approx(0) for amount in result["final_compartments_mol"]["condensate_trap"].values())
    assert result["scientific_firewall"]["quantitative_yield_claim_allowed"] is False


def test_closed_flask_co2_o2_collision_intake_is_not_a_yield_validation() -> None:
    payload = _closed_flask_phase_test_request()
    profile = Path(__file__).resolve().parents[1] / "data" / "benchmarks" / "stribling_miller_measurement_profile.json"
    payload["species"].append({
        "id": "O2", "name": "Oxygen", "composition": {"O": 2},
        "charge": 0, "smiles": "O=O", "phase_behavior": "permanent_gas",
    })
    payload["initial_inventory"] = [
        {"species_id": "H2O", "compartment": "boiler_liquid", "amount": {"value": 2.775, "unit": "mol"}},
        *[
            {"species_id": species, "compartment": "gas_loop", "amount": {"value": amount, "unit": "mol"}}
            for species, amount in (
                ("N2", 0.03676), ("H2", 0.0004), ("CO2", 0.002),
                ("O2", 0.00004), ("H2O", 0.0008),
            )
        ],
    ]
    payload["spark"] = json.loads(profile.read_text(encoding="utf-8"))["spark"]
    payload["spark"]["boltzmann"]["collision_set"] = (
        "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26"
    )
    payload["apparatus"]["discharge_zone_volume"] = {
        "value": 0.05, "unit": "L", "source": "test-only active-zone hypothesis",
        "method": "not measured", "basis": "hypothetical",
    }
    payload["apparatus"]["discharge_duration"] = {"value": 60, "unit": "s"}
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    assert result["material_ledger"]["closed"] is True
    assert result["plasma_model"]["self_consistent_eedf_solved"] is True
    assert result["scientific_firewall"]["quantitative_yield_claim_allowed"] is False
    assert result["cycles"][0]["plasma_snapshot"]["collision_set"]["dataset_id"] == (
        "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26"
    )
    fate = result["plasma_model"]["cycle_summaries"][0]["o1d_fate_diagnostic"]
    assert fate["status"] == "conditional_fate_resolved"
    assert 0 < fate["physical_ground_O_fraction"] < 1
    assert fate["products_injected_into_solver"] is False
    later_fate = result["plasma_model"]["cycle_summaries"][1]["o1d_fate_diagnostic"]
    assert later_fate["status"] == "unavailable"
    assert "CO" not in later_fate["unreviewed_bath_species"]
    assert "H2O2" in later_fate["unreviewed_bath_species"]
    assert 0 < later_fate["unreviewed_bath_mole_fraction"] < 1e-3
    excluded = result["cycles"][0]["heavy_species_snapshot"][
        "collision_energy_allocation"
    ]["excluded_fragment_channels"]
    assert any(
        item["product_label"] == "CO(X) + O(1D)"
        and item["reason"] == "unresolved_excited_state"
        for item in excluded
    )


def _completed_result(*, initial_hcn: float = 0.0, final_hcn: float = 2e-7) -> dict:
    return {
        "species_identity": [
            {"id": "CH4", "composition": {"C": 1, "H": 4}},
            {"id": "HCN", "composition": {"C": 1, "H": 1, "N": 1}},
        ],
        "material_ledger": {
            "closed": True,
            "initial_mol": {"CH4": 0.02, "HCN": initial_hcn},
            "added_mol": {},
            "sampled_mol": {},
        },
        "final_compartments_mol": {
            "gas_loop": {"HCN": 0.25 * final_hcn},
            "boiler_liquid": {"HCN": 0.5 * final_hcn},
            "collected_condensate": {"HCN": 0.25 * final_hcn},
            "condensate_trap": {},
            "wall": {},
        },
        "energy_ledger": {
            "spark_mode": "boltzmann_eedf",
            "deposited_energy_used_J": 2.0,
            "deposited_energy_basis": "declared_electrical_estimate",
        },
    }


def test_hcn_formation_separates_unseeded_creation_from_fate_and_yield() -> None:
    report = evaluate_hcn_formation(_completed_result())
    assert report["status"] == "unseeded_hcn_formed"
    assert report["net_hcn_generated_mol"] == pytest.approx(2e-7)
    assert report["final_hcn_aqueous_mol"] == pytest.approx(1.5e-7)
    assert report["net_hcn_per_feed_carbon"] == pytest.approx(1e-5)
    assert report["final_hcn_total_nmol_per_declared_J"] == pytest.approx(100.0)
    assert report["net_hcn_generated_nmol_per_declared_J"] == pytest.approx(100.0)
    assert report["experimental_detection_or_yield_validated"] is False
    assert report["scientific_firewall"]["reference_hcn_value_used_for_simulation"] is False
    assert report["scientific_firewall"]["preexisting_gri_mechanism_contains_hcn_reactions"] is True
    assert report["scientific_firewall"]["de_novo_product_identity_discovery_claimed"] is False


def test_hcn_formation_does_not_call_seeded_hcn_a_discovery() -> None:
    report = evaluate_hcn_formation(
        _completed_result(initial_hcn=1e-7, final_hcn=2e-7)
    )
    assert report["status"] == "initial_or_added_hcn_present"
    assert report["net_hcn_generated_mol"] == pytest.approx(1e-7)


def test_hcn_formation_rejects_invalid_fate() -> None:
    result = _completed_result()
    result["final_compartments_mol"]["gas_loop"]["HCN"] = float("nan")
    with pytest.raises(ValueError, match="non-finite"):
        evaluate_hcn_formation(result)


def test_aqueous_hcn_exposure_is_endpoint_quadrature_not_hydrolysis() -> None:
    steps = [
        {"cycle": 1, "substep": 1, "start_time_s": 0.0, "end_time_s": 10.0,
         "compartments_mol": {"boiler_liquid": {
             "HCN": 1e-6, "CN-": 2e-7, "H2O": 10.0}}},
        {"cycle": 1, "substep": 2, "start_time_s": 10.0, "end_time_s": 30.0,
         "compartments_mol": {"boiler_liquid": {
             "HCN": 2e-6, "CN-": 3e-7, "H2O": 10.0}}},
    ]
    original = [step["compartments_mol"]["boiler_liquid"].copy() for step in steps]
    report = evaluate_aqueous_hcn_exposure(
        steps, {"boiler_temperature_K": 373.15, "condenser_temperature_K": 280.0}
    )
    boiler = report["by_compartment"]["boiler_liquid"]
    assert report["step_count"] == 2
    assert report["total_duration_s"] == pytest.approx(30.0)
    assert boiler["molecular_hcn_mol_s"] == pytest.approx(5e-5)
    assert boiler["cyanide_family_mol_s"] == pytest.approx(5.8e-5)
    assert boiler["water_present_s"] == pytest.approx(30.0)
    assert boiler["pH_min"] is not None
    assert report["step_snapshots"][0]["compartments"]["boiler_liquid"]["aqueous_free_formaldehyde_mol"] == 0
    assert report["step_snapshots"][0]["compartments"]["boiler_liquid"]["liquid_volume_L"] > 0
    assert report["aqueous_hydrolysis_calculated"] is False
    assert report["prediction_inventory_mutated"] is False
    assert [step["compartments_mol"]["boiler_liquid"] for step in steps] == original


def test_aqueous_hcn_exposure_rejects_missing_time_interval() -> None:
    with pytest.raises(ValueError, match="contiguous"):
        evaluate_aqueous_hcn_exposure(
            [{"cycle": 1, "start_time_s": 1.0, "end_time_s": 2.0,
              "compartments_mol": {}}],
            {"boiler_temperature_K": 373.15, "condenser_temperature_K": 280.0},
        )


def test_hcn_formation_counts_ionized_cyanide_without_calling_it_molecular_hcn() -> None:
    completed = _completed_result(final_hcn=1e-7)
    completed["final_compartments_mol"]["boiler_liquid"]["CN-"] = 1e-8
    report = evaluate_hcn_formation(completed)
    assert report["net_hcn_generated_mol"] == pytest.approx(1.1e-7)
    assert report["final_hcn_total_mol"] == pytest.approx(1e-7)
    assert report["final_cyanide_ion_total_mol"] == pytest.approx(1e-8)
    assert report["final_hcn_plus_cyanide_total_nmol_per_declared_J"] == pytest.approx(55.0)
    assert report["final_hcn_plus_cyanide_aqueous_nmol_per_declared_J"] == pytest.approx(42.5)
    assert report["net_hcn_generated_basis"] == "HCN plus CN- conserved cyanide family"


def test_cyanide_electrode_comparator_uses_free_cyanide_family() -> None:
    inventory = {"HCN": 2e-7, "CN-": 3e-8, "CH2O": 4e-7}
    assert _analytical_amount_mol("HCN", inventory) == pytest.approx(2.3e-7)
    assert _analytical_amount_mol("CH2O", inventory) == pytest.approx(4e-7)


def test_hcn_resolution_summary_is_relative_to_fine_result_and_not_a_yield_claim() -> None:
    def report(step: int, total: float, aqueous: float) -> dict:
        return {
            "status": "unseeded_hcn_formation_gate_passed",
            "cycles": 1,
            "substeps_per_cycle": step,
            "cases": {"spark_on": {
                "net_hcn_generated_mol": total,
                "final_hcn_aqueous_mol": aqueous,
                "final_hcn_plus_cyanide_aqueous_mol": aqueous,
                "charged_branch_final_hcn_plus_cyanide_total_mol": total * 0.8,
                "charged_branch_final_hcn_plus_cyanide_aqueous_mol": aqueous * 0.8,
            }},
        }

    result = summarize_resolution_sweep([
        report(2, 2.2e-7, 1.1e-7),
        report(4, 2e-7, 1e-7),
    ])
    comparison = result["adjacent_refinements"][0]["spark_on"]
    assert comparison["net_hcn_generated_mol"]["relative_change_from_fine"] == pytest.approx(0.1)
    assert comparison["final_hcn_aqueous_mol"]["relative_change_from_fine"] == pytest.approx(0.1)
    assert comparison["final_hcn_plus_cyanide_aqueous_mol"]["relative_change_from_fine"] == pytest.approx(0.1)
    assert comparison["charged_branch_final_hcn_plus_cyanide_total_mol"]["relative_change_from_fine"] == pytest.approx(0.1)
    assert result["all_input_controls_passed"] is True
    assert result["scientific_firewall"]["full_long_duration_convergence_established"] is False
    assert result["scientific_firewall"]["experimental_magnitude_validated"] is False
    with pytest.raises(ValueError, match="increase strictly"):
        summarize_resolution_sweep([report(4, 1, 1), report(2, 1, 1)])
    incomplete = report(4, 2e-7, 1e-7)
    del incomplete["cases"]["spark_on"]["charged_branch_final_hcn_plus_cyanide_total_mol"]
    with pytest.raises(ValueError, match="missing charged_branch_final_hcn_plus_cyanide_total_mol"):
        summarize_resolution_sweep([report(2, 2.2e-7, 1.1e-7), incomplete])


@pytest.mark.skipif(
    not (
        Path(__file__).resolve().parents[1]
        / "data/collision/private/IST_Lisbon_H2_2026-09-23/Cross section.txt"
    ).exists(),
    reason="Opt-in LXCat contributor files are private and not distributed",
)
def test_hcn_gate_fails_closed_when_n2_product_states_are_unresolved() -> None:
    report = run_gate(cycles=1, substeps=2)
    assert report["status"] == "formation_gate_not_passed"
    assert report["cases"]["spark_on"]["status"] == "no_hcn_above_numerical_floor"
    assert report["cases"]["spark_on"]["net_hcn_generated_mol"] <= report["cases"]["spark_on"]["numerical_reporting_floor_mol"]
    assert (
        report["cases"]["spark_on"]["charged_branch_final_hcn_plus_cyanide_total_mol"]
        > report["cases"]["spark_on"]["numerical_reporting_floor_mol"]
    )
    assert report["cases"]["spark_on"]["charged_to_primary_final_cyanide_ratio"] is None
    assert report["cases"]["spark_on"]["final_cyanide_ion_total_mol"] >= 0
    carbonate = report["cases"]["spark_on"]["final_liquid_carbonate_family_mol"]
    assert carbonate["CO2"] > 0
    assert carbonate["HCO3-"] > 0
    assert carbonate["CO3--"] >= 0
    ammonia = report["cases"]["spark_on"]["final_liquid_ammonia_family_mol"]
    assert ammonia["NH3"] > 0
    assert ammonia["NH4+"] > 0
    assert math.isfinite(
        report["cases"]["spark_on"]["final_liquid_pH_ideal_aqueous_model"]["boiler_liquid"]
    )
    assert report["cases"]["spark_on"]["final_hcn_plus_cyanide_aqueous_mol"] >= report["cases"]["spark_on"]["final_hcn_aqueous_mol"]
    for control in ("spark_off", "no_methane", "no_nitrogen"):
        assert report["cases"][control]["status"] == "no_hcn_above_numerical_floor"
    for control in ("no_methane", "no_nitrogen"):
        assert (
            report["cases"][control]["charged_branch_final_hcn_plus_cyanide_total_mol"]
            <= report["cases"][control]["numerical_reporting_floor_mol"]
        )
        assert report["cases"][control]["charged_to_primary_final_cyanide_ratio"] is None
    assert all(case["material_ledger_closed"] for case in report["cases"].values())
    assert all(
        case["material_ledger_maximum_absolute_residual_mol"] < 1e-12
        for case in report["cases"].values()
    )
    assert all(
        case["charged_branch_material_ledger_closed"] is True
        and case["charged_branch_material_ledger_maximum_absolute_residual_mol"] < 1e-11
        and case["charged_branch_failed_substep_count"] == 0
        for name, case in report["cases"].items() if name != "spark_off"
    )
    assert 0 < report["methane_response"]["half_to_full_hcn_ratio"] < 1
    assert 0 < report["methane_response"]["charged_half_to_full_cyanide_ratio"] < 1
    assert report["scientific_firewall"]["charged_branch_negative_controls_below_floor"] is True
    assert report["scientific_firewall"]["modern_experiment_magnitude_validated"] is False
    assert report["cases"]["spark_on"]["complete_thermal_energy_balance_closed"] is False
    assert report["cases"]["spark_on"]["gas_stage_thermochemistry_coverage_complete"] is True
    assert math.isfinite(
        report["cases"]["spark_on"]["gas_stage_ideal_gas_internal_energy_change_J"]
    )


def test_selected_spark_on_run_cannot_masquerade_as_a_control_gate() -> None:
    report = run_gate(cycles=1, substeps=2, case_names=("spark_on",))
    assert report["status"] == "partial_formation_diagnostic_only"
    assert report["executed_cases"] == ["spark_on"]
    assert report["scientific_firewall"]["all_negative_controls_executed"] is False
    assert report["methane_response"]["half_to_full_hcn_ratio"] is None


def test_control_only_subset_cannot_masquerade_as_a_formation_gate() -> None:
    report = run_gate(cycles=1, substeps=2, case_names=("spark_off",))
    assert report["status"] == "partial_formation_diagnostic_only"
    assert report["executed_cases"] == ["spark_off"]
    assert report["scientific_firewall"]["spark_on_case_executed"] is False
    assert report["scientific_firewall"]["all_negative_controls_executed"] is False
    assert report["methane_response"]["half_to_full_hcn_ratio"] is None


@pytest.mark.skipif(
    not (
        Path(__file__).resolve().parents[1]
        / "data/collision/private/IST_Lisbon_H2_2026-09-23/Cross section.txt"
    ).exists(),
    reason="Opt-in LXCat contributor files are private and not distributed",
)
def test_nitrogen_free_charged_control_remains_exactly_nitrogen_free_across_cycles() -> None:
    report = run_gate(cycles=3, substeps=2, case_names=("spark_on", "no_nitrogen"))
    control = report["cases"]["no_nitrogen"]
    assert control["charged_branch_failed_substep_count"] == 0
    assert control["charged_branch_material_ledger_closed"] is True
    assert control["charged_branch_material_ledger_residual_mol"]["N"] == 0.0
    assert control["net_hcn_generated_mol"] == 0.0


def test_selected_cases_require_a_nonempty_unique_known_set() -> None:
    for selected in ((), ("spark_on", "spark_on"), ("spark_on", "unknown")):
        with pytest.raises(ValueError, match="Selected cases"):
            run_gate(cycles=1, substeps=2, case_names=selected)


def test_aggregation_keeps_earlier_charged_substep_failure() -> None:
    cycles = _aggregate_integration_steps([
        {"start_time_s": 0.0, "end_time_s": 1.0,
         "solver_status": "propagated_pavan_charged_shadow_failed", "solver_error": "synthetic failure"},
        {"start_time_s": 1.0, "end_time_s": 2.0,
         "solver_status": "correlated_pavan_charge_balanced_shadow_complete"},
    ], 2)
    assert cycles[0]["solver_status"] == "correlated_pavan_charge_balanced_shadow_complete"
    assert cycles[0]["solver_failed_in_any_substep"] is True
    assert cycles[0]["integration_solver_statuses"][0].endswith("_failed")
    assert cycles[0]["integration_solver_errors"] == ["synthetic failure", None]
