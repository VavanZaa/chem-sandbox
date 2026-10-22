"""The total N2 dissociation rate never becomes an invented product branch."""

import pytest
from scipy.constants import R

from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.n2d_fate import n2d_fate_diagnostic
from app.prebiotic.schema import PrebioticApparatusRequest
from scripts.audit_n2_source_energy_support import audit_barth_n2_source_energy_support
from scripts.audit_n2d_feed_fate import audit_barth_initial_n2d_fate
from scripts.probe_barth_feed_pathway import make_probe_payload


def _channel(report: dict, species: str) -> dict:
    return next(item for item in report["channels"] if item["bath_species"] == species)


def test_dry_one_percent_methane_room_temperature_n2d_limit_is_conditional() -> None:
    total_mol = 101_325.0 * 1e-3 / (R * 298.15)
    fate = n2d_fate_diagnostic(
        {"CH4": 0.01 * total_mol, "N2": 0.99 * total_mol},
        gas_volume_L=1.0,
        gas_temperature_K=298.15,
    )
    assert fate["status"] == "conditional_first_collision_fate_resolved"
    assert _channel(fate, "CH4")["conditional_first_collision_fraction"] == pytest.approx(
        0.7038535985
    )
    assert _channel(fate, "CH4")["rate_factor_extreme_fraction_range"][0] > 0.5
    assert fate["conditional_thermalized_lifetime_s"] == pytest.approx(7.1486466e-7)
    assert fate["products_injected_into_solver"] is False
    assert fate["product_branching_resolved"] is False
    hot_fate = n2d_fate_diagnostic(
        {"CH4": 0.01 * total_mol, "N2": 0.99 * total_mol},
        gas_volume_L=1.0,
        gas_temperature_K=550.0,
    )
    assert hot_fate["status"] == "unavailable"
    assert "channels" not in hot_fate


def test_barth_initial_n2d_audit_is_water_dominated_but_not_a_yield() -> None:
    report = audit_barth_initial_n2d_fate()
    fate = report["n2d_fate"]
    assert fate["status"] == "conditional_first_collision_fate_resolved"
    assert report["feed_card_contains_hcn_assay"] is False
    assert report["electron_dissociation_state_branch_assigned"] is False
    assert fate["products_injected_into_solver"] is False
    assert fate["product_branching_resolved"] is False
    assert fate["quantitative_hcn_yield_claim_allowed"] is False
    assert 0.85 < _channel(fate, "H2O")["conditional_first_collision_fraction"] < 0.95
    assert fate["conditional_thermalized_lifetime_s"] < 1e-6
    assert sum(row["conditional_first_collision_fraction"] for row in fate["channels"]) == pytest.approx(1)
    lower, central, upper = (
        _channel(fate, "H2O")["rate_factor_extreme_fraction_range"][0],
        _channel(fate, "H2O")["conditional_first_collision_fraction"],
        _channel(fate, "H2O")["rate_factor_extreme_fraction_range"][1],
    )
    assert 0 < lower < central < upper < 1
    assert fate["rate_factor_extremes_are_confidence_intervals"] is False


def test_barth_n2_source_energy_support_scan_is_not_product_branching() -> None:
    report = audit_barth_n2_source_energy_support()
    assert report["feed_card_contains_hcn_assay"] is False
    assert report["source_state_branch_inferred"] is False
    assert report["conditional_ground_atom_energetic_window_used_as_source"] is False
    assert report["hcn_yield_predicted"] is False
    fields = report["field_sweep"]
    assert [row["reduced_electric_field_Td"] for row in fields] == [50, 100, 150, 200]
    assert all(row["eedf_converged"] is True for row in fields)
    assert all(row["product_state_branching_resolved"] is False for row in fields)
    assert all(
        row["tabata_analytic_fit_diagnostic"]["used_for_eedf_or_chemical_source"] is False
        for row in fields
    )
    rates = [row["total_n2_neutral_dissociation_rate_coefficient_m3_s"] for row in fields]
    assert rates == sorted(rates)
    below = [
        row["product_state_energy_support"]["fraction_of_computed_rate_below_measured_window"]
        for row in fields
    ]
    assert below == sorted(below, reverse=True)
    assert below[0] > 0.99 and 0.5 < below[-1] < 0.6
    ground_window = [
        row["ground_atom_energetic_window"]["rate_fraction_below_that_threshold"]
        for row in fields
    ]
    assert ground_window == sorted(ground_window, reverse=True)
    assert 0.8 < ground_window[0] < 0.9
    assert 0.1 < ground_window[-1] < 0.13
    assert all(
        row["ground_atom_energetic_window"]["first_neutral_channel_without_ground_N4S_threshold_eV"]
        == 14.526
        and row["ground_atom_energetic_window"]["ground_atom_or_partner_injected_into_chemistry"]
        is False
        for row in fields
    )
    fit_at_declared_field = fields[-1]["tabata_analytic_fit_diagnostic"]
    assert 1.08 < fit_at_declared_field["tabular_linear_over_analytic_fit_ratio"] < 1.10
    assert 0.50 < fit_at_declared_field["product_state_energy_support"][
        "fraction_of_computed_rate_below_measured_window"
    ] < 0.53
    grids = report["grid_sweep_at_declared_hypothetical_field"]
    assert grids[-1]["product_state_energy_support"][
        "measurement_window_extends_above_computed_grid"
    ] is False
    assert grids[-1]["total_n2_neutral_dissociation_rate_coefficient_m3_s"] == pytest.approx(
        grids[0]["total_n2_neutral_dissociation_rate_coefficient_m3_s"], rel=0.01,
    )
    assert grids[-1]["ground_atom_energetic_window"]["rate_fraction_below_that_threshold"] == pytest.approx(
        grids[0]["ground_atom_energetic_window"]["rate_fraction_below_that_threshold"], abs=1e-3,
    )


def test_n2d_audit_rejects_unreviewed_bath_without_renormalizing() -> None:
    report = n2d_fate_diagnostic(
        {"H2O": 0.01, "unknown": 0.01},
        gas_volume_L=1.0, gas_temperature_K=298.15,
    )
    assert report["status"] == "unavailable"
    assert report["branch_fractions_reported"] is False
    assert report["unreviewed_bath_species"] == ["unknown"]
    assert report["unreviewed_rate_for_equal_reviewed_hazard_cm3_molecule_s"]["unknown"] > 0
    assert report["equal_hazard_rate_is_measured_or_bounded"] is False
    assert "channels" not in report


def test_evaluated_atomic_oxygen_removal_does_not_assign_oxygen_product_state() -> None:
    report = n2d_fate_diagnostic(
        {"H2O": 0.01, "O": 1e-5},
        gas_volume_L=1.0, gas_temperature_K=298.15,
    )
    assert report["status"] == "conditional_first_collision_fate_resolved"
    oxygen = _channel(report, "O")
    assert oxygen["rate_coefficient_cm3_molecule_s"] == pytest.approx(1.4e-12)
    assert oxygen["multiplicative_rate_error_factor_298K"] == 2.0
    assert oxygen["kind"] == "oxygen_product_state_unresolved"
    assert report["product_branching_resolved"] is False
    assert report["products_injected_into_solver"] is False


def test_n2d_audit_is_source_bounded_and_validates_inputs() -> None:
    unavailable = n2d_fate_diagnostic(
        {"H2O": 0.01}, gas_volume_L=1.0, gas_temperature_K=400.0,
    )
    assert unavailable["status"] == "unavailable"
    assert "outside reviewed" in unavailable["reason"]
    with pytest.raises(ValueError, match="positive finite gas volume"):
        n2d_fate_diagnostic({"H2O": 0.01}, gas_volume_L=0, gas_temperature_K=298.15)
    with pytest.raises(ValueError, match="finite and nonnegative"):
        n2d_fate_diagnostic({"H2O": -0.01}, gas_volume_L=1, gas_temperature_K=298.15)


def test_plasma_snapshot_reports_n2d_fate_without_injecting_unspecified_state() -> None:
    payload = make_probe_payload(
        duration_s=60, substeps=2, gas_water_kla_per_s=1e-4,
        gas_water_solver="coupled_bdf", primary_source_forcing="inventory_limited",
        formaldehyde_hydration_model="winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    )
    result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
    snapshot = result["plasma_model"]["cycle_summaries"][0]
    fate = snapshot["n2d_fate_diagnostic"]
    assert fate["status"] == "conditional_first_collision_fate_resolved"
    assert fate["products_injected_into_solver"] is False
    assert snapshot["source_terms"]
    assert any(
        source["target"] == "N2" and source.get("product_state_unresolved")
        for source in snapshot["source_terms"]
    )
    assert all(
        extent["target"] != "N2"
        for extent in result["plasma_model"]["heavy_species_cycle_summaries"][0]["source_extents"]
    )
    assert result["material_ledger"]["closed"] is True
