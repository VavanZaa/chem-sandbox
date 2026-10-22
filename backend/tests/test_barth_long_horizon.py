"""Guard the feed-only long-horizon stress runner against assay/fidelity claims."""

import math
import warnings

import cantera as ct
import pytest
from scipy.constants import Avogadro

from app.prebiotic.heavy_gas import (
    HEBRARD_2012_H2CN_N_TO_CH2_N2_RATE_CM3_MOLECULE_S,
    HEBRARD_2012_H2CN_N_TO_HCN_NH_RATE_CM3_MOLECULE_S,
    HEBRARD_2012_H2CN_H_TO_HCN_H2_RATE_CM3_MOLECULE_S,
    HEBRARD_2012_H2CN_H_TO_HNC_H2_RATE_CM3_MOLECULE_S,
    HEBRARD_2012_HNC_H_TO_HCN_H_PREEXP_CM3_MOLECULE_S,
    HEBRARD_2012_HNC_H_TO_HCN_H_BARRIER_K,
    LOW_T_H2CN_H_TO_HCN_RATE_CM3_MOLECULE_S,
    _gri_phase_with_diagnostic_rates,
)
from scripts.probe_barth_feed_pathway import run_probe
from scripts.probe_barth_long_horizon import run_long_horizon


def test_long_horizon_rejects_outside_reported_duration_or_oversized_cycle() -> None:
    with pytest.raises(ValueError, match="no longer than the published spark duration"):
        run_long_horizon(duration_min=926, cycles=16)
    with pytest.raises(ValueError, match="each cycle at or below one hour"):
        run_long_horizon(duration_min=120, cycles=1)
    with pytest.raises(ValueError, match="one altered-feed negative control"):
        run_long_horizon(
            duration_min=1, cycles=1,
            remove_carbon_for_control=True,
            remove_nitrogen_for_control=True,
        )
    with pytest.raises(ValueError, match="finite and positive"):
        run_long_horizon(duration_min=1, cycles=1, spark_coupling_factor=0)
    with pytest.raises(ValueError, match="physical fraction one"):
        run_long_horizon(duration_min=1, cycles=1, spark_coupling_factor=4)
    with pytest.raises(ValueError, match="finite and nonnegative"):
        run_long_horizon(duration_min=1, cycles=1, gri_direct_hcn_rate_factor=-1)
    with pytest.raises(ValueError, match="finite and nonnegative"):
        run_long_horizon(duration_min=1, cycles=1, gri_hcn_h_association_rate_factor=float("nan"))
    with pytest.raises(ValueError, match="finite and nonnegative"):
        run_long_horizon(duration_min=1, cycles=1, charged_hcn_h_association_rate_factor=-1)
    with pytest.raises(ValueError, match="finite and nonnegative"):
        run_long_horizon(duration_min=1, cycles=1, charged_hcnn_o_to_hcn_rate_factor=float("nan"))
    with pytest.raises(ValueError, match="finite and nonnegative"):
        run_long_horizon(duration_min=1, cycles=1, charged_hcnn_h_to_ch2_n2_rate_factor=-1)
    with pytest.raises(ValueError, match="finite and nonnegative"):
        run_long_horizon(duration_min=1, cycles=1, charged_ch_water_rate_factor=-1)
    with pytest.raises(ValueError, match="finite and nonnegative"):
        run_long_horizon(duration_min=1, cycles=1, charged_ch_water_rate_factor=float("nan"))
    with pytest.raises(ValueError, match="finite and nonnegative"):
        run_long_horizon(duration_min=1, cycles=1, charged_ch_n2_rate_factor=-1)


def test_diagnostic_gri_rate_stress_is_isolated_from_default_mechanism() -> None:
    baseline = _gri_phase_with_diagnostic_rates(())
    stressed = _gri_phase_with_diagnostic_rates(((275, 0.0), (236, 1e-5)))
    literature = _gri_phase_with_diagnostic_rates(
        ((275, 0.0),), include_literature_h2cn_h_to_hcn=True,
    )
    branched = _gri_phase_with_diagnostic_rates(
        ((275, 0.0),),
        include_literature_h2cn_h_to_hcn=True,
        include_hebrard_h2cn_n_branches=True,
    )
    isomers = _gri_phase_with_diagnostic_rates(
        (), include_hebrard_h2cn_h_isomer_branches=True,
    )
    isomerized = _gri_phase_with_diagnostic_rates(
        (), include_hebrard_h2cn_h_isomer_branches=True,
        include_hebrard_hnc_h_isomerization=True,
    )
    composition = {"H2": 0.36, "CH4": 0.27, "NH3": 0.27, "H2O": 0.10}
    baseline.TPX = 298.15, ct.one_atm, composition
    stressed.TPX = 298.15, ct.one_atm, composition
    literature.TPX = 298.15, ct.one_atm, composition
    branched.TPX = 298.15, ct.one_atm, composition
    isomers.TPX = 298.15, ct.one_atm, composition
    isomerized.TPX = 298.15, ct.one_atm, composition
    assert stressed.forward_rate_constants[275] == 0
    assert stressed.forward_rate_constants[236] == pytest.approx(
        baseline.forward_rate_constants[236] * 1e-5
    )
    assert baseline.forward_rate_constants[275] > 0
    assert literature.n_reactions == baseline.n_reactions + 1
    assert literature.reaction_equations()[-1] == "H + H2CN <=> H2 + HCN"
    assert literature.forward_rate_constants[-1] * 1000 / Avogadro == pytest.approx(
        LOW_T_H2CN_H_TO_HCN_RATE_CM3_MOLECULE_S
    )
    assert baseline.n_reactions == 325
    assert branched.n_reactions == baseline.n_reactions + 2
    assert branched.reaction_equations()[-1] == "H2CN + N <=> HCN + NH"
    assert branched.forward_rate_constants[237] * 1000 / Avogadro == pytest.approx(
        HEBRARD_2012_H2CN_N_TO_CH2_N2_RATE_CM3_MOLECULE_S
    )
    assert branched.forward_rate_constants[-1] * 1000 / Avogadro == pytest.approx(
        HEBRARD_2012_H2CN_N_TO_HCN_NH_RATE_CM3_MOLECULE_S
    )
    assert "HNC" not in baseline.species_names
    assert "HNC" in isomers.species_names
    assert isomers.n_reactions == baseline.n_reactions + 2
    assert isomers.reaction_equations()[-2:] == [
        "H + H2CN <=> H2 + HCN", "H + H2CN <=> H2 + HNC",
    ]
    assert isomers.forward_rate_constants[-2] * 1000 / Avogadro == pytest.approx(
        HEBRARD_2012_H2CN_H_TO_HCN_H2_RATE_CM3_MOLECULE_S
    )
    assert isomers.forward_rate_constants[-1] * 1000 / Avogadro == pytest.approx(
        HEBRARD_2012_H2CN_H_TO_HNC_H2_RATE_CM3_MOLECULE_S
    )
    assert isomerized.n_reactions == isomers.n_reactions + 1
    assert isomerized.reaction_equations()[-1] == "H + HNC <=> H + HCN"
    assert isomerized.forward_rate_constants[-1] * 1000 / Avogadro == pytest.approx(
        HEBRARD_2012_HNC_H_TO_HCN_H_PREEXP_CM3_MOLECULE_S
        * math.exp(-HEBRARD_2012_HNC_H_TO_HCN_H_BARRIER_K / 298.15)
    )
    assert isomerized.reverse_rate_constants[-1] > 0
    with pytest.raises(ValueError, match="requires the opt-in HNC species branch"):
        _gri_phase_with_diagnostic_rates(
            (), include_hebrard_hnc_h_isomerization=True,
        )
    with pytest.raises(ValueError, match="either Pearce or Hebrard"):
        _gri_phase_with_diagnostic_rates(
            (), include_literature_h2cn_h_to_hcn=True,
            include_hebrard_h2cn_h_isomer_branches=True,
        )
    with pytest.raises(ValueError, match="Unsupported diagnostic GRI reaction index"):
        _gri_phase_with_diagnostic_rates(((0, 0.0),))
    with pytest.raises(ValueError, match="must be unique"):
        _gri_phase_with_diagnostic_rates(((275, 0.0), (275, 1.0)))


def test_current_low_contact_full_horizon_rejects_out_of_domain_hydration() -> None:
    # This used to be reported as pH >= 5 before the N2 fragment-state
    # correction. The current model must not silently reuse that old result.
    with pytest.raises(ValueError, match="hydration left measured domain"):
        run_long_horizon(
            duration_min=925, cycles=16, substeps_per_cycle=12,
            gas_water_kla_per_s=1e-5,
            formaldehyde_hydration_model="winkelman_2002_source_bounded",
        )


def test_one_minute_long_horizon_matches_feed_only_probe_without_assay() -> None:
    long = run_long_horizon(
        duration_min=1, cycles=1, substeps_per_cycle=2,
        trace_hcn_chemistry=True,
    )
    short = run_probe(
        duration_s=60, cycles=1, substeps=2,
        gas_water_kla_per_s=1e-4, gas_water_solver="coupled_bdf",
        primary_source_forcing="inventory_limited",
        formaldehyde_hydration_model="winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    )
    assert long["experimental_hcn_assay_loaded"] is False
    assert long["status"] == "long_horizon_numerical_screen_passed_not_yield_validation"
    assert long["total_integration_steps"] == 2
    assert long["primary_material_ledger_closed"] is True
    assert long["charged_material_ledger_closed"] is True
    assert long["charged_failed_substeps"] == 0
    assert long["charged_hcn_endpoint_rate_trace"] is None
    assert long["charged_hcn_integrated_pathway_trace"] is None
    assert long["charged_hcn_integrated_pathway_summary"] is None
    assert long["charged_hcnn_integrated_pathway_summary"] is None
    assert long["charged_hcnn_net_fate_diagnostic"] is None
    assert long["charged_ch_integrated_pathway_summary"] is None
    assert long["charged_electron_source_summary"] is None
    overlap = long["charged_mechanism_overlap"]
    assert overlap["pavan_hcn_reactions_all_reused_from_gri"] is True
    assert overlap["independent_hcn_kinetic_validation"] is False
    assert overlap["pavan_ncn_species_present"] is False
    assert long["spark_coupling_factor_from_hypothetical_baseline"] == 1
    assert long["gri_diagnostic_reaction_rate_factors"] == {
        "CH3_plus_N_to_HCN_plus_H2": 1.0,
        "H_plus_HCN_to_H2CN": 1.0,
    }
    assert long["gri_rate_stress_is_a_fitted_or_validated_mechanism"] is False
    kinetics = long["unmodified_gri_h_hcn_initial_feed_reference_audit"]
    assert kinetics["gri_over_tsang_estimated_high_pressure_limit"] > 100
    assert kinetics["reviewed_collider_terms_are_complete_mixture_rate"] is False
    assert kinetics["experimental_hcn_assay_loaded"] is False
    ch_n2 = long["unmodified_gri_ch_n2_initial_feed_reference_audit"]
    assert ch_n2["source_temperature_in_range"] is True
    assert 0.8 < ch_n2["gri_over_brownsword_argon_low_pressure_limit"] < 1.2
    assert ch_n2["argon_low_pressure_comparison_validates_wet_feed_effective_rate"] is False
    assert ch_n2["association_rate_validates_hcn_yield"] is False
    assert ch_n2["experimental_hcn_assay_loaded"] is False
    ch_water = long["unmodified_gri_ch_water_initial_feed_reference_audit"]
    assert ch_water["gri_over_blitz_ch_water_rate"] > 2
    assert 0.2 < ch_water["gri_initial_water_over_n2_ch_loss_hazard"] < 0.3
    assert ch_water["gri_CH2O_product_assignment_independently_verified"] is False
    assert ch_water["experimental_hcn_assay_loaded"] is False
    assert long["ch_water_reference_audit_matches_active_charged_rate_profile"] is True
    assert long["ch_n2_reference_audit_matches_active_charged_rate_profile"] is True
    assert long["h_hcn_reference_audit_matches_active_gri_rate_profile"] is True
    assert long["charged_shadow_uses_gri_rate_stress"] is False
    assert long["declared_electrical_coupling_efficiency"] == pytest.approx(0.3)
    assert long["deposited_energy_basis"] == "declared_electrical_estimate"
    assert long["complete_thermal_energy_balance_closed"] is False
    thermal = long["conditional_water_thermal_capacity_screen"]
    assert thermal["feed_card_contains_hcn_assay"] is False
    assert thermal["complete_thermal_energy_balance_closed"] is False
    assert thermal["actual_heat_rejection_inferred_J"] is None
    assert thermal["declared_total_deposited_energy_J"] == pytest.approx(
        long["declared_electrical_estimate_J"]
    )
    assert long["hcn_gas_reaction_trace"]["status"] == "below_numerical_trace_floor"
    assert long["hcn_gas_reaction_trace"]["top_formation_reactions"] == []
    assert long["hcn_gas_reaction_trace"]["method_substep_counts"] == {
        "5-point Gauss quadrature on adaptive BDF intervals": 1,
        "2-point Gauss quadrature on adaptive BDF intervals": 1,
    }
    source = long["source_diagnostic"]
    assert source["thermochemical_floor_status"] == "accepted_source_thermochemical_floor_passed"
    assert source["ordinary_fragment_enthalpy_deficit_count"] == 0
    assert source["n2_parent_event_capacity_deficit_count"] == 0
    assert source["n2_parent_minimum_capacity_to_outcome_ratio"] is None
    assert source["accepted_neutral_extent_mol_by_target"].get("N2", 0.0) == 0
    assert source["accepted_neutral_extent_mol_by_target"]["CO2"] == 0
    assert source["excluded_CO2_to_CO_X_plus_O1D_channel_count"] > 0
    assert source["excluded_CO2_to_CO_X_plus_O1D_raw_opportunities_mol"] > 0
    assert source["excluded_opportunities_are_product_moles"] is False
    assert source["n2d_fate_status_counts"]["conditional_first_collision_fate_resolved"] > 0
    assert source["n2d_source_state_and_products_propagated"] is False
    assert source["n2d_water_conditional_first_collision_fraction_range"][0] > 0.8
    assert source["n2d_total_outcome_uncoupled_opportunities_mol"] > 0
    assert source["n2d_cosby_product_state_measurement_energy_window_eV"] == [18.5, 148.5]
    assert source["n2d_ground_N4S_energetic_necessity_window_eV"] == [0.0, 14.526]
    assert source["n2d_ground_N4S_energetic_window_used_as_fragment_source"] is False
    assert 0.1 < source["n2d_uncoupled_rate_fraction_in_conditional_ground_N4S_window"] < 0.13
    assert 0.1 < source["n2d_analytic_fit_uncoupled_rate_fraction_in_conditional_ground_N4S_window"] < 0.14
    assert source["n2d_electron_energy_window_fractions_are_state_branching"] is False
    assert sum(source["n2d_uncoupled_rate_fraction_by_cosby_product_state_energy_window"].values()) == pytest.approx(1)
    assert source["n2d_analytic_fit_used_to_propagate_chemistry"] is False
    assert source["n2d_analytic_fit_uncoupled_opportunities_mol"] > 0
    assert sum(source["n2d_analytic_fit_uncoupled_rate_fraction_by_cosby_energy_window"].values()) == pytest.approx(1)
    assert source["o1d_fate_status_counts"]["conditional_fate_resolved"] > 0
    assert source["o1d_fate_status_counts"]["unavailable"] > 0
    assert long["primary_aqueous_HCN_plus_CN_minus_mol"] == pytest.approx(
        short["final_aqueous_HCN_plus_CN_minus_mol"], rel=1e-9, abs=1e-25,
    )
    assert long["charged_aqueous_HCN_plus_CN_minus_mol"] == pytest.approx(
        short["charged_final_aqueous_HCN_plus_CN_minus_mol"], rel=1e-9, abs=1e-25,
    )


def test_charged_hcn_endpoint_rate_trace_is_opt_in_and_not_integrated_extent() -> None:
    report = run_long_horizon(
        duration_min=1, cycles=1, substeps_per_cycle=2,
        trace_charged_hcn_final_rates=True,
    )
    snapshots = report["charged_hcn_endpoint_rate_trace"]
    assert len(snapshots) == 2
    assert [item["substep"] for item in snapshots] == [1, 2]
    assert all("not integrated reaction extents" in item["rate_basis"] for item in snapshots)
    assert all(item["hcn"] is not None for item in snapshots)
    assert report["experimental_hcn_assay_loaded"] is False
    assert report["charged_material_ledger_closed"] is True


def test_charged_integrated_hcn_trace_reconciles_without_changing_prediction() -> None:
    baseline = run_long_horizon(duration_min=60, cycles=1, substeps_per_cycle=8)
    traced = run_long_horizon(
        duration_min=60, cycles=1, substeps_per_cycle=8,
        trace_charged_hcn_integrated=True,
    )
    summary = traced["charged_hcn_integrated_pathway_summary"]
    assert summary["status"] == "quadrature_reconciled"
    assert summary["substeps_traced"] == 8
    assert summary["positive_signed_net_reaction_hcn_mol"] > 1e-11
    assert summary["negative_signed_net_reaction_hcn_mol"] > 1e-11
    assert summary["opposing_signed_reaction_contributions_over_abs_observed_net"] > 1
    assert summary["opposing_contribution_ratio_is_yield_uncertainty_multiplier"] is False
    assert abs(summary["total_quadrature_closure_residual_mol"]) < 1e-18
    assert summary["top_formation_reactions"][0]["equation"] == (
        "HCNN + O <=> HCN + NO"
    )
    assert summary["top_consumption_reactions"][0]["equation"] == (
        "H + HCN (+M) <=> H2CN (+M)"
    )
    precursor = traced["charged_hcnn_integrated_pathway_summary"]
    assert precursor["status"] == "quadrature_reconciled"
    assert precursor["top_formation_reactions"][0]["equation"] == (
        "CH + N2 (+M) <=> HCNN (+M)"
    )
    ch = traced["charged_ch_integrated_pathway_summary"]
    assert ch["status"] == "quadrature_reconciled"
    assert ch["top_formation_reactions"][0]["equation"] == "CH + H2 <=> CH2 + H"
    assert ch["top_consumption_reactions"][0]["equation"] == (
        "CH + N2 (+M) <=> HCNN (+M)"
    )
    ch2o = traced["charged_ch2o_integrated_pathway_summary"]
    assert ch2o["status"] == "quadrature_reconciled"
    assert abs(ch2o["total_quadrature_closure_residual_mol"]) < 1e-14
    assert ch2o["top_formation_reactions"][0]["equation"] == (
        "H + HCO (+M) <=> CH2O (+M)"
    )
    assert traced["charged_ch2oh_integrated_pathway_summary"]["top_formation_reactions"][0]["equation"] == (
        "CH2O + H (+M) <=> CH2OH (+M)"
    )
    assert traced["charged_ch3o_integrated_pathway_summary"]["top_formation_reactions"][0]["equation"] == (
        "CH2O + H (+M) <=> CH3O (+M)"
    )
    co = traced["charged_co_integrated_pathway_summary"]
    assert co["status"] == "quadrature_reconciled"
    assert any(
        item["equation"] == "CO2+ + e => CO + O" and item["formation_mol"] > 1e-6
        for item in co["top_formation_reactions"]
    )
    assert traced["charged_hco_integrated_pathway_summary"]["top_formation_reactions"][0]["equation"] == (
        "HCO + M <=> CO + H + M"
    )
    assert any(
        item["diagnostic"]["species"]["CH2O"]["coupled_methanediol_delta_mol"] > 0
        for item in traced["charged_hcn_integrated_pathway_trace"]
    )
    electron_sources = traced["charged_electron_source_summary"]
    assert electron_sources["status"] == "integrated_accepted_electron_source_extents"
    assert electron_sources["is_product_yield"] is False
    assert electron_sources["channels"][0]["target"] == "H2O"
    fate = traced["charged_hcnn_net_fate_diagnostic"]
    assert fate["status"] == "model_internal_signed_net_fate_resolved"
    assert fate["is_forward_branching_probability"] is False
    assert fate["independent_kinetic_validation"] is False
    routes = {item["reaction_index"]: item for item in fate["selected_routes"]}
    assert routes[260]["fraction_of_all_signed_net_hcnn_consumption"] > 0.85
    assert 0.005 < routes[257]["fraction_of_all_signed_net_hcnn_consumption"] < 0.02
    assert routes[256]["signed_net_hcnn_consumption_mol"] > (
        routes[257]["signed_net_hcnn_consumption_mol"]
    )
    assert routes[257]["signed_net_hcnn_consumption_mol"] == pytest.approx(
        summary["top_formation_reactions"][0]["formation_mol"], rel=1e-9,
    )
    refined = run_long_horizon(
        duration_min=60, cycles=1, substeps_per_cycle=16,
        trace_charged_hcn_integrated=True,
    )
    refined_routes = {
        item["reaction_index"]: item
        for item in refined["charged_hcnn_net_fate_diagnostic"]["selected_routes"]
    }
    for index in (256, 257, 260):
        assert refined_routes[index]["fraction_of_all_signed_net_hcnn_consumption"] == pytest.approx(
            routes[index]["fraction_of_all_signed_net_hcnn_consumption"], abs=1e-4,
        )
    assert refined["charged_hcnn_integrated_pathway_summary"]["status"] == "quadrature_reconciled"
    assert refined["experimental_hcn_assay_loaded"] is False
    assert traced["charged_aqueous_HCN_plus_CN_minus_mol"] == pytest.approx(
        baseline["charged_aqueous_HCN_plus_CN_minus_mol"], rel=1e-9, abs=1e-20,
    )
    assert traced["charged_material_ledger_closed"] is True
    assert traced["experimental_hcn_assay_loaded"] is False
    assert traced["charged_aqueous_HCN_plus_CN_minus_mol"] < traced["numerical_reporting_floor_mol"]


def test_charged_rate_stress_is_isolated_and_restored() -> None:
    kwargs = dict(duration_min=60, cycles=1, substeps_per_cycle=8)
    baseline = run_long_horizon(**kwargs)
    stress = run_long_horizon(**kwargs, charged_hcn_h_association_rate_factor=0.5)
    restored = run_long_horizon(**kwargs)
    assert stress["charged_aqueous_HCN_plus_CN_minus_mol"] > (
        baseline["charged_aqueous_HCN_plus_CN_minus_mol"]
    )
    assert stress["primary_aqueous_HCN_plus_CN_minus_mol"] == pytest.approx(
        baseline["primary_aqueous_HCN_plus_CN_minus_mol"], abs=1e-45,
    )
    assert restored["charged_aqueous_HCN_plus_CN_minus_mol"] == pytest.approx(
        baseline["charged_aqueous_HCN_plus_CN_minus_mol"], rel=1e-9, abs=1e-20,
    )
    assert stress["charged_diagnostic_reaction_rate_factors"]["H_plus_HCN_to_H2CN"] == 0.5
    assert stress["charged_rate_stress_is_a_fitted_or_validated_mechanism"] is False
    assert all(item["charged_material_ledger_closed"] for item in (baseline, stress, restored))
    assert all(item["experimental_hcn_assay_loaded"] is False for item in (baseline, stress, restored))


def test_dominant_charged_hcnn_sink_stress_is_reversible_and_assay_free() -> None:
    kwargs = dict(duration_min=60, cycles=1, substeps_per_cycle=8)
    baseline = run_long_horizon(**kwargs)
    stress = run_long_horizon(**kwargs, charged_hcnn_h_to_ch2_n2_rate_factor=0.5)
    restored = run_long_horizon(**kwargs)
    assert stress["charged_aqueous_HCN_plus_CN_minus_mol"] > (
        baseline["charged_aqueous_HCN_plus_CN_minus_mol"]
    )
    assert stress["primary_aqueous_HCN_plus_CN_minus_mol"] == pytest.approx(
        baseline["primary_aqueous_HCN_plus_CN_minus_mol"], abs=1e-45,
    )
    assert restored["charged_aqueous_HCN_plus_CN_minus_mol"] == pytest.approx(
        baseline["charged_aqueous_HCN_plus_CN_minus_mol"], rel=1e-9, abs=1e-20,
    )
    assert stress["charged_diagnostic_reaction_rate_factors"]["H_plus_HCNN_to_CH2_plus_N2"] == 0.5
    assert stress["charged_rate_stress_is_a_fitted_or_validated_mechanism"] is False
    assert all(item["charged_material_ledger_closed"] for item in (baseline, stress, restored))
    assert all(item["experimental_hcn_assay_loaded"] is False for item in (baseline, stress, restored))


def test_charged_ch_water_competition_stress_reallocates_ch_and_restores() -> None:
    from scripts.audit_ch_water_competition import audit_ch_water_competition

    factor = 1 / audit_ch_water_competition()["gri_over_blitz_ch_water_rate"]
    kwargs = dict(
        duration_min=60, cycles=1, substeps_per_cycle=8,
        trace_charged_hcn_integrated=True,
    )
    baseline = run_long_horizon(**kwargs)
    stress = run_long_horizon(**kwargs, charged_ch_water_rate_factor=factor)
    restored = run_long_horizon(duration_min=60, cycles=1, substeps_per_cycle=8)

    def net_ch_loss(report: dict, equation: str) -> float:
        return next(
            item["consumption_mol"]
            for item in report["charged_ch_integrated_pathway_summary"]["signed_reaction_contributions"]
            if item["equation"] == equation
        )

    assert net_ch_loss(stress, "CH + H2O <=> CH2O + H") < net_ch_loss(
        baseline, "CH + H2O <=> CH2O + H",
    )
    assert net_ch_loss(stress, "CH + N2 (+M) <=> HCNN (+M)") > net_ch_loss(
        baseline, "CH + N2 (+M) <=> HCNN (+M)",
    )
    assert stress["charged_diagnostic_reaction_rate_factors"]["CH_plus_H2O_to_CH2O_plus_H"] == pytest.approx(factor)
    assert stress["ch_water_reference_audit_matches_active_charged_rate_profile"] is False
    assert stress["charged_aqueous_HCN_plus_CN_minus_mol"] < stress["numerical_reporting_floor_mol"]
    assert stress["charged_rate_stress_is_a_fitted_or_validated_mechanism"] is False
    assert stress["primary_aqueous_HCN_plus_CN_minus_mol"] == pytest.approx(
        baseline["primary_aqueous_HCN_plus_CN_minus_mol"], abs=1e-45,
    )
    assert restored["charged_aqueous_HCN_plus_CN_minus_mol"] == pytest.approx(
        baseline["charged_aqueous_HCN_plus_CN_minus_mol"], rel=1e-9, abs=1e-20,
    )
    assert all(item["charged_material_ledger_closed"] for item in (baseline, stress, restored))
    assert all(item["experimental_hcn_assay_loaded"] is False for item in (baseline, stress, restored))


def test_charged_ch_n2_association_stress_changes_hcnn_and_restores() -> None:
    kwargs = dict(duration_min=60, cycles=1, substeps_per_cycle=8)
    baseline = run_long_horizon(**kwargs, trace_charged_hcn_integrated=True)
    stress = run_long_horizon(
        **kwargs, trace_charged_hcn_integrated=True,
        charged_ch_n2_rate_factor=0.5,
    )
    restored = run_long_horizon(**kwargs)

    def net_ch_n2(report: dict) -> float:
        return next(
            item["consumption_mol"]
            for item in report["charged_ch_integrated_pathway_summary"]["signed_reaction_contributions"]
            if item["equation"] == "CH + N2 (+M) <=> HCNN (+M)"
        )

    assert net_ch_n2(stress) < net_ch_n2(baseline)
    assert stress["charged_diagnostic_reaction_rate_factors"]["CH_plus_N2_to_HCNN"] == 0.5
    assert stress["ch_n2_reference_audit_matches_active_charged_rate_profile"] is False
    assert stress["charged_rate_stress_is_a_fitted_or_validated_mechanism"] is False
    assert restored["charged_aqueous_HCN_plus_CN_minus_mol"] == pytest.approx(
        baseline["charged_aqueous_HCN_plus_CN_minus_mol"], rel=1e-9, abs=1e-20,
    )
    assert all(item["charged_material_ledger_closed"] for item in (baseline, stress, restored))
    assert all(item["experimental_hcn_assay_loaded"] is False for item in (baseline, stress, restored))


@pytest.mark.parametrize("stressed_rate", (
    "charged_hcn_h_association_rate_factor",
    "charged_hcnn_h_to_ch2_n2_rate_factor",
    "charged_ch_water_rate_factor",
    "charged_ch_n2_rate_factor",
))
def test_charged_rate_stress_keeps_carbon_and_nitrogen_negative_controls(
    stressed_rate: str,
) -> None:
    kwargs = dict(
        duration_min=60, cycles=1, substeps_per_cycle=8,
        **{stressed_rate: 0.5},
    )
    for removed in ("remove_carbon_for_control", "remove_nitrogen_for_control"):
        report = run_long_horizon(**kwargs, **{removed: True})
        assert report["charged_aqueous_HCN_plus_CN_minus_mol"] == 0.0
        assert report["primary_material_ledger_closed"] is True
        assert report["charged_material_ledger_closed"] is True
        assert report["experimental_hcn_assay_loaded"] is False


def test_diagnostic_gri_stress_threads_through_apparatus_only_on_request() -> None:
    report = run_long_horizon(
        duration_min=1, cycles=1, substeps_per_cycle=2,
        gri_direct_hcn_rate_factor=0.0,
    )
    assert report["status"] == "long_horizon_numerical_screen_passed_not_yield_validation"
    assert report["gri_diagnostic_reaction_rate_factors"]["CH3_plus_N_to_HCN_plus_H2"] == 0
    assert report["charged_shadow_uses_gri_rate_stress"] is False
    assert report["primary_material_ledger_closed"] is True
    assert report["charged_material_ledger_closed"] is True


def test_literature_h2cn_h_diagnostic_does_not_mutate_charged_shadow() -> None:
    report = run_long_horizon(
        duration_min=1, cycles=1, substeps_per_cycle=2,
        gri_direct_hcn_rate_factor=0.0,
        include_literature_h2cn_h_to_hcn=True,
        include_hebrard_h2cn_n_branches=True,
    )
    assert report["status"] == "long_horizon_numerical_screen_passed_not_yield_validation"
    assert report["literature_h2cn_h_to_hcn_in_primary_diagnostic"] is True
    assert "Pearce" in report["literature_h2cn_h_to_hcn_rate_source"]
    assert report["hebrard_h2cn_n_branches_in_primary_diagnostic"] is True
    assert "Hebrard" in report["hebrard_h2cn_n_branch_source"]
    assert report["charged_shadow_uses_gri_rate_stress"] is False
    assert report["primary_material_ledger_closed"] is True
    assert report["charged_material_ledger_closed"] is True


def test_hnc_isomer_branch_is_opt_in_and_not_a_yield_claim() -> None:
    report = run_long_horizon(
        duration_min=1, cycles=1, substeps_per_cycle=2,
        include_hebrard_h2cn_h_isomer_branches=True,
        trace_hcn_chemistry=True, trace_h2cn_chemistry=True,
    )
    assert report["status"] == "long_horizon_numerical_screen_passed_not_yield_validation"
    assert report["experimental_hcn_assay_loaded"] is False
    assert report["hnc_gas_water_fate_reviewed"] is False
    assert report["hebrard_h2cn_h_isomer_branches_in_primary_diagnostic"] is True
    assert report["hnc_gas_reaction_trace"]["species"] == "HNC"
    assert report["primary_material_ledger_closed"] is True
    assert report["primary_gas_ideal_energy_species_covered_each_step"] is True
    assert report["complete_thermal_energy_balance_closed"] is False
    assert report["charged_material_ledger_closed"] is True
    with_isomerization = run_long_horizon(
        duration_min=1, cycles=1, substeps_per_cycle=2,
        include_hebrard_h2cn_h_isomer_branches=True,
        include_hebrard_hnc_h_isomerization=True,
    )
    assert with_isomerization["status"] == "long_horizon_numerical_screen_passed_not_yield_validation"
    assert with_isomerization["experimental_hcn_assay_loaded"] is False
    assert with_isomerization["hebrard_hnc_h_isomerization_in_primary_diagnostic"] is True
    assert "uncertainty factor 10" in with_isomerization["hebrard_hnc_h_isomerization_source"]
    assert with_isomerization["primary_material_ledger_closed"] is True
    assert with_isomerization["charged_material_ledger_closed"] is True
    with pytest.raises(ValueError, match="requires the opt-in HNC species branch"):
        run_long_horizon(
            duration_min=1, cycles=1, substeps_per_cycle=2,
            include_hebrard_hnc_h_isomerization=True,
        )
    with pytest.raises(ValueError, match="either Pearce or Hebrard"):
        run_long_horizon(
            duration_min=1, cycles=1, substeps_per_cycle=2,
            include_literature_h2cn_h_to_hcn=True,
            include_hebrard_h2cn_h_isomer_branches=True,
        )


def test_short_carbon_removed_horizon_stays_below_reporting_floor() -> None:
    report = run_long_horizon(
        duration_min=1, cycles=1, substeps_per_cycle=2,
        remove_carbon_for_control=True,
    )
    assert report["status"] == "long_horizon_numerical_screen_passed_not_yield_validation"
    assert report["primary_aqueous_HCN_plus_CN_minus_mol"] <= report["numerical_reporting_floor_mol"]
    assert report["charged_aqueous_HCN_plus_CN_minus_mol"] <= report["numerical_reporting_floor_mol"]
    assert report["unmodified_gri_h_hcn_initial_feed_reference_audit"] is None
    assert report["unmodified_gri_ch_n2_initial_feed_reference_audit"] is None
    assert report["unmodified_gri_ch_water_initial_feed_reference_audit"] is None
    assert report["ch_water_reference_audit_matches_active_charged_rate_profile"] is False
    assert report["ch_n2_reference_audit_matches_active_charged_rate_profile"] is False
    assert report["h_hcn_reference_audit_matches_active_gri_rate_profile"] is False


def test_one_hour_coupled_solver_has_no_runtime_jacobian_warning() -> None:
    # The old adaptive finite-difference Jacobian overflows on this wet run
    # despite closed endpoint ledgers. Keep this a numerical regression, not
    # a reference-product or experimental-yield target.
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        report = run_long_horizon(
            duration_min=60, cycles=1, substeps_per_cycle=48,
            trace_hcn_chemistry=True,
            trace_h2cn_chemistry=True,
        )
    assert report["status"] == "long_horizon_numerical_screen_passed_not_yield_validation"
    assert report["primary_material_ledger_closed"] is True
    assert report["charged_material_ledger_closed"] is True
    trace = report["hcn_gas_reaction_trace"]
    assert trace["status"] == "below_numerical_trace_floor"
    assert trace["substeps_traced"] == 48
    assert trace["top_formation_reactions"] == []
    assert trace["top_consumption_reactions"] == []
    radical_trace = report["h2cn_gas_reaction_trace"]
    assert radical_trace["species"] == "H2CN"
    assert radical_trace["status"] == "below_numerical_trace_floor"
    assert radical_trace["substeps_traced"] == 48
    assert radical_trace["top_formation_reactions"] == []
    assert radical_trace["top_consumption_reactions"] == []


def test_one_hour_literature_h2cn_branch_is_audited_without_yield_claim() -> None:
    report = run_long_horizon(
        duration_min=60, cycles=1, substeps_per_cycle=48,
        gri_direct_hcn_rate_factor=0.0,
        gri_hcn_h_association_rate_factor=1e-5,
        include_literature_h2cn_h_to_hcn=True,
        include_hebrard_h2cn_n_branches=True,
        trace_hcn_chemistry=True,
        trace_h2cn_chemistry=True,
    )
    assert report["status"] == "long_horizon_numerical_screen_passed_not_yield_validation"
    assert report["experimental_hcn_assay_loaded"] is False
    assert report["primary_material_ledger_closed"] is True
    assert report["charged_material_ledger_closed"] is True
    assert report["unmodified_gri_h_hcn_initial_feed_reference_audit"] is not None
    assert report["h_hcn_reference_audit_matches_active_gri_rate_profile"] is False
    hcn = report["hcn_gas_reaction_trace"]
    h2cn = report["h2cn_gas_reaction_trace"]
    assert hcn["status"] == "below_numerical_trace_floor"
    assert hcn["top_formation_reactions"] == []
    assert h2cn["status"] == "below_numerical_trace_floor"
