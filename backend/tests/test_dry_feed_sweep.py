"""The measured feed trend is held out of the dry-flow solver."""

import pytest

from evaluation.reference_only.dry_hcn_2024 import evaluate_frozen_dry_feed_sweep
from scripts.audit_dry_feed_sweep import METHANE_FEED_FRACTIONS, audit_dry_feed_sweep


def test_feed_sweep_conserves_and_fails_hcn_ordering_but_matches_ethane_trend() -> None:
    sweep = audit_dry_feed_sweep()
    assert sweep["measured_product_assay_loaded"] is False
    assert sweep["n2_excitation_opportunities_are_product_yields"] is False
    assert sweep["selectivity_basis_matches_paper_exactly"] is False
    assert [row["feed_CH4_mole_fraction"] for row in sweep["rows"]] == list(
        METHANE_FEED_FRACTIONS
    )
    for row in sweep["rows"]:
        assert row["ch4_conversion_fraction"] == pytest.approx(
            row["ch4_electron_fragment_source_fraction_of_feed"]
            + row["ch4_net_neutral_chemistry_consumption_fraction_of_feed"],
            abs=1e-10,
        )
        assert row["ch4_net_neutral_chemistry_consumption_fraction_of_feed"] < 0
        assert row["electron_density_was_power_capped"] is True
        assert row["material_ledger_closed"] is True
        assert row["collision_energy_upper_bound_closed"] is True
        assert row["n2_product_state_source_excluded"] is True
        assert row["uncoupled_n2_a3_excitation_opportunities_per_feed_CH4"] > 0
        assert row["uncoupled_n2_aprime_excitation_opportunities_per_feed_CH4"] > 0
        assert row["hcn_amount_below_1e-12_mol_reporting_floor"] is True
        assert row["hcn_carbon_yield_fraction_below_reporting_floor"] < row[
            "c2h6_carbon_yield_fraction"
        ]
    report = evaluate_frozen_dry_feed_sweep(sweep)
    assert report["feed_points_match_reference"] is True
    assert report["directional_trend_matches"] is True
    assert report["h2_directional_trend_matches"] is True
    assert report["ch4_conversion_directional_trend_matches"] is False
    assert report["ch4_conversion_low_feed_coarse_plot_bound_matches"] is False
    assert report["ch4_conversion_high_feed_coarse_plot_bound_matches"] is True
    assert report["ch4_conversion_coarse_plot_bounds_both_match"] is False
    assert report["ch4_conversion_plot_bounds_are_statistical_uncertainty"] is False
    assert report["numeric_experimental_conversions_digitized"] is False
    assert report["status"] == "one_or_more_directional_trends_mismatch_under_unmatched_scenario"
    assert [row["ch4_electron_fragment_source_fraction_of_feed"] for row in sweep["rows"]] == pytest.approx(
        [0.2139856823, 0.2133890971, 0.2085080596], rel=1e-5
    )
    assert [row["ch4_net_neutral_chemistry_consumption_fraction_of_feed"] for row in sweep["rows"]] == pytest.approx(
        [-0.0425968335, -0.0410738902, -0.0332654602], rel=1e-5
    )
    assert all(row["resolved_fragment_J"] < row["declared_discharge_J"] for row in sweep["rows"])
    assert [row["uncoupled_n2_a3_excitation_opportunities_per_feed_CH4"] for row in sweep["rows"]] == pytest.approx(
        [22.42084978, 2.21981027, 0.29628825], rel=1e-5
    )
    assert report["like_for_like_quantitative_comparison_ready"] is False
    assert report["numeric_experimental_selectivities_digitized"] is False
    assert report["reference_used_by_solver"] is False
    assert len(report["frozen_prediction_sha256"]) == 64
    assert list(report["model_c2h6_selectivity_fraction_by_feed"].values()) == pytest.approx(
        [0.5331828914, 0.7396627624, 0.8358488448], rel=1e-5
    )
    assert list(report["model_h2_selectivity_fraction_by_feed"].values()) == pytest.approx(
        [0.3593332185, 0.3010887320, 0.2736481529], rel=1e-5
    )
    assert list(report["model_ch4_conversion_fraction_by_feed"].values()) == pytest.approx(
        [0.1713888487, 0.1723152070, 0.1752425993], rel=1e-5
    )


def test_feed_conversion_failure_persists_at_paper_average_field_estimate() -> None:
    sweep = audit_dry_feed_sweep(reduced_electric_field_Td=189.0)
    assert sweep["reduced_electric_field_Td_scenario"] == 189.0
    assert sweep["measured_product_assay_loaded"] is False
    report = evaluate_frozen_dry_feed_sweep(sweep)
    assert report["ch4_conversion_directional_trend_matches"] is False
    assert report["ch4_conversion_low_feed_coarse_plot_bound_matches"] is False
    assert report["ch4_conversion_high_feed_coarse_plot_bound_matches"] is True
    assert report["directional_trend_matches"] is True
    assert report["h2_directional_trend_matches"] is True
    assert [row["ch4_conversion_fraction"] for row in sweep["rows"]] == pytest.approx(
        [0.1666425182, 0.1676824568, 0.1709923434], rel=1e-5
    )
    assert all(row["material_ledger_closed"] for row in sweep["rows"])
    assert all(row["collision_energy_upper_bound_closed"] for row in sweep["rows"])
    assert all(row["n2_product_state_source_excluded"] for row in sweep["rows"])


def test_optical_field_estimate_reverses_slope_but_not_low_feed_magnitude() -> None:
    # 435 Td is a paper optical inference, not a verified local streamer field.
    sweep = audit_dry_feed_sweep(reduced_electric_field_Td=435.0)
    assert sweep["measured_product_assay_loaded"] is False
    conversions = [row["ch4_conversion_fraction"] for row in sweep["rows"]]
    assert conversions == pytest.approx(
        [0.1768312725, 0.1745949792, 0.1736584410], rel=1e-5
    )
    assert conversions == sorted(conversions, reverse=True)
    assert [row["ch4_electron_fragment_source_fraction_of_feed"] for row in sweep["rows"]] == pytest.approx(
        [0.2416441692, 0.2401731267, 0.2279019599], rel=1e-5
    )
    assert [row["ch4_net_neutral_chemistry_consumption_fraction_of_feed"] for row in sweep["rows"]] == pytest.approx(
        [-0.0648128967, -0.0655781475, -0.0542435189], rel=1e-5
    )
    report = evaluate_frozen_dry_feed_sweep(sweep)
    assert report["ch4_conversion_directional_trend_matches"] is True
    assert report["ch4_conversion_low_feed_coarse_plot_bound_matches"] is False
    assert report["ch4_conversion_coarse_plot_bounds_both_match"] is False
    assert report["like_for_like_quantitative_comparison_ready"] is False
    assert all(row["material_ledger_closed"] for row in sweep["rows"])
    assert all(row["collision_energy_upper_bound_closed"] for row in sweep["rows"])
    assert all(row["n2_product_state_source_excluded"] for row in sweep["rows"])


def test_optional_methane_pathway_trace_closes_neutral_consumption() -> None:
    sweep = audit_dry_feed_sweep(trace_methane_neutral_pathways=True)
    assert sweep["trace_methane_neutral_pathways"] is True
    for row in sweep["rows"]:
        diagnostic = row["methane_neutral_pathways"]
        assert diagnostic["not_a_validated_mechanism"] is True
        assert abs(diagnostic["quadrature_closure_residual_mol"]) < 1e-12
        assert diagnostic["integrated_chemical_net_mol"] == pytest.approx(
            -row["ch4_net_neutral_chemistry_consumption_fraction_of_feed"]
            * row["inlet_CH4_mol"],
            abs=1e-12,
        )
        assert diagnostic["reaction_contributions"][0]["equation"] == (
            "CH3 + H (+M) <=> CH4 (+M)"
        )
    assert [row["ch4_conversion_fraction"] for row in sweep["rows"]] == pytest.approx(
        [0.1713888487, 0.1723152070, 0.1752425993], rel=1e-5
    )
