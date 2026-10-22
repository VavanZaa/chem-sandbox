"""Held-out conversion and selectivity imply a coarse ethane-yield screen."""

import pytest

from scripts.evaluate_dry_joint_ethane_yield_bands import (
    run_heldout_joint_ethane_yield_check,
)


def test_joint_ethane_yield_screen_fails_without_hcn_or_solver_leakage() -> None:
    report = run_heldout_joint_ethane_yield_check()
    assert len(report["frozen_prediction_sha256"]) == 64
    assert report["experimental_hcn_values_loaded"] is False
    assert report["reference_used_by_solver"] is False
    assert report["derived_bands_are_measured_product_yields"] is False
    assert report["derived_bands_are_experimental_confidence_intervals"] is False
    assert report["like_for_like_quantitative_comparison_ready"] is False
    assert len(report["rows"]) == 7
    assert all(not row["model_inside_derived_visual_band"] for row in report["rows"])
    assert all(row["material_ledger_closed"] for row in report["rows"])
    assert all(row["collision_energy_upper_bound_closed"] for row in report["rows"])
    assert all(row["n2_product_state_source_excluded"] for row in report["rows"])
    one_percent = next(
        row for row in report["rows"] if row["feed_CH4_mole_fraction"] == 0.01
    )
    assert one_percent["derived_plot_read_C2H6_feed_carbon_yield_lower_fraction"] == pytest.approx(0.0117)
    assert one_percent["derived_plot_read_C2H6_feed_carbon_yield_upper_fraction"] == pytest.approx(0.0637)
    assert one_percent["model_C2H6_feed_carbon_yield_fraction"] == pytest.approx(0.1274551420)
    assert one_percent["model_C2H6_feed_carbon_yield_fraction"] > 2 * one_percent[
        "derived_plot_read_C2H6_feed_carbon_yield_upper_fraction"
    ]
