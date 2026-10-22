"""Figure 1D conversion bands are external, not discharge-model inputs."""

import pytest

from evaluation.reference_only.dry_conversion_plot_bands import (
    evaluate_frozen_dry_residence_plot_bands,
)
from scripts.audit_dry_residence_sweep import (
    NOMINAL_RESIDENCE_TIMES_S,
    audit_dry_residence_sweep,
)


@pytest.fixture(scope="module")
def figure_1d_sweep() -> dict:
    return audit_dry_residence_sweep()


def test_four_residence_predictions_have_closed_ledgers(figure_1d_sweep: dict) -> None:
    assert [row["nominal_residence_s"] for row in figure_1d_sweep["rows"]] == list(
        NOMINAL_RESIDENCE_TIMES_S
    )
    assert figure_1d_sweep["measured_product_assay_loaded"] is False
    assert all(row["material_ledger_closed"] for row in figure_1d_sweep["rows"])
    assert all(row["collision_energy_upper_bound_closed"] for row in figure_1d_sweep["rows"])
    assert all(row["n2_product_state_source_excluded"] for row in figure_1d_sweep["rows"])
    assert [row["ch4_conversion_fraction"] for row in figure_1d_sweep["rows"]] == pytest.approx(
        [0.51856165, 0.30657148, 0.21637286, 0.16521672], rel=2e-5,
    )


def test_frozen_figure_1d_visual_bands_all_fail_without_fitting(figure_1d_sweep: dict) -> None:
    report = evaluate_frozen_dry_residence_plot_bands(figure_1d_sweep)
    assert report["visual_band_matches"] == 0
    assert report["visual_band_point_count"] == 4
    assert all(not row["model_inside_visual_band"] for row in report["comparison_rows"])
    assert all(row["minimum_absolute_fraction_gap_to_visual_band"] > 0.2 for row in report["comparison_rows"])
    assert report["visual_bands_are_experimental_confidence_intervals"] is False
    assert report["plot_read_intervals_are_exact_numeric_assays"] is False
    assert report["like_for_like_quantitative_comparison_ready"] is False
    assert report["reference_read_after_prediction_frozen"] is True
    assert report["reference_used_by_solver"] is False
    assert report["experimental_hcn_values_loaded"] is False
    assert len(report["frozen_prediction_sha256"]) == 64
