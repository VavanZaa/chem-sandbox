"""Held-out hydrogen bars remain an evaluator-only non-HCN observable."""

import pytest

from evaluation.reference_only.dry_hydrogen_plot_bands import (
    evaluate_frozen_dry_hydrogen_plot_bands,
)
from scripts.audit_dry_feed_sweep import (
    FIGURE_1A_METHANE_FEED_FRACTIONS,
    audit_dry_feed_sweep,
)


@pytest.fixture(scope="module")
def figure_1b_sweep() -> dict:
    return audit_dry_feed_sweep(
        methane_feed_fractions=FIGURE_1A_METHANE_FEED_FRACTIONS,
    )


def test_hydrogen_curve_is_frozen_after_assay_free_prediction(figure_1b_sweep: dict) -> None:
    assert figure_1b_sweep["measured_product_assay_loaded"] is False
    report = evaluate_frozen_dry_hydrogen_plot_bands(figure_1b_sweep)
    assert report["visual_band_matches"] == 0
    assert report["visual_band_point_count"] == 7
    assert all(not row["model_inside_visual_band"] for row in report["comparison_rows"])
    assert all(row["minimum_absolute_fraction_gap_to_visual_band"] > 0.1 for row in report["comparison_rows"])
    assert all(row["material_ledger_closed"] for row in report["comparison_rows"])
    assert all(row["collision_energy_upper_bound_closed"] for row in report["comparison_rows"])
    assert all(row["n2_product_state_source_excluded"] for row in report["comparison_rows"])
    assert report["visual_bands_are_experimental_confidence_intervals"] is False
    assert report["model_selectivity_basis_matches_paper_exactly"] is False
    assert report["like_for_like_quantitative_comparison_ready"] is False
    assert report["reference_read_after_prediction_frozen"] is True
    assert report["reference_used_by_solver"] is False
    assert report["experimental_hcn_values_loaded"] is False
    assert len(report["frozen_prediction_sha256"]) == 64
