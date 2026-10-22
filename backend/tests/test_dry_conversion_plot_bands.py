"""Figure 1A visual conversion bands stay outside the dry predictive solver."""

import pytest

from evaluation.reference_only.dry_conversion_plot_bands import (
    evaluate_frozen_dry_conversion_plot_bands,
)
from scripts.audit_dry_feed_sweep import (
    FIGURE_1A_METHANE_FEED_FRACTIONS,
    audit_dry_feed_sweep,
)


@pytest.fixture(scope="module")
def figure_1a_sweep() -> dict:
    return audit_dry_feed_sweep(
        methane_feed_fractions=FIGURE_1A_METHANE_FEED_FRACTIONS,
    )


def test_all_seven_predictions_are_computed_before_plot_evaluation(
    figure_1a_sweep: dict,
) -> None:
    assert [row["feed_CH4_mole_fraction"] for row in figure_1a_sweep["rows"]] == list(
        FIGURE_1A_METHANE_FEED_FRACTIONS
    )
    assert figure_1a_sweep["measured_product_assay_loaded"] is False
    assert all(row["material_ledger_closed"] for row in figure_1a_sweep["rows"])
    assert all(row["collision_energy_upper_bound_closed"] for row in figure_1a_sweep["rows"])
    assert all(row["n2_product_state_source_excluded"] for row in figure_1a_sweep["rows"])
    assert [row["ch4_conversion_fraction"] for row in figure_1a_sweep["rows"]] == pytest.approx(
        [0.17138885, 0.17148156, 0.17166764, 0.17231521, 0.17340490, 0.17434991, 0.17524260],
        rel=2e-5,
    )


def test_frozen_figure_1a_comparison_is_a_mismatch_not_a_fit(figure_1a_sweep: dict) -> None:
    report = evaluate_frozen_dry_conversion_plot_bands(figure_1a_sweep)
    assert report["visual_band_matches"] == 1
    assert report["visual_band_point_count"] == 7
    assert [row["model_inside_visual_band"] for row in report["comparison_rows"]] == [
        False, False, False, False, False, False, True,
    ]
    assert report["comparison_rows"][0]["minimum_absolute_fraction_gap_to_visual_band"] > 0.57
    assert report["visual_bands_are_experimental_confidence_intervals"] is False
    assert report["plot_read_intervals_are_exact_numeric_assays"] is False
    assert report["like_for_like_quantitative_comparison_ready"] is False
    assert report["reference_read_after_prediction_frozen"] is True
    assert report["reference_used_by_solver"] is False
    assert report["experimental_hcn_values_loaded"] is False
    assert len(report["frozen_prediction_sha256"]) == 64


@pytest.mark.parametrize("fractions", [(), (0.01, 0.005), (0.0,), (1.0,)])
def test_invalid_feed_grids_rejected_before_solving(fractions: tuple[float, ...]) -> None:
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_dry_feed_sweep(methane_feed_fractions=fractions)
