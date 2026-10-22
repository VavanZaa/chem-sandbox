"""Held-out Figure 1B molecular-hydrogen selectivity visual comparison."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


REFERENCE_PATH = (
    Path(__file__).resolve().parents[2]
    / "data" / "benchmarks"
    / "kamarinopoulou_2024_figure1b_hydrogen_selectivity_plot_bands_v1.json"
)


def evaluate_frozen_dry_hydrogen_plot_bands(sweep: dict[str, Any]) -> dict[str, Any]:
    frozen_json = json.dumps(sweep, sort_keys=True, separators=(",", ":"), allow_nan=False)
    frozen_sha256 = hashlib.sha256(frozen_json.encode("utf-8")).hexdigest()
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    measured_bands = reference["hydrogen_selectivity_plot_bands"]
    predictions = sweep["rows"]
    if [row["feed_CH4_mole_fraction"] for row in predictions] != [
        row["feed_CH4_mole_fraction"] for row in measured_bands
    ]:
        raise ValueError("Prediction feed points do not match the frozen Figure 1B x positions")

    comparisons = []
    for predicted, measured in zip(predictions, measured_bands):
        value = float(predicted["h2_hydrogen_selectivity_fraction"])
        low = float(measured["lower_fraction"])
        high = float(measured["upper_fraction"])
        if not (0 <= low <= measured["visual_central_fraction"] <= high <= 1):
            raise ValueError("The frozen plot-read band is not ordered within [0, 1]")
        comparisons.append({
            "feed_CH4_mole_fraction": measured["feed_CH4_mole_fraction"],
            "model_H2_hydrogen_selectivity_fraction": value,
            "plot_read_central_fraction": measured["visual_central_fraction"],
            "plot_read_lower_fraction": low,
            "plot_read_upper_fraction": high,
            "model_inside_visual_band": low <= value <= high,
            "minimum_absolute_fraction_gap_to_visual_band": max(low - value, value - high, 0.0),
            "material_ledger_closed": predicted["material_ledger_closed"],
            "collision_energy_upper_bound_closed": predicted["collision_energy_upper_bound_closed"],
            "n2_product_state_source_excluded": predicted["n2_product_state_source_excluded"],
        })
    return {
        "benchmark_id": reference["benchmark_id"],
        "source_doi": reference["source"]["doi"],
        "frozen_prediction_sha256": frozen_sha256,
        "comparison_rows": comparisons,
        "visual_band_matches": sum(row["model_inside_visual_band"] for row in comparisons),
        "visual_band_point_count": len(comparisons),
        "visual_bands_are_experimental_confidence_intervals": reference["extraction"][
            "bands_are_experimental_confidence_intervals"
        ],
        "model_selectivity_basis_matches_paper_exactly": sweep["selectivity_basis_matches_paper_exactly"],
        "like_for_like_quantitative_comparison_ready": False,
        "unmatched_conditions": [
            "actual_in_reactor_temperature_profile",
            "actual_local_streamer_exposure_and_gas_exchange",
            "flow_reference_and_residence_convention",
            "wall_deposited_carbon_and_complete_stable_product_basis",
        ],
        "status": (
            "all_plot_read_bands_match_under_unmatched_scenario"
            if all(row["model_inside_visual_band"] for row in comparisons)
            else "one_or_more_plot_read_bands_fail_under_unmatched_scenario"
        ),
        "reference_read_after_prediction_frozen": True,
        "reference_used_by_solver": False,
        "experimental_hcn_values_loaded": False,
    }
