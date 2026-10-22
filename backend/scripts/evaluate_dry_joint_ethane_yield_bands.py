"""Derive held-out C2H6 feed-carbon yield bands from Figure 1A readings.

The assay-free feed predictions are frozen before either non-HCN plot card is
loaded. Multiplying conversion by C2H6 carbon selectivity is an algebraic
comparison, not a new measured yield or a fit to product observations.
"""
from __future__ import annotations

import json
import math
from typing import Any

from scripts.audit_dry_feed_sweep import (
    FIGURE_1A_METHANE_FEED_FRACTIONS,
    audit_dry_feed_sweep,
)


def run_heldout_joint_ethane_yield_check() -> dict[str, Any]:
    sweep = audit_dry_feed_sweep(
        methane_feed_fractions=FIGURE_1A_METHANE_FEED_FRACTIONS,
    )
    # The reference-only modules are imported only after the solver returns.
    from evaluation.reference_only.dry_conversion_plot_bands import (
        evaluate_frozen_dry_conversion_plot_bands,
    )
    from evaluation.reference_only.dry_ethane_plot_bands import (
        evaluate_frozen_dry_ethane_plot_bands,
    )

    conversion = evaluate_frozen_dry_conversion_plot_bands(sweep)
    ethane = evaluate_frozen_dry_ethane_plot_bands(sweep)
    if conversion["frozen_prediction_sha256"] != ethane["frozen_prediction_sha256"]:
        raise ValueError("The two reference cards did not evaluate the same frozen prediction")
    if conversion["source_doi"] != ethane["source_doi"]:
        raise ValueError("The two Figure 1A reference cards cite different sources")
    if conversion["visual_bands_are_experimental_confidence_intervals"] or ethane[
        "visual_bands_are_experimental_confidence_intervals"
    ]:
        raise ValueError("Plot-reading bands must not be recast as confidence intervals")

    rows = []
    for prediction, converted, selected in zip(
        sweep["rows"], conversion["comparison_rows"], ethane["comparison_rows"], strict=True
    ):
        feed = prediction["feed_CH4_mole_fraction"]
        if feed != converted["feed_CH4_mole_fraction"] or feed != selected[
            "feed_CH4_mole_fraction"
        ]:
            raise ValueError("Conversion and ethane plot rows are not aligned")
        model_yield = prediction["c2h6_carbon_yield_fraction"]
        model_product = (
            prediction["ch4_conversion_fraction"]
            * prediction["c2h6_carbon_selectivity_fraction"]
        )
        if not math.isclose(model_yield, model_product, rel_tol=1e-8, abs_tol=1e-12):
            raise ValueError("Model ethane yield does not equal conversion times selectivity")
        lower = (
            converted["plot_read_lower_fraction"]
            * selected["plot_read_lower_fraction"]
        )
        upper = (
            converted["plot_read_upper_fraction"]
            * selected["plot_read_upper_fraction"]
        )
        rows.append({
            "feed_CH4_mole_fraction": feed,
            "derived_plot_read_C2H6_feed_carbon_yield_lower_fraction": lower,
            "derived_plot_read_C2H6_feed_carbon_yield_upper_fraction": upper,
            "model_C2H6_feed_carbon_yield_fraction": model_yield,
            "model_inside_derived_visual_band": lower <= model_yield <= upper,
            "minimum_absolute_fraction_gap_to_derived_visual_band": max(
                lower - model_yield, model_yield - upper, 0.0
            ),
            "material_ledger_closed": prediction["material_ledger_closed"],
            "collision_energy_upper_bound_closed": prediction[
                "collision_energy_upper_bound_closed"
            ],
            "n2_product_state_source_excluded": prediction[
                "n2_product_state_source_excluded"
            ],
        })
    return {
        "status": "derived_non_hcn_plot_read_yield_screen_not_like_for_like_validation",
        "source_doi": conversion["source_doi"],
        "frozen_prediction_sha256": conversion["frozen_prediction_sha256"],
        "conversion_benchmark_id": conversion["benchmark_id"],
        "ethane_benchmark_id": ethane["benchmark_id"],
        "rows": rows,
        "derived_bands_are_measured_product_yields": False,
        "derived_bands_are_experimental_confidence_intervals": False,
        "experimental_hcn_values_loaded": False,
        "reference_used_by_solver": False,
        "like_for_like_quantitative_comparison_ready": False,
        "unmatched_conditions": sorted(set(
            conversion["unmatched_conditions"] + ethane["unmatched_conditions"]
        )),
    }


if __name__ == "__main__":
    print(json.dumps(run_heldout_joint_ethane_yield_check(), indent=2))
