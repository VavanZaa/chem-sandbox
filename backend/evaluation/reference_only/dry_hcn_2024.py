"""Held-out qualitative dry-flow HCN comparison, never a solver dependency."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


REFERENCE_PATH = (
    Path(__file__).resolve().parents[2]
    / "data" / "benchmarks" / "kamarinopoulou_2024_dry_flow_qualitative_v1.json"
)


def evaluate_frozen_dry_flow_prediction(prediction: dict[str, Any]) -> dict[str, Any]:
    """Read the reference only after the caller has produced its prediction."""
    frozen_json = json.dumps(prediction, sort_keys=True, separators=(",", ":"), allow_nan=False)
    prediction_sha256 = hashlib.sha256(frozen_json.encode("utf-8")).hexdigest()
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    conditions = prediction["conditions"]
    outlet = prediction["outlet_packet_diagnostics"]
    yields = outlet["product_carbon_yield_fraction_by_species"]
    expected = reference["nominal_experimental_condition"]
    nominal_feed_power_residence_matches = (
        conditions["initial_CH4_mole_fraction"] == expected["feed_CH4_mole_fraction"]
        and conditions["nominal_residence_s"] == expected["nominal_residence_s"]
        and conditions["average_discharge_power_W"] == expected["average_discharge_power_W"]
    )
    temperature_proxy_value_matches = (
        conditions["uniform_gas_temperature_K_scenario"]
        == reference["thermometry_provenance"]["reported_steady_downstream_probe_temperature_C"] + 273.15
    )
    model_hcn_exceeds_ethane = yields["HCN"] > yields["C2H6"]
    unmatched_requirements = []
    if not nominal_feed_power_residence_matches:
        unmatched_requirements.append("nominal_feed_power_or_residence")
    if not temperature_proxy_value_matches:
        unmatched_requirements.append("downstream_probe_temperature_proxy_value")
    unmatched_requirements.append("reactor_temperature_profile_and_thermometry_run_match")
    if not prediction["flow_temperature_basis_diagnostic"]["comparison_basis_reconciled"]:
        unmatched_requirements.append("flow_reference_and_temperature_profile")
    if not prediction["scientific_firewall"]["flow_reactor_spatial_or_microdischarge_model_solved"]:
        unmatched_requirements.append("microdischarge_exposure_and_reactor_heterogeneity")
    if "N2_dissociation_outcome_000" in prediction["heavy_chemistry"]["excluded_fragment_channels_by_id"]:
        unmatched_requirements.append("n2_fragment_product_states")
    if not reference["held_out_observation"]["numeric_yields_digitized"]:
        unmatched_requirements.append("numeric_outlet_observations_and_uncertainties")
    like_for_like = not unmatched_requirements
    return {
        "benchmark_id": reference["benchmark_id"],
        "frozen_prediction_sha256": prediction_sha256,
        "source_doi": reference["source"]["doi"],
        "nominal_feed_power_residence_matches": nominal_feed_power_residence_matches,
        "uniform_temperature_equals_downstream_probe_value": temperature_proxy_value_matches,
        "thermometry_run_power_W": reference["thermometry_provenance"]["measurement_run_power_W"],
        "nominal_benchmark_power_W": expected["average_discharge_power_W"],
        "reactor_temperature_profile_measured": reference["thermometry_provenance"][
            "measured_inside_reactor_spatial_temperature_profile"
        ],
        "like_for_like_quantitative_comparison_ready": like_for_like,
        "unmatched_requirements": unmatched_requirements,
        "held_out_observation": reference["held_out_observation"]["comparison"],
        "model_hcn_carbon_yield_fraction": yields["HCN"],
        "model_c2h6_carbon_yield_fraction": yields["C2H6"],
        "model_hcn_exceeds_ethane": model_hcn_exceeds_ethane,
        "directional_ordering_matches": model_hcn_exceeds_ethane,
        "status": (
            "directional_match_under_unmatched_scenario"
            if model_hcn_exceeds_ethane else "directional_mismatch_under_unmatched_scenario"
        ) if not like_for_like else (
            "directional_match" if model_hcn_exceeds_ethane else "directional_mismatch"
        ),
        "reference_read_after_prediction_frozen": True,
        "reference_used_by_solver": False,
        "numeric_experimental_yields_digitized": reference["held_out_observation"]["numeric_yields_digitized"],
    }


def evaluate_frozen_dry_feed_sweep(sweep: dict[str, Any]) -> dict[str, Any]:
    """Compare a fully computed sweep with held-out conversion and selectivity trends."""
    frozen_json = json.dumps(sweep, sort_keys=True, separators=(",", ":"), allow_nan=False)
    prediction_sha256 = hashlib.sha256(frozen_json.encode("utf-8")).hexdigest()
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    observed = reference["held_out_feed_trend"]
    observed_conversion = reference["held_out_methane_conversion_feed_trend"]
    observed_conversion_bounds = reference["held_out_methane_conversion_plot_bounds"]
    observed_hydrogen = reference["held_out_hydrogen_feed_trend"]
    rows = sweep["rows"]
    fractions = [row["feed_CH4_mole_fraction"] for row in rows]
    selectivities = [row["c2h6_carbon_selectivity_fraction"] for row in rows]
    conversions = [row["ch4_conversion_fraction"] for row in rows]
    hydrogen_selectivities = [row["h2_hydrogen_selectivity_fraction"] for row in rows]
    conditions_match = fractions == observed["feed_CH4_mole_fractions"]
    directional_match = bool(
        conditions_match
        and len(selectivities) >= 2
        and all(a < b for a, b in zip(selectivities, selectivities[1:]))
    )
    hydrogen_conditions_match = fractions == observed_hydrogen["feed_CH4_mole_fractions"]
    hydrogen_directional_match = bool(
        hydrogen_conditions_match
        and len(hydrogen_selectivities) >= 2
        and all(a > b for a, b in zip(hydrogen_selectivities, hydrogen_selectivities[1:]))
    )
    conversion_conditions_match = fractions == observed_conversion[
        "feed_CH4_mole_fractions"
    ]
    conversion_directional_match = bool(
        conversion_conditions_match
        and len(conversions) >= 2
        and all(a > b for a, b in zip(conversions, conversions[1:]))
    )
    conversion_by_feed = dict(zip(fractions, conversions))
    low_feed = observed_conversion_bounds["feed_CH4_mole_fraction_low"]
    high_feed = observed_conversion_bounds["feed_CH4_mole_fraction_high"]
    low_feed_coarse_bound_matches = bool(
        low_feed in conversion_by_feed
        and conversion_by_feed[low_feed]
        > observed_conversion_bounds["low_feed_conversion_fraction_strictly_above"]
    )
    high_feed_coarse_bound_matches = bool(
        high_feed in conversion_by_feed
        and conversion_by_feed[high_feed]
        < observed_conversion_bounds["high_feed_conversion_fraction_strictly_below"]
    )
    return {
        "benchmark_id": reference["benchmark_id"],
        "frozen_prediction_sha256": prediction_sha256,
        "source_doi": reference["source"]["doi"],
        "held_out_observation": observed["comparison"],
        "feed_points_match_reference": conditions_match,
        "model_c2h6_selectivity_fraction_by_feed": dict(zip(fractions, selectivities)),
        "directional_trend_matches": directional_match,
        "held_out_h2_observation": observed_hydrogen["comparison"],
        "model_h2_selectivity_fraction_by_feed": dict(zip(fractions, hydrogen_selectivities)),
        "h2_directional_trend_matches": hydrogen_directional_match,
        "held_out_ch4_conversion_observation": observed_conversion["comparison"],
        "model_ch4_conversion_fraction_by_feed": conversion_by_feed,
        "ch4_conversion_directional_trend_matches": conversion_directional_match,
        "ch4_conversion_low_feed_coarse_plot_bound_matches": low_feed_coarse_bound_matches,
        "ch4_conversion_high_feed_coarse_plot_bound_matches": high_feed_coarse_bound_matches,
        "ch4_conversion_coarse_plot_bounds_both_match": (
            low_feed_coarse_bound_matches and high_feed_coarse_bound_matches
        ),
        "ch4_conversion_plot_bounds_are_statistical_uncertainty": False,
        "selectivity_basis_matches_paper_exactly": sweep[
            "selectivity_basis_matches_paper_exactly"
        ],
        "like_for_like_quantitative_comparison_ready": False,
        "status": (
            "all_directional_trends_match_under_unmatched_scenario"
            if (
                directional_match and hydrogen_directional_match
                and conversion_directional_match
                and low_feed_coarse_bound_matches and high_feed_coarse_bound_matches
            )
            else "one_or_more_directional_trends_mismatch_under_unmatched_scenario"
        ),
        "reference_read_after_prediction_frozen": True,
        "reference_used_by_solver": False,
        "numeric_experimental_selectivities_digitized": observed[
            "numeric_selectivities_digitized"
        ],
        "numeric_experimental_conversions_digitized": observed_conversion[
            "numeric_conversions_digitized"
        ],
    }


def evaluate_frozen_dry_residence_sweep(sweep: dict[str, Any]) -> dict[str, Any]:
    """Hold Figure 1D conversion observations outside the predictive solver."""
    frozen_json = json.dumps(sweep, sort_keys=True, separators=(",", ":"), allow_nan=False)
    prediction_sha256 = hashlib.sha256(frozen_json.encode("utf-8")).hexdigest()
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    observed = reference["held_out_residence_conversion"]
    rows = sweep["rows"]
    times = [row["nominal_residence_s"] for row in rows]
    conversions = [row["ch4_conversion_fraction"] for row in rows]
    conditions_match = times == observed["nominal_residence_times_s"]
    direction_matches = bool(
        conditions_match and len(conversions) >= 2
        and all(a > b for a, b in zip(conversions, conversions[1:]))
    )
    long_bound_matches = bool(
        conditions_match
        and conversions[0] > observed["long_residence_mean_conversion_strictly_above"]
    )
    short_bound_matches = bool(
        conditions_match
        and conversions[-1] > observed["short_residence_mean_conversion_strictly_above"]
    )
    return {
        "benchmark_id": reference["benchmark_id"],
        "frozen_prediction_sha256": prediction_sha256,
        "source_doi": reference["source"]["doi"],
        "held_out_observation": observed["comparison"],
        "nominal_residence_points_match_reference": conditions_match,
        "model_ch4_conversion_fraction_by_nominal_residence_s": dict(zip(times, conversions)),
        "conversion_directional_trend_matches": direction_matches,
        "long_residence_coarse_plot_bound_matches": long_bound_matches,
        "short_residence_coarse_plot_bound_matches": short_bound_matches,
        "coarse_plot_bounds_are_statistical_uncertainty": False,
        "like_for_like_quantitative_comparison_ready": False,
        "status": (
            "direction_and_coarse_bounds_match_under_unmatched_scenario"
            if direction_matches and long_bound_matches and short_bound_matches
            else "direction_or_coarse_bounds_mismatch_under_unmatched_scenario"
        ),
        "reference_read_after_prediction_frozen": True,
        "reference_used_by_solver": False,
        "numeric_experimental_conversions_digitized": observed[
            "numeric_conversions_digitized"
        ],
    }
