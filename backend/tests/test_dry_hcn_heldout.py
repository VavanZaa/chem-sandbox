"""The published observation is evaluator-only, never a solver input."""
from __future__ import annotations

import json

from evaluation.reference_only.dry_hcn_2024 import (
    REFERENCE_PATH,
    evaluate_frozen_dry_flow_prediction,
)
from scripts.run_dry_flow_hcn_probe import run_dry_flow_probe


def test_heldout_observation_is_sourced_and_not_an_input() -> None:
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    assert reference["source"]["doi"] == "10.1126/sciadv.adl4246"
    assert reference["held_out_observation"]["numeric_yields_digitized"] is False
    prediction = run_dry_flow_probe(source_forcing="evolving_eedf", substeps=16)
    report = evaluate_frozen_dry_flow_prediction(prediction)
    assert report["nominal_feed_power_residence_matches"] is True
    assert report["uniform_temperature_equals_downstream_probe_value"] is True
    assert report["thermometry_run_power_W"] == 28.6
    assert report["nominal_benchmark_power_W"] == 28.1
    assert report["reactor_temperature_profile_measured"] is False
    assert "reactor_temperature_profile_and_thermometry_run_match" in report[
        "unmatched_requirements"
    ]
    assert report["like_for_like_quantitative_comparison_ready"] is False
    assert report["model_hcn_carbon_yield_fraction"] < report["model_c2h6_carbon_yield_fraction"]
    assert report["status"] == "directional_mismatch_under_unmatched_scenario"
    assert report["reference_used_by_solver"] is False
    assert len(report["frozen_prediction_sha256"]) == 64
    assert prediction["source"]["measured_hcn_used_as_solver_input"] is False
    assert prediction["outlet_packet_diagnostics"]["measured_product_values_used_as_inputs"] is False
