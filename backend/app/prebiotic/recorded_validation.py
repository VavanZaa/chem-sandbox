"""Read-only display of recorded prebiotic development diagnostics.

This module is imported by the API route only. The apparatus solver and blind
candidate generator do not import these records or use their reference values.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


RECORDS = Path(__file__).resolve().parents[2] / "data" / "validation" / "prebiotic_latest"


def _load(name: str) -> dict[str, Any] | None:
    path = RECORDS / name
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else None


def recorded_prebiotic_validation() -> dict[str, Any]:
    runs = [
        item
        for item in (
            _load("drift_power_capped_168h_summary.json"),
            _load("continuous_overlap_168h_summary.json"),
            _load("time_averaged_source_168h_summary.json"),
            _load("cold_trap_holdback_168h_summary.json"),
            _load("flow_limited_capture_168h_summary.json"),
            _load("coupled_flow_transport_168h_summary.json"),
        )
        if item is not None
    ]
    convergence = [
        item
        for item in (
            _load("drift_power_capped_convergence_24h.json"),
            _load("timestep_convergence_24h.json"),
            _load("time_averaged_source_convergence_24h.json"),
            _load("cold_trap_holdback_convergence_24h.json"),
            _load("flow_limited_capture_convergence_24h.json"),
            _load("coupled_flow_transport_convergence_24h.json"),
        )
        if item is not None
    ]
    sensitivity = [
        item for item in (
            _load("coupled_flow_contact_sensitivity_24h.json"),
        ) if item is not None
    ]
    return {
        "status": "recorded_scenarios_not_validated_reproduction",
        "runs": runs,
        "convergence_studies": convergence,
        "sensitivity_studies": sensitivity,
        "scientific_firewall": {
            "reference_yields_used_for_discovery": False,
            "reference_yields_used_for_solver_fitting": False,
            "quantitative_yield_claim_allowed": False,
        },
    }
