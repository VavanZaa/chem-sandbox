"""Held-out Miller-Urey comparison loaded only after discovery is frozen."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from scipy.constants import R


DATA_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "benchmarks"
    / "miller_urey_1953_v1.json"
)
PROFILE_PATH = DATA_PATH.with_name("miller_urey_1953_profile.json")
WATER_MOLAR_VOLUME_ML = 18.01528
PA_PER_MMHG = 133.32236842105263


def historical_material_time_profile() -> dict[str, Any]:
    """Return the versioned profile; unknown electrical inputs stay hypothetical."""

    return json.loads(PROFILE_PATH.read_text(encoding="utf-8"))


def evaluate_historical_benchmark(
    discovery: dict[str, Any],
    discovery_sha256: str,
    initial_mol: dict[str, float],
    final_liquid_mol: dict[str, float],
    gas_volume_L: float,
    duration_s: float,
) -> dict[str, Any]:
    """Compare a frozen result with held-out observations without feedback."""

    benchmark = json.loads(DATA_PATH.read_text(encoding="utf-8"))
    conditions = benchmark["apparatus_conditions"]
    discovered_smiles = {
        candidate["smiles"]
        for record in discovery.get("formula_records", [])
        for candidate in record.get("candidates", [])
    }
    reference_temperature_K = 298.15
    gas_volume_m3 = gas_volume_L * 1e-3
    inferred_pressures = {
        species: (
            float(initial_mol.get(species, 0.0))
            * R
            * reference_temperature_K
            / gas_volume_m3
            / PA_PER_MMHG
        )
        if gas_volume_m3 > 0
        else 0.0
        for species in conditions["gas_partial_pressures_mmHg"]
    }
    water_mL = float(initial_mol.get("H2O", 0.0)) * WATER_MOLAR_VOLUME_ML

    def close(actual: float, expected: float, tolerance: float = 0.1) -> bool:
        return math.isclose(actual, expected, rel_tol=tolerance, abs_tol=0.0)

    condition_checks = {
        "water_amount": {
            "actual": water_mL,
            "target": conditions["water_mL"],
            "unit": "mL",
            "matched": close(water_mL, float(conditions["water_mL"])),
        },
        "duration": {
            "actual": duration_s,
            "target": conditions["duration_s"],
            "unit": "s",
            "matched": close(duration_s, float(conditions["duration_s"]), 0.01),
        },
    }
    for species, target in conditions["gas_partial_pressures_mmHg"].items():
        condition_checks[f"partial_pressure_{species}"] = {
            "actual": inferred_pressures[species],
            "target": float(target),
            "unit": "mmHg",
            "matched": close(inferred_pressures[species], float(target)),
        }

    observations = []
    for observation in benchmark["evaluation_only_observations"]:
        candidate_hits = sorted(set(observation["acceptable_smiles"]) & discovered_smiles)
        precursors = {
            species: float(final_liquid_mol.get(species, 0.0))
            for species in observation["precursor_requirements"]
        }
        precursor_supported = bool(precursors) and all(value > 0 for value in precursors.values())
        observations.append(
            {
                "name": observation["name"],
                "candidate_recovered": bool(candidate_hits),
                "candidate_smiles": candidate_hits,
                "captured_precursors_mol": precursors,
                "precursor_supported": precursor_supported,
                "pathway_status": (
                    "candidate_and_precursors_only"
                    if candidate_hits and precursor_supported
                    else "candidate_only"
                    if candidate_hits
                    else "not_recovered"
                ),
                "source": observation["source"],
            }
        )

    conditions_matched = all(item["matched"] for item in condition_checks.values())
    enumeration_coverage = sum(item["candidate_recovered"] for item in observations) / len(
        observations
    )
    remaining_blockers = [
        "The historical discharge current and deposited energy are not supplied by the benchmark sources used here.",
        "Enumeration coverage only asks whether an identity occurs anywhere in the broad candidate set; it is not a ranked prediction, pathway calculation, or analytical detection prediction.",
        "Alpha- and beta-alanine still lack a complete evidence-backed aqueous formation network.",
    ]
    if not conditions_matched:
        remaining_blockers.append(
            "The current scenario must match the one-week material and pressure profile before historical comparison is valid."
        )
    return {
        "status": (
            "historical_profile_partial_not_reproduced"
            if conditions_matched
            else "development_partial_not_reproduced"
        ),
        "benchmark_id": benchmark["benchmark_id"],
        "discovery_artifact_sha256": discovery_sha256,
        "historical_conditions_matched": conditions_matched,
        "condition_checks": condition_checks,
        "observation_comparison": observations,
        "enumeration_coverage": enumeration_coverage,
        "enumeration_coverage_is_prediction_score": False,
        "validated_reproduction": False,
        "remaining_blockers": remaining_blockers,
        "sources": {
            "apparatus": conditions["source"],
            "protocol": conditions["protocol_source"],
            "mechanism": benchmark["mechanism_source"],
        },
        "scientific_firewall": {
            "benchmark_loaded_after_discovery_frozen": True,
            "benchmark_observations_available_to_discovery": False,
            "benchmark_yields_used_for_fitting": False,
            "evaluation_feedback_to_simulation": False,
        },
    }
