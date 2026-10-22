from __future__ import annotations

from typing import Any
from collections import Counter


def build_statistics(
    graph: dict[str, Any],
    *,
    total_runs: int,
    valid_states: int,
    failed_optimizations: int,
    rejected_candidate_states: int = 0,
    outcomes: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    co2_nodes = [
        node for node in graph["nodes"] if "CO2" in {piece.strip() for piece in node["formula"].split("+")}
    ]
    co2_runs = len({run for node in co2_nodes for run in node.get("run_indices", [])})
    co2_terminal_runs = len({run for node in co2_nodes for run in node.get("endpoint_run_indices", [])})
    outcomes = outcomes or []
    return {
        "search_runs": total_runs,
        "accounted_runs": len(outcomes),
        "run_status_counts": dict(Counter(item["status"] for item in outcomes)),
        "endpoint_reason_counts": dict(Counter(endpoint["reason"] for item in outcomes for endpoint in item["endpoints"])),
        "candidate_evaluations": sum(item["evaluations"] for item in outcomes),
        "initialization_evaluations": sum(item.get("initialization_evaluations", 0) for item in outcomes),
        "valid_optimized_states": valid_states,
        "failed_optimizations": failed_optimizations,
        "rejected_candidate_states": rejected_candidate_states,
        "unique_states": len(graph["nodes"]),
        "candidate_pathways": len(graph["edges"]),
        "co2_containing_runs": co2_runs,
        "co2_run_frequency": co2_runs / max(total_runs, 1),
        "co2_terminal_runs": co2_terminal_runs,
        "co2_terminal_frequency": co2_terminal_runs / max(total_runs, 1),
        "disclaimer": "Discovery and branch-stop frequencies can overlap across states. Neither is a yield or physical branching probability.",
    }
