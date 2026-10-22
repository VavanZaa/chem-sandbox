r"""Read-only integrity check for a saved v0.2 run. No product target is asserted.

Usage from backend: ..\.venv\Scripts\python.exe scripts\verify_saved_result.py <job-id>
"""
import json
import math
from pathlib import Path
import re
import sys


def verify(job_id: str) -> dict:
    if not re.fullmatch(r"[A-Za-z0-9_-]{1,64}", job_id):
        raise ValueError("Invalid job id")
    backend = Path(__file__).resolve().parents[1]
    result = json.loads((backend / "data" / "runs" / job_id / "result.json").read_text(encoding="utf-8"))
    runs = result["search_settings"]["runs"]
    assert len(result["run_outcomes"]) == runs
    assert {item["run"] for item in result["run_outcomes"]} == set(range(runs))
    assert result["statistics"]["accounted_runs"] == runs
    assert sum(result["statistics"]["run_status_counts"].values()) == runs
    node_ids = {node["id"] for node in result["graph"]["nodes"]}
    for run in result["run_outcomes"]:
        assert run["endpoints"]
        assert run["evaluations"] <= result["search_settings"]["max_evaluations_per_run"]
        assert all(item["state_id"] is None or item["state_id"] in node_ids for item in run["endpoints"])
        assert sum(count for reason, count in run["attempt_counts"].items() if reason != "not_evaluated_budget") == run["evaluations"]
    for node in result["graph"]["nodes"]:
        assert math.isfinite(node["energy_ev"])
        assert len(node["run_indices"]) == node["runs_seen"] <= runs
        assert len(node["endpoint_run_indices"]) == node["terminal_runs"] <= runs
        assert (backend / "data" / node["structure_reference"]).is_file()
    for error in result["errors"]:
        if error.get("structure_reference"):
            assert (backend / "data" / error["structure_reference"]).is_file()
    return {"job_id": job_id, "verified": True, "statistics": result["statistics"], "runtime_seconds": result["runtime_seconds"]}


if __name__ == "__main__":
    print(json.dumps(verify(sys.argv[1]), indent=2))
