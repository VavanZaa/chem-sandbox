"""Cross-engine discovery agreement without converting model votes to probability."""
from __future__ import annotations

from collections import defaultdict


def build_discovery_comparison(engine_results, excluded, failures=None):
    failures = failures or []
    states = defaultdict(lambda: {"engines": [], "labels": set(), "formulas": set()})
    paths = defaultdict(lambda: {"engines": [], "labels": set(), "delta_energy_ev_by_engine": {}})
    summaries = []
    for engine_id, result, result_job_id in engine_results:
        summaries.append({"engine": engine_id, "status": result["status"], "model": result["model"],
                          "result_job_id": result_job_id, "states": len(result["graph"]["nodes"]),
                          "paths": len(result["graph"]["edges"]), "runtime_seconds": result["runtime_seconds"],
                          "failed_optimizations": result["statistics"]["failed_optimizations"]})
        node_by_id = {node["id"]: node for node in result["graph"]["nodes"]}
        for node in result["graph"]["nodes"]:
            record = states[node["comparison_key"]]
            record["engines"].append(engine_id)
            record["labels"].add(node["label"])
            record["formulas"].add(node["formula"])
        for edge in result["graph"]["edges"]:
            record = paths[edge["comparison_key"]]
            record["engines"].append(engine_id)
            source, target = node_by_id[edge["source"]], node_by_id[edge["target"]]
            record["labels"].add(f"{source['label']} → {target['label']}")
            record["delta_energy_ev_by_engine"][engine_id] = edge["delta_energy_ev"]
    successful = len(engine_results)
    state_rows = [{"comparison_key": key, "labels": sorted(value["labels"]),
                   "formulas": sorted(value["formulas"]), "engines": value["engines"],
                   "model_count": len(value["engines"]),
                   "agreement_fraction": len(value["engines"]) / successful}
                  for key, value in states.items()]
    path_rows = [{"comparison_key": key, "labels": sorted(value["labels"]),
                  "engines": value["engines"], "model_count": len(value["engines"]),
                  "agreement_fraction": len(value["engines"]) / successful,
                  "delta_energy_ev_by_engine": value["delta_energy_ev_by_engine"]}
                 for key, value in paths.items()]
    state_rows.sort(key=lambda row: (-row["model_count"], row["labels"]))
    path_rows.sort(key=lambda row: (-row["model_count"], row["labels"]))
    return {"mode": "independent_searches", "successful_engine_count": successful,
            "engine_runs": summaries, "failed_engines": failures, "excluded_engines": excluded,
            "states": state_rows, "paths": path_rows, "calibrated_confidence": None,
            "warnings": [
                "Agreement fraction is the fraction of successful compatible model searches that discovered a connectivity; it is not a reaction probability or confidence level.",
                "Each engine performs an independent bounded search with the same starting inputs, seed and budgets. Different relaxed geometries can lead to different proposed branches.",
                "A missing state or path means it was not found within this bounded search, not that the chemistry is impossible.",
                "Models can share reference methods, chemical data and systematic bias; agreement is not independent experimental confirmation.",
                "Only within-engine energy differences for matching atom inventories are shown. Absolute energies are never compared across models or AIMNet families.",
            ]}
