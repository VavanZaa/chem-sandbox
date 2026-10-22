"""Read-only intake audit for a locally obtained electron-collision set.

This deliberately does not install the candidate or change the active solver.
The owner must obtain and retain the full contributor file under its own terms;
the report contains metadata and channel labels, not the cross-section tables.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from bolos import parser

from app.prebiotic.boltzmann import _coverage
from app.prebiotic.collision_intake import (
    normalize_lxcat_bolsig_text,
    normalize_lxcat_processes,
)


def audit_candidate_collision_set(
    path: Path, target: str, *, normalize_lxcat_export: bool = False
) -> dict:
    path = path.resolve(strict=True)
    if not path.is_file():
        raise ValueError("The candidate collision path is not a file")
    raw_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    normalization_log: list[dict] = []
    source_text = path.read_text(encoding="utf-8")
    if normalize_lxcat_export:
        source_text, normalization_log = normalize_lxcat_bolsig_text(source_text, target)
    processes = parser.parse(io.StringIO(source_text))
    if normalize_lxcat_export:
        processes, process_changes = normalize_lxcat_processes(processes, target)
        normalization_log.extend(process_changes)
    selected = [process for process in processes if process.get("target") == target]
    if not selected:
        raise ValueError(f"The candidate contains no {target} collision processes")

    issues = []
    for index, process in enumerate(selected):
        data = np.asarray(process.get("data", []), dtype=float)
        if data.ndim != 2 or data.shape[1] != 2 or data.shape[0] < 2:
            issues.append(f"{target} process {index} has no usable two-column energy grid")
            continue
        energy, sigma = data.T
        if not np.all(np.isfinite(data)):
            issues.append(f"{target} process {index} contains a nonfinite value")
        if np.any(energy < 0) or np.any(np.diff(energy) <= 0):
            issues.append(f"{target} process {index} has an invalid energy grid")
        if np.any(sigma < 0):
            issues.append(f"{target} process {index} has a negative cross section")

    coverage = _coverage(selected, (target,))["targets"][target]
    momentum_kinds = Counter(
        process["kind"]
        for process in selected
        if process["kind"] in {"ELASTIC", "EFFECTIVE", "MOMENTUM"}
    )
    effective_count = momentum_kinds["EFFECTIVE"] + momentum_kinds["MOMENTUM"]
    if momentum_kinds["ELASTIC"] and effective_count:
        issues.append(
            f"{target} mixes ELASTIC with EFFECTIVE/MOMENTUM; BOLOS cannot load both"
        )
    if momentum_kinds["ELASTIC"] > 1 or effective_count > 1:
        issues.append(f"{target} has multiple momentum-transfer curves for one target")
    labels = list(dict.fromkeys(coverage["explicit_neutral_fragment_labels"]))
    unresolved = list(dict.fromkeys(coverage["unresolved_dissociation_labels"]))
    suggested_only = [
        str(process.get("product", ""))
        for process in selected
        if "suggested dissociation channel" in str(process.get("comment", "")).lower()
        and str(process.get("product", "")) not in labels
    ]
    if not coverage["complete"]:
        issues.append(f"{target} momentum, inelastic, or ionization input is missing")
    if not labels:
        issues.append("No explicit neutral-fragment labels were found")
    if unresolved:
        issues.append("One or more dissociative energy-loss labels have unresolved fragments")
    return {
        "status": "structural_intake_only_not_swarm_or_chemistry_validation",
        "candidate_sha256": raw_hash,
        "candidate_file_name": path.name,
        "target": target,
        "process_count": len(selected),
        "kind_counts": dict(sorted(Counter(process["kind"] for process in selected).items())),
        "lxcat_export_normalized_in_memory": normalize_lxcat_export,
        "normalization_log": normalization_log,
        "eedf_input_coverage": coverage["complete"],
        "explicit_neutral_fragment_labels": labels,
        "unresolved_dissociation_labels": unresolved,
        "suggested_dissociation_channel_labels_not_verified_as_products": suggested_only,
        "structural_issues": issues,
        "passes_preliminary_structure_check": not issues,
        "candidate_data_copied_or_installed": False,
        "active_solver_changed": False,
        "next_gates": [
            "Review provenance, contributor terms, channel identities and overlapping loss accounting manually.",
            f"Replace the entire {target} EEDF collision set; do not append overlapping loss curves to an existing set.",
            f"Compare the replacement against independent measured {target} swarm and transport data.",
            "Only then couple charge-balanced fragment source terms and retest wet apparatus predictions.",
        ],
    }


def audit_candidate_methane_collision_set(path: Path) -> dict:
    return audit_candidate_collision_set(path, "CH4")


def main() -> None:
    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("candidate", type=Path, help="Local BOLSIG+ complete-set text file")
    command.add_argument("--target", default="CH4", help="Collision target to audit (default: CH4)")
    command.add_argument(
        "--normalize-lxcat-export",
        action="store_true",
        help="Audit and normalize clear LXCat label annotations, products and duplicate rows in memory",
    )
    args = command.parse_args()
    print(
        json.dumps(
            audit_candidate_collision_set(
                args.candidate,
                args.target,
                normalize_lxcat_export=args.normalize_lxcat_export,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
