"""Run the balanced reaction suite on one or more installed GPU/CPU engines."""
from __future__ import annotations

import argparse
import copy
import csv
import gc
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from ase.io import write

from evaluation.reference_only.reaction_suite import (
    ALL_ENGINES,
    REACTIONS,
    SPECIES,
    build_geometry,
    inferred_connectivity,
)
from app.models.registry import make_engine
from evaluation.protocol import freeze_prediction_artifact, unlock_evaluation

EV_TO_KJ_MOL = 96.4853321233


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cuda", choices=("cuda", "cpu", "auto"))
    parser.add_argument("--engines", nargs="+", default=list(ALL_ENGINES))
    parser.add_argument("--cases", nargs="+", default=[case.id for case in REACTIONS])
    parser.add_argument("--fmax", type=float, default=0.08)
    parser.add_argument("--steps", type=int, default=120)
    parser.add_argument("--seed", type=int, default=20260922)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "data" / "benchmarks" / "reaction_suite_v01")
    parser.add_argument("--summarize-existing", type=Path, help="Rebuild CSV/Markdown and aggregate metrics from an existing JSON result")
    return parser.parse_args()


def release_engine(engine):
    del engine
    gc.collect()
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    except Exception:
        pass


def add_comparison_summary(payload: dict) -> None:
    valid = [row for row in payload["results"] if row.get("predicted_kj_mol") is not None]
    engine_summary = []
    for engine in payload.get("requested_engines", ALL_ENGINES):
        rows = [row for row in payload["results"] if row["engine"] == engine]
        referenced = [row for row in rows if row.get("raw_deviation_kj_mol") is not None]
        absolute = sorted(abs(row["raw_deviation_kj_mol"]) for row in referenced)
        engine_summary.append({
            "engine": engine,
            "supported_cases": sum(row.get("status") != "unsupported" for row in rows),
            "valid_cases": sum(row.get("predicted_kj_mol") is not None for row in rows),
            "reference_cases": len(referenced),
            "sign_agreement_cases": sum(row.get("sign_agreement") is True for row in referenced),
            "diagnostic_mae_kj_mol": None if not absolute else sum(absolute) / len(absolute),
            "diagnostic_median_abs_deviation_kj_mol": None if not absolute else absolute[len(absolute) // 2],
        })
    case_summary = []
    for case in payload["cases"]:
        rows = [row for row in valid if row["reaction_id"] == case["id"]]
        predictions = [row["predicted_kj_mol"] for row in rows]
        referenced = [row for row in rows if row.get("raw_deviation_kj_mol") is not None]
        closest = min(referenced, key=lambda row: abs(row["raw_deviation_kj_mol"])) if referenced else None
        case_summary.append({
            "reaction_id": case["id"],
            "reaction_name": case["name"],
            "models_with_result": len(rows),
            "prediction_min_kj_mol": min(predictions) if predictions else None,
            "prediction_max_kj_mol": max(predictions) if predictions else None,
            "cross_model_spread_kj_mol": max(predictions) - min(predictions) if predictions else None,
            "sign_agreement_models": sum(row.get("sign_agreement") is True for row in referenced),
            "referenced_models": len(referenced),
            "closest_model_diagnostic": None if closest is None else closest["engine"],
            "closest_raw_deviation_kj_mol": None if closest is None else closest["raw_deviation_kj_mol"],
        })
    payload["comparison_summary"] = {"by_engine": engine_summary, "by_case": case_summary}


def markdown_report(payload: dict) -> str:
    lines = [
        "# Reaction benchmark v0.1",
        "",
        f"Generated: {payload['generated_at']}",
        "",
        "Model values are optimized isolated-molecule electronic potential-energy differences. NIST references are gas-phase ΔH°298 values and include thermal/ZPE effects, so raw deviations are diagnostic only. Catalysts are not explicit unless stated; endpoint thermochemistry cannot validate catalytic barriers.",
        "",
        "| Case | Engine | ΔE model (kJ/mol) | Reference | Raw deviation | Converged | Connectivity | Time (s) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["results"]:
        if row.get("status") == "unsupported":
            lines.append(
                f"| {row['reaction_name']} | {row['engine']} | unsupported | "
                f"{'—' if row.get('reference_kj_mol') is None else f'{row['reference_kj_mol']:.1f}'} | — | — | — | — |"
            )
            continue
        predicted = "—" if row.get("predicted_kj_mol") is None else f"{row['predicted_kj_mol']:.1f}"
        reference = "—" if row.get("reference_kj_mol") is None else f"{row['reference_kj_mol']:.1f}"
        deviation = "—" if row.get("raw_deviation_kj_mol") is None else f"{row['raw_deviation_kj_mol']:+.1f}"
        lines.append(
            f"| {row['reaction_name']} | {row['engine']} | {predicted} | {reference} | {deviation} | "
            f"{'yes' if row.get('all_converged') else 'no'} | {'preserved' if row.get('connectivity_preserved') else 'changed'} | {row.get('wall_seconds', 0):.1f} |"
        )
    lines.extend(["", "## Engine summary", "",
                  "The deviation columns below are diagnostics against ΔH°298, not a formal electronic-energy benchmark.", "",
                  "| Engine | Valid / supported | Sign agreement | Diagnostic MAE (kJ/mol) | Median absolute deviation |",
                  "|---|---:|---:|---:|---:|"])
    for row in payload["comparison_summary"]["by_engine"]:
        mae = "—" if row["diagnostic_mae_kj_mol"] is None else f"{row['diagnostic_mae_kj_mol']:.1f}"
        median = "—" if row["diagnostic_median_abs_deviation_kj_mol"] is None else f"{row['diagnostic_median_abs_deviation_kj_mol']:.1f}"
        lines.append(f"| {row['engine']} | {row['valid_cases']} / {row['supported_cases']} | "
                     f"{row['sign_agreement_cases']} / {row['reference_cases']} | {mae} | {median} |")
    lines.extend(["", "## Cross-model spread", "",
                  "| Case | Models | Spread (kJ/mol) | Closest diagnostic model | Closest raw deviation |",
                  "|---|---:|---:|---:|---:|"])
    for row in payload["comparison_summary"]["by_case"]:
        spread = "—" if row["cross_model_spread_kj_mol"] is None else f"{row['cross_model_spread_kj_mol']:.1f}"
        closest = row["closest_model_diagnostic"] or "—"
        deviation = "—" if row["closest_raw_deviation_kj_mol"] is None else f"{row['closest_raw_deviation_kj_mol']:+.1f}"
        lines.append(f"| {row['reaction_name']} | {row['models_with_result']} | {spread} | {closest} | {deviation} |")
    lines.extend(["", "## Interpretation", ""])
    lines.extend(f"- {warning}" for warning in payload["warnings"])
    return "\n".join(lines) + "\n"


def write_plots(payload: dict, output: Path) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import SymLogNorm

    engines = [row["engine"] for row in payload["comparison_summary"]["by_engine"]]
    referenced_cases = [case for case in payload["cases"] if case["reference_kj_mol"] is not None]
    deviations = np.full((len(referenced_cases), len(engines)), np.nan)
    lookup = {(row["reaction_id"], row["engine"]): row for row in payload["results"]}
    for case_index, case in enumerate(referenced_cases):
        for engine_index, engine in enumerate(engines):
            value = lookup.get((case["id"], engine), {}).get("raw_deviation_kj_mol")
            if value is not None:
                deviations[case_index, engine_index] = value
    finite = np.abs(deviations[np.isfinite(deviations)])
    extent = max(100.0, float(finite.max(initial=100.0)))
    figure, axis = plt.subplots(figsize=(9.8, 5.8), constrained_layout=True)
    image = axis.imshow(deviations, cmap="coolwarm", norm=SymLogNorm(linthresh=10, vmin=-extent, vmax=extent), aspect="auto")
    axis.set_xticks(range(len(engines)), engines)
    axis.set_yticks(range(len(referenced_cases)), [case["name"] for case in referenced_cases])
    axis.set_title("Raw model ΔE deviation from NIST gas-phase ΔH°298")
    axis.set_xlabel("Engine")
    axis.set_ylabel("Balanced reaction")
    for row in range(deviations.shape[0]):
        for column in range(deviations.shape[1]):
            value = deviations[row, column]
            axis.text(column, row, "N/A" if np.isnan(value) else f"{value:+.0f}", ha="center", va="center", fontsize=8,
                      color="black" if np.isnan(value) or abs(value) < 100 else "white")
    colorbar = figure.colorbar(image, ax=axis, shrink=0.86)
    colorbar.set_label("model ΔE − reference ΔH°298 (kJ/mol; symlog color scale)")
    heatmap = output / "diagnostic_deviation.png"
    figure.savefig(heatmap, dpi=180)
    plt.close(figure)

    case_summary = payload["comparison_summary"]["by_case"]
    figure, axis = plt.subplots(figsize=(9.8, 5.5), constrained_layout=True)
    names = [row["reaction_name"] for row in case_summary]
    spreads = [row["cross_model_spread_kj_mol"] for row in case_summary]
    bars = axis.barh(range(len(names)), spreads, color="#38bdf8")
    axis.set_yticks(range(len(names)), names)
    axis.invert_yaxis()
    axis.set_xlabel("max − min model prediction (kJ/mol)")
    axis.set_title("Cross-model reaction-energy spread")
    axis.set_xlim(0, max(1.0, max(spreads)) * 1.12)
    axis.grid(axis="x", alpha=0.25)
    axis.bar_label(bars, labels=[f"{value:.1f}" for value in spreads], padding=4, fontsize=8)
    spread_plot = output / "cross_model_spread.png"
    figure.savefig(spread_plot, dpi=180)
    plt.close(figure)
    return [heatmap.name, spread_plot.name]


def prediction_only_payload(payload: dict) -> dict:
    """Project computed values into a reference-free artifact before scoring."""
    projected = {
        key: copy.deepcopy(value)
        for key, value in payload.items()
        if key not in {"comparison_summary", "artifacts", "evaluation_protocol"}
    }
    for case in projected.get("cases", []):
        for key in ("reference_kj_mol", "reference_note", "references"):
            case.pop(key, None)
    for row in projected.get("results", []):
        for key in ("reference_kj_mol", "raw_deviation_kj_mol", "sign_agreement"):
            row.pop(key, None)
    return projected


def reference_payload(payload: dict) -> dict:
    return {
        "cases": [
            {
                "id": case["id"],
                "reference_kj_mol": case.get("reference_kj_mol"),
                "reference_note": case.get("reference_note"),
                "references": case.get("references", []),
            }
            for case in payload.get("cases", [])
        ]
    }


def persist(payload: dict, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    frozen = freeze_prediction_artifact(prediction_only_payload(payload), output / "prediction.frozen.json")
    unlocked = unlock_evaluation(
        output / "prediction.frozen.json", reference_payload(payload), output / "evaluation.unlocked.json"
    )
    payload["evaluation_protocol"] = {
        "status": unlocked["status"],
        "prediction_sha256": frozen["prediction_sha256"],
        "reference_sha256": unlocked["reference_sha256"],
        "prediction_artifact": "prediction.frozen.json",
        "evaluation_artifact": "evaluation.unlocked.json",
    }
    add_comparison_summary(payload)
    payload["artifacts"] = write_plots(payload, output)
    (output / "latest.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (output / "latest.md").write_text(markdown_report(payload), encoding="utf-8")
    columns = ["reaction_id", "reaction_name", "category", "engine", "status", "predicted_kj_mol", "reference_kj_mol",
               "raw_deviation_kj_mol", "sign_agreement", "all_converged", "connectivity_preserved", "max_force_ev_a", "wall_seconds"]
    with (output / "latest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(payload["results"])


def main():
    args = parse_args()
    if args.summarize_existing:
        payload = json.loads(args.summarize_existing.read_text(encoding="utf-8"))
        persist(payload, args.output)
        print(json.dumps({"event": "summary_rebuilt", "output": str(args.output)}), flush=True)
        return
    unknown_engines = set(args.engines) - set(ALL_ENGINES)
    selected = [case for case in REACTIONS if case.id in args.cases]
    if unknown_engines or len(selected) != len(set(args.cases)):
        raise SystemExit(f"Unknown engines/cases: engines={sorted(unknown_engines)}, cases={sorted(set(args.cases) - {case.id for case in selected})}")

    payload = {
        "schema_version": "0.1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "device": args.device,
        "requested_engines": args.engines,
        "settings": {"fmax_ev_a": args.fmax, "max_steps": args.steps, "seed": args.seed},
        "cases": [],
        "species_results": [],
        "results": [],
        "warnings": [
            "NIST ΔH°298 and model electronic ΔE are not identical observables; no ZPE, thermal, entropy, solvent or pressure correction was applied.",
            "A shared sign or close numerical agreement is not a calibrated probability or proof of mechanistic accuracy.",
            "MACE, ANI and AIMNet errors are correlated; model agreement does not create independent experimental evidence.",
            "Catalysts primarily change activation barriers. Endpoint energies with an omitted catalyst do not benchmark catalytic activity or selectivity.",
        ],
    }
    for case in selected:
        payload["cases"].append({
            "id": case.id, "name": case.name, "category": case.category, "difficulty": case.difficulty,
            "stoichiometry": case.stoichiometry, "reference_kj_mol": case.reference_kj_mol,
            "reference_note": case.reference_note, "references": case.references, "catalyst_note": case.catalyst_note,
            "eligible_engines": case.eligible_engines,
        })

    geometry_cache = {species_id: build_geometry(SPECIES[species_id], args.seed) for case in selected for species_id in case.stoichiometry}
    structures_root = args.output / "structures"
    for engine_id in args.engines:
        engine_cases = [case for case in selected if engine_id in case.eligible_engines]
        for case in selected:
            if engine_id not in case.eligible_engines:
                payload["results"].append({
                    "reaction_id": case.id, "reaction_name": case.name, "category": case.category, "engine": engine_id,
                    "status": "unsupported", "reason": "outside declared element/model domain", "reference_kj_mol": case.reference_kj_mol,
                })
        if not engine_cases:
            continue
        print(json.dumps({"event": "engine_start", "engine": engine_id, "cases": len(engine_cases)}), flush=True)
        engine = make_engine(engine_id, args.device)
        required_species = sorted({species_id for case in engine_cases for species_id in case.stoichiometry})
        energies: dict[str, float] = {}
        species_rows: dict[str, dict] = {}
        for species_id in required_species:
            atoms, expected_connectivity, formula = geometry_cache[species_id]
            started = time.perf_counter()
            try:
                result = engine.optimize(atoms, fmax=args.fmax, steps=args.steps)
                elapsed = time.perf_counter() - started
                actual_connectivity = inferred_connectivity(result.atoms)
                row = {
                    "engine": engine_id, "species_id": species_id, "species_name": SPECIES[species_id].name,
                    "formula": formula, "energy_ev": result.energy_ev, "converged": result.converged,
                    "steps": result.steps, "max_force_ev_a": result.max_force,
                    "connectivity_preserved": actual_connectivity == expected_connectivity,
                    "expected_bonds": sorted(expected_connectivity), "inferred_bonds": sorted(actual_connectivity),
                    "wall_seconds": elapsed, "status": "ok", "metadata": result.metadata,
                }
                energies[species_id] = result.energy_ev
                species_rows[species_id] = row
                payload["species_results"].append(row)
                destination = structures_root / engine_id
                destination.mkdir(parents=True, exist_ok=True)
                write(destination / f"{species_id}.xyz", result.atoms)
                print(json.dumps({"event": "species_done", "engine": engine_id, "species": species_id,
                                  "converged": result.converged, "connectivity": row["connectivity_preserved"],
                                  "seconds": round(elapsed, 3)}), flush=True)
            except Exception as error:
                elapsed = time.perf_counter() - started
                row = {"engine": engine_id, "species_id": species_id, "species_name": SPECIES[species_id].name,
                       "formula": formula, "status": "failed", "error": f"{type(error).__name__}: {error}", "wall_seconds": elapsed}
                species_rows[species_id] = row
                payload["species_results"].append(row)
                print(json.dumps({"event": "species_failed", **row}), flush=True)

        for case in engine_cases:
            started = time.perf_counter()
            missing = [species_id for species_id in case.stoichiometry if species_id not in energies]
            if missing:
                row = {"reaction_id": case.id, "reaction_name": case.name, "category": case.category, "engine": engine_id,
                       "status": "failed", "reason": f"species failed: {', '.join(missing)}", "reference_kj_mol": case.reference_kj_mol}
            else:
                delta_ev = sum(coefficient * energies[species_id] for species_id, coefficient in case.stoichiometry.items())
                predicted = delta_ev * EV_TO_KJ_MOL
                all_converged = all(species_rows[species_id]["converged"] for species_id in case.stoichiometry)
                connectivity = all(species_rows[species_id]["connectivity_preserved"] for species_id in case.stoichiometry)
                reference = case.reference_kj_mol
                row = {
                    "reaction_id": case.id, "reaction_name": case.name, "category": case.category, "engine": engine_id,
                    "status": "ok" if all_converged and connectivity else "diagnostic_only",
                    "predicted_ev": delta_ev, "predicted_kj_mol": predicted, "reference_kj_mol": reference,
                    "raw_deviation_kj_mol": None if reference is None else predicted - reference,
                    "sign_agreement": None if reference is None else (predicted == 0) or ((predicted < 0) == (reference < 0)),
                    "all_converged": all_converged, "connectivity_preserved": connectivity,
                    "max_force_ev_a": max(species_rows[species_id]["max_force_ev_a"] for species_id in case.stoichiometry),
                    "wall_seconds": sum(species_rows[species_id]["wall_seconds"] for species_id in case.stoichiometry),
                    "catalyst_note": case.catalyst_note,
                }
            payload["results"].append(row)
            print(json.dumps({"event": "reaction_done", **row}), flush=True)
        release_engine(engine)
        persist(payload, args.output)
    persist(payload, args.output)
    print(json.dumps({"event": "complete", "output": str(args.output), "rows": len(payload["results"])}), flush=True)


if __name__ == "__main__":
    main()
