from __future__ import annotations

import hashlib
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any


REFERENCE_KEYS = frozenset(
    {
        "answer",
        "expected_pathway",
        "expected_product",
        "known_yield",
        "literature_yield",
        "raw_deviation_kj_mol",
        "reference_kj_mol",
        "reference_note",
        "references",
        "sign_agreement",
        "target_product",
    }
)


def canonical_sha256(payload: Any) -> str:
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _find_reference_keys(value: Any, prefix: str = "$") -> list[str]:
    hits: list[str] = []
    if isinstance(value, dict):
        for key, nested in value.items():
            path = f"{prefix}.{key}"
            if key in REFERENCE_KEYS:
                hits.append(path)
            hits.extend(_find_reference_keys(nested, path))
    elif isinstance(value, list):
        for index, nested in enumerate(value):
            hits.extend(_find_reference_keys(nested, f"{prefix}[{index}]"))
    return hits


def freeze_prediction_artifact(payload: dict[str, Any], destination: Path) -> dict[str, Any]:
    """Write an immutable pre-reference artifact for later evaluation.

    The caller must pass prediction-only content. Any reference-shaped field is
    rejected so a benchmark cannot accidentally hash an already-unlocked file.
    """

    hits = _find_reference_keys(payload)
    if hits:
        raise ValueError(f"Reference fields are forbidden before prediction freeze: {', '.join(hits[:8])}")
    frozen = {
        "protocol": "evaluation-unlock-v1",
        "status": "prediction_frozen",
        "frozen_at": datetime.now(UTC).isoformat(),
        "prediction": payload,
    }
    frozen["prediction_sha256"] = canonical_sha256(payload)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(frozen, indent=2, ensure_ascii=False), encoding="utf-8")
    return frozen


def unlock_evaluation(
    frozen_artifact: Path,
    references: dict[str, Any],
    destination: Path,
) -> dict[str, Any]:
    frozen = json.loads(frozen_artifact.read_text(encoding="utf-8"))
    prediction = frozen.get("prediction")
    expected_hash = frozen.get("prediction_sha256")
    actual_hash = canonical_sha256(prediction)
    if frozen.get("status") != "prediction_frozen" or expected_hash != actual_hash:
        raise ValueError("Frozen prediction artifact is missing, altered, or invalid")
    unlocked = {
        "protocol": "evaluation-unlock-v1",
        "status": "references_unlocked_post_run",
        "unlocked_at": datetime.now(UTC).isoformat(),
        "prediction_sha256": expected_hash,
        "reference_sha256": canonical_sha256(references),
        "prediction": prediction,
        "references": references,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(unlocked, indent=2, ensure_ascii=False), encoding="utf-8")
    return unlocked
