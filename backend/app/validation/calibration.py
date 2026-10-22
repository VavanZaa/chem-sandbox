from __future__ import annotations

import hashlib
import json
import math
from datetime import UTC, datetime
from typing import Any

import numpy as np

from .schema import CalibrationRequest, HeldOutEvaluationRequest, PredictionSetRequest


def _canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _prediction_payload(frozen) -> dict[str, Any]:
    return {
        "artifact_id": frozen.artifact_id,
        "model_id": frozen.model_id,
        "output_name": frozen.output_name,
        "unit": frozen.unit,
        "predictions": [item.model_dump(mode="json") for item in frozen.predictions],
    }


def _verify_frozen(frozen) -> None:
    if _canonical_hash(_prediction_payload(frozen)) != frozen.prediction_sha256:
        raise ValueError("Frozen prediction payload was altered after hashing")


def freeze_predictions(request: PredictionSetRequest) -> dict[str, Any]:
    payload = request.model_dump(mode="json")
    return {
        **payload,
        "protocol": "post-run-calibration-v1",
        "status": "prediction_frozen",
        "frozen_at": datetime.now(UTC).isoformat(),
        "prediction_sha256": _canonical_hash(payload),
    }


def fit_calibration(request: CalibrationRequest) -> dict[str, Any]:
    _verify_frozen(request.frozen)
    predictions = {item.condition_id: item for item in request.frozen.predictions}
    unknown = [item.condition_id for item in request.observations if item.condition_id not in predictions]
    if unknown:
        raise ValueError(f"Calibration observations have no frozen prediction: {unknown[:8]}")
    x = np.asarray([predictions[item.condition_id].mean for item in request.observations], dtype=float)
    y = np.asarray([item.value for item in request.observations], dtype=float)
    sigma = np.asarray(
        [predictions[item.condition_id].standard_deviation for item in request.observations], dtype=float
    )
    weights = 1.0 / np.square(sigma)
    if request.method == "offset_only":
        slope = 1.0
        intercept = float(np.average(y - x, weights=weights))
    else:
        design = np.column_stack([np.ones_like(x), x])
        weighted_design = design * np.sqrt(weights)[:, None]
        weighted_y = y * np.sqrt(weights)
        coefficients, *_ = np.linalg.lstsq(weighted_design, weighted_y, rcond=None)
        intercept, slope = map(float, coefficients)
    residual = y - (intercept + slope * x)
    calibrated_sigma = np.maximum(abs(slope) * sigma, np.finfo(float).tiny)
    uncertainty_scale = max(float(np.sqrt(np.mean(np.square(residual / calibrated_sigma)))), 1e-6)
    training_ids = [item.condition_id for item in request.observations]
    observations_payload = [item.model_dump(mode="json") for item in request.observations]
    artifact = {
        "protocol": "post-run-calibration-v1",
        "status": "calibration_fitted",
        "prediction_sha256": request.frozen.prediction_sha256,
        "model_id": request.frozen.model_id,
        "output_name": request.frozen.output_name,
        "unit": request.frozen.unit,
        "method": request.method,
        "intercept": intercept,
        "slope": slope,
        "uncertainty_scale": uncertainty_scale,
        "training_ids": training_ids,
        "observation_sha256": _canonical_hash(observations_payload),
    }
    artifact["artifact_sha256"] = _canonical_hash(artifact)
    return artifact


def _verify_calibration(calibration) -> None:
    payload = calibration.model_dump(mode="json")
    digest = payload.pop("artifact_sha256")
    if _canonical_hash(payload) != digest:
        raise ValueError("Calibration artifact was altered after fitting")


def _metrics(rows: list[dict[str, float]]) -> dict[str, Any]:
    errors = np.asarray([row["mean"] - row["observed"] for row in rows], dtype=float)
    standardized = np.asarray(
        [(row["observed"] - row["mean"]) / row["standard_deviation"] for row in rows], dtype=float
    )
    coverages = {}
    for label, z in (("50", 0.6744897501960817), ("80", 1.2815515655446004), ("95", 1.959963984540054)):
        coverages[label] = float(np.mean(np.abs(standardized) <= z))
    nll = 0.5 * np.mean(
        np.log(2 * math.pi * np.square([row["standard_deviation"] for row in rows]))
        + np.square(standardized)
    )
    return {
        "count": len(rows),
        "mae": float(np.mean(np.abs(errors))),
        "rmse": float(np.sqrt(np.mean(np.square(errors)))),
        "bias": float(np.mean(errors)),
        "gaussian_mean_nll": float(nll),
        "interval_coverage": coverages,
        "mean_standardized_residual": float(np.mean(standardized)),
        "rms_standardized_residual": float(np.sqrt(np.mean(np.square(standardized)))),
    }


def evaluate_held_out(request: HeldOutEvaluationRequest) -> dict[str, Any]:
    _verify_frozen(request.frozen)
    predictions = {item.condition_id: item for item in request.frozen.predictions}
    unknown = [item.condition_id for item in request.observations if item.condition_id not in predictions]
    if unknown:
        raise ValueError(f"Held-out observations have no frozen prediction: {unknown[:8]}")
    if request.calibration is not None:
        _verify_calibration(request.calibration)
        if request.calibration.prediction_sha256 != request.frozen.prediction_sha256:
            raise ValueError("Calibration artifact belongs to a different frozen prediction set")
        overlap = set(request.calibration.training_ids) & {item.condition_id for item in request.observations}
        if overlap:
            raise ValueError(f"Held-out evaluation overlaps calibration ids: {sorted(overlap)[:8]}")
    raw_rows = []
    calibrated_rows = []
    for observation in request.observations:
        prediction = predictions[observation.condition_id]
        raw_rows.append({
            "mean": prediction.mean,
            "standard_deviation": prediction.standard_deviation,
            "observed": observation.value,
        })
        if request.calibration is not None:
            calibrated_rows.append({
                "mean": request.calibration.intercept + request.calibration.slope * prediction.mean,
                "standard_deviation": max(
                    abs(request.calibration.slope)
                    * request.calibration.uncertainty_scale
                    * prediction.standard_deviation,
                    np.finfo(float).tiny,
                ),
                "observed": observation.value,
            })
    observation_payload = [item.model_dump(mode="json") for item in request.observations]
    return {
        "status": "held_out_evaluation_completed",
        "prediction_sha256": request.frozen.prediction_sha256,
        "observation_sha256": _canonical_hash(observation_payload),
        "raw_metrics": _metrics(raw_rows),
        "calibrated_metrics": None if request.calibration is None else _metrics(calibrated_rows),
        "held_out_ids": [item.condition_id for item in request.observations],
        "leakage_check": "passed",
        "boundary": "Observations were supplied only after the prediction artifact was frozen; they do not enter discovery or mechanism generation.",
    }
