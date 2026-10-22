from copy import deepcopy

from fastapi.testclient import TestClient

from app.main import app
from app.validation.calibration import evaluate_held_out, fit_calibration, freeze_predictions
from app.validation.schema import (
    CalibrationRequest,
    HeldOutEvaluationRequest,
    PredictionSetRequest,
)


def prediction_payload() -> dict:
    return {
        "artifact_id": "campaign-001",
        "model_id": "frozen-model-sha256:example",
        "output_name": "yield_fraction",
        "unit": "dimensionless",
        "predictions": [
            {"condition_id": "cal-1", "mean": 0.20, "standard_deviation": 0.04},
            {"condition_id": "cal-2", "mean": 0.40, "standard_deviation": 0.04},
            {"condition_id": "cal-3", "mean": 0.60, "standard_deviation": 0.05},
            {"condition_id": "test-1", "mean": 0.30, "standard_deviation": 0.04},
            {"condition_id": "test-2", "mean": 0.50, "standard_deviation": 0.05},
        ],
    }


def fitted_artifacts() -> tuple[dict, dict]:
    frozen = freeze_predictions(PredictionSetRequest.model_validate(prediction_payload()))
    calibration = fit_calibration(CalibrationRequest.model_validate({
        "frozen": frozen,
        "method": "affine",
        "observations": [
            {"condition_id": "cal-1", "value": 0.28},
            {"condition_id": "cal-2", "value": 0.49},
            {"condition_id": "cal-3", "value": 0.68},
        ],
    }))
    return frozen, calibration


def test_freeze_calibrate_and_heldout_evaluate_without_leakage() -> None:
    frozen, calibration = fitted_artifacts()
    result = evaluate_held_out(HeldOutEvaluationRequest.model_validate({
        "frozen": frozen,
        "calibration": calibration,
        "observations": [
            {"condition_id": "test-1", "value": 0.38},
            {"condition_id": "test-2", "value": 0.59},
        ],
    }))
    assert result["leakage_check"] == "passed"
    assert result["calibrated_metrics"]["count"] == 2
    assert result["calibrated_metrics"]["mae"] < result["raw_metrics"]["mae"]


def test_tampered_prediction_and_training_overlap_are_rejected() -> None:
    frozen, calibration = fitted_artifacts()
    tampered = deepcopy(frozen)
    tampered["predictions"][0]["mean"] = 0.99
    response = TestClient(app).post("/api/validation/evaluate", json={
        "frozen": tampered,
        "observations": [{"condition_id": "test-1", "value": 0.38}],
    })
    assert response.status_code == 422
    overlap = TestClient(app).post("/api/validation/evaluate", json={
        "frozen": frozen,
        "calibration": calibration,
        "observations": [{"condition_id": "cal-1", "value": 0.28}],
    })
    assert overlap.status_code == 422
    assert "overlaps calibration ids" in overlap.text


def test_validation_endpoints() -> None:
    client = TestClient(app)
    frozen_response = client.post("/api/validation/freeze", json=prediction_payload())
    assert frozen_response.status_code == 200
    frozen = frozen_response.json()
    calibrated = client.post("/api/validation/calibrate", json={
        "frozen": frozen,
        "method": "offset_only",
        "observations": [
            {"condition_id": "cal-1", "value": 0.28},
            {"condition_id": "cal-2", "value": 0.49},
        ],
    })
    assert calibrated.status_code == 200, calibrated.text
    assert len(calibrated.json()["artifact_sha256"]) == 64
