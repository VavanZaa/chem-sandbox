from __future__ import annotations

from typing import Any, Literal

from pydantic import Field, model_validator

from app.operations.schema import StrictModel


class PredictionRecord(StrictModel):
    condition_id: str = Field(min_length=1, max_length=120)
    mean: float
    standard_deviation: float = Field(gt=0)
    metadata: dict[str, Any] = Field(default_factory=dict)


class PredictionSetRequest(StrictModel):
    artifact_id: str = Field(min_length=1, max_length=120)
    model_id: str = Field(min_length=1, max_length=200)
    output_name: str = Field(min_length=1, max_length=120)
    unit: str = Field(min_length=1, max_length=40)
    predictions: list[PredictionRecord] = Field(min_length=1, max_length=100000)

    @model_validator(mode="after")
    def unique_conditions(self) -> "PredictionSetRequest":
        ids = [item.condition_id for item in self.predictions]
        if len(ids) != len(set(ids)):
            raise ValueError("Prediction condition ids must be unique")
        return self


class FrozenPredictionSet(PredictionSetRequest):
    protocol: Literal["post-run-calibration-v1"]
    status: Literal["prediction_frozen"]
    frozen_at: str
    prediction_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class Observation(StrictModel):
    condition_id: str = Field(min_length=1, max_length=120)
    value: float


class CalibrationRequest(StrictModel):
    frozen: FrozenPredictionSet
    observations: list[Observation] = Field(min_length=2, max_length=100000)
    method: Literal["affine", "offset_only"] = "affine"

    @model_validator(mode="after")
    def unique_observations(self) -> "CalibrationRequest":
        ids = [item.condition_id for item in self.observations]
        if len(ids) != len(set(ids)):
            raise ValueError("Calibration observation ids must be unique")
        return self


class CalibrationArtifact(StrictModel):
    protocol: Literal["post-run-calibration-v1"]
    status: Literal["calibration_fitted"]
    prediction_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    model_id: str
    output_name: str
    unit: str
    method: Literal["affine", "offset_only"]
    intercept: float
    slope: float
    uncertainty_scale: float = Field(gt=0)
    training_ids: list[str]
    observation_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    artifact_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class HeldOutEvaluationRequest(StrictModel):
    frozen: FrozenPredictionSet
    observations: list[Observation] = Field(min_length=1, max_length=100000)
    calibration: CalibrationArtifact | None = None

    @model_validator(mode="after")
    def unique_observations(self) -> "HeldOutEvaluationRequest":
        ids = [item.condition_id for item in self.observations]
        if len(ids) != len(set(ids)):
            raise ValueError("Held-out observation ids must be unique")
        return self

