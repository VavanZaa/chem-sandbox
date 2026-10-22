from __future__ import annotations

from typing import Literal

from pydantic import Field, model_validator

from .schema import OperationValidationRequest, StrictModel


class CanteraMechanismSource(StrictModel):
    builtin: Literal["gri30.yaml", "h2o2.yaml"] | None = None
    yaml_text: str | None = Field(default=None, min_length=20, max_length=5_000_000)
    phase_name: str | None = Field(default=None, max_length=120)
    source: str = Field(min_length=3, max_length=1000)
    allow_declared_subset_of_cantera_mechanism: Literal[True]

    @model_validator(mode="after")
    def exactly_one_source(self) -> "CanteraMechanismSource":
        if (self.builtin is None) == (self.yaml_text is None):
            raise ValueError("Supply exactly one Cantera builtin mechanism or YAML text")
        return self


class CanteraSimulationRequest(StrictModel):
    definition: OperationValidationRequest
    cantera: CanteraMechanismSource
    batch_constraint: Literal["constant_volume", "constant_pressure"] = "constant_volume"
    points: int = Field(default=101, ge=2, le=501)
    reference_species: str | None = Field(default=None, max_length=80)
    target_species: list[str] = Field(default_factory=list, max_length=20)
    relative_tolerance: float = Field(default=1e-8, gt=0, le=1e-3)
    absolute_tolerance: float = Field(default=1e-15, gt=0, le=1e-6)
    inventory_state_tolerance: float = Field(default=0.05, gt=0, le=0.5)
