from __future__ import annotations

from typing import Literal

from pydantic import Field, model_validator

from app.operations.schema import QuantityInput, StrictModel
from app.thermochemistry.schema import EvidenceQuantity


class PressureRatePoint(StrictModel):
    pressure: QuantityInput
    forward_multiplier: float = Field(gt=0, le=1e12)
    reverse_multiplier: float = Field(gt=0, le=1e12)


class PressureDependence(StrictModel):
    model: Literal["none", "tabulated_log_interpolation"] = "none"
    points: list[PressureRatePoint] = Field(default_factory=list, max_length=100)
    source: str | None = Field(default=None, max_length=500)
    basis: Literal["computed", "experimental", "literature", "hypothetical"] | None = None

    @model_validator(mode="after")
    def validate_pressure_data(self) -> "PressureDependence":
        if self.model == "none" and self.points:
            raise ValueError("Pressure points require tabulated_log_interpolation")
        if self.model != "none":
            if len(self.points) < 2 or not self.source or self.basis is None:
                raise ValueError("Pressure dependence requires at least two sourced points")
            pressures = [item.pressure.as_unit("Pa") for item in self.points]
            if any(value <= 0 for value in pressures) or len(set(pressures)) != len(pressures):
                raise ValueError("Pressure points must be positive and unique")
        return self


class RateArtifactRequest(StrictModel):
    reaction_id: str = Field(min_length=1, max_length=80)
    temperature: QuantityInput
    pressure: QuantityInput | None = None
    forward_order: int = Field(ge=1, le=3)
    reverse_order: int = Field(ge=1, le=3)
    standard_concentration: QuantityInput = Field(
        default_factory=lambda: QuantityInput(value=1000.0, unit="mol/m^3")
    )
    forward_activation_gibbs: EvidenceQuantity
    reverse_activation_gibbs: EvidenceQuantity
    reaction_gibbs: EvidenceQuantity
    transmission_coefficient: float = Field(default=1.0, gt=0, le=1.0)
    detailed_balance_tolerance: QuantityInput = Field(
        default_factory=lambda: QuantityInput(value=1.0, unit="kJ/mol")
    )
    pressure_dependence: PressureDependence = Field(default_factory=PressureDependence)
    impact_weight: float = Field(default=1.0, ge=0, le=1)
