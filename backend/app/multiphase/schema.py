from __future__ import annotations

from typing import Literal

from pydantic import Field, model_validator

from app.operations.schema import QuantityInput, StrictModel
from app.thermochemistry.schema import EvidenceQuantity


class GasLiquidTransferRequest(StrictModel):
    solute_id: str = Field(min_length=1, max_length=80)
    mode: Literal["closed_headspace", "open_fixed_partial_pressure"]
    temperature: QuantityInput
    total_pressure: QuantityInput
    gas_volume: QuantityInput
    liquid_volume: QuantityInput
    duration: QuantityInput
    points: int = Field(default=101, ge=2, le=1001)
    initial_gas_amount: QuantityInput
    initial_dissolved_amount: QuantityInput
    henry_hcp: EvidenceQuantity
    kla: EvidenceQuantity
    ambient_partial_pressure: QuantityInput | None = None

    @model_validator(mode="after")
    def validate_mode(self) -> "GasLiquidTransferRequest":
        if self.mode == "open_fixed_partial_pressure" and self.ambient_partial_pressure is None:
            raise ValueError("Open gas transfer requires ambient_partial_pressure")
        if self.mode == "closed_headspace" and self.ambient_partial_pressure is not None:
            raise ValueError("Closed headspace does not accept an ambient partial pressure")
        return self


class ElectrolyteComponent(StrictModel):
    species: str = Field(min_length=1, max_length=120)
    amount: QuantityInput


class ElectrolyteEquilibriumRequest(StrictModel):
    database_kind: Literal["phreeqc", "supcrt", "nasa"]
    database_name: str = Field(min_length=2, max_length=200)
    aqueous_species: list[str] = Field(min_length=1, max_length=500)
    gas_species: list[str] = Field(default_factory=list, max_length=200)
    mineral_species: list[str] = Field(default_factory=list, max_length=200)
    activity_model: Literal["ideal", "davies", "debye_huckel", "pitzer", "phreeqc"]
    temperature: QuantityInput
    pressure: QuantityInput
    components: list[ElectrolyteComponent] = Field(min_length=1, max_length=500)
    fixed_pH: float | None = Field(default=None, ge=-2, le=16)
    fixed_gas_fugacity: dict[str, QuantityInput] = Field(default_factory=dict, max_length=20)
    requested_species: list[str] = Field(default_factory=list, max_length=500)
    source: str = Field(min_length=3, max_length=1000)
