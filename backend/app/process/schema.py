from __future__ import annotations

import re
from typing import Literal

from pydantic import Field, model_validator

from app.operations.schema import QuantityInput, StrictModel
from app.thermochemistry.schema import EvidenceQuantity


COMPONENT_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.+-]{0,63}$")


class ProcessComponent(StrictModel):
    id: str = Field(pattern=COMPONENT_ID.pattern)
    name: str = Field(min_length=1, max_length=120)
    elements: dict[str, int] = Field(min_length=1)
    charge: int = Field(default=0, ge=-20, le=20)

    @model_validator(mode="after")
    def validate_elements(self) -> "ProcessComponent":
        if any(not re.fullmatch(r"[A-Z][a-z]?", key) or value <= 0 for key, value in self.elements.items()):
            raise ValueError("Element keys must be symbols and counts must be positive integers")
        return self


class ProcessReaction(StrictModel):
    id: str = Field(min_length=1, max_length=80)
    stoichiometry: dict[str, float] = Field(min_length=2, max_length=100)
    limiting_species: str = Field(pattern=COMPONENT_ID.pattern)
    evidence_source: str = Field(min_length=3, max_length=500)

    @model_validator(mode="after")
    def validate_reaction(self) -> "ProcessReaction":
        values = list(self.stoichiometry.values())
        if any(abs(value) <= 1e-15 for value in values):
            raise ValueError("Reaction coefficients cannot be zero")
        if not any(value < 0 for value in values) or not any(value > 0 for value in values):
            raise ValueError("Reaction requires reactants and products")
        if self.stoichiometry.get(self.limiting_species, 0.0) >= 0:
            raise ValueError("limiting_species must be a reactant")
        return self


class RecyclePolicy(StrictModel):
    fraction: EvidenceQuantity
    optimize: bool = False
    lower_bound: float = Field(default=0.0, ge=0, lt=1)
    upper_bound: float = Field(default=0.95, gt=0, lt=1)
    grid_points: int = Field(default=101, ge=11, le=1001)

    @model_validator(mode="after")
    def validate_bounds(self) -> "RecyclePolicy":
        if self.lower_bound >= self.upper_bound:
            raise ValueError("Recycle optimization lower_bound must be below upper_bound")
        return self


class ProcessFlowsheetRequest(StrictModel):
    id: str = Field(min_length=1, max_length=80)
    components: list[ProcessComponent] = Field(min_length=2, max_length=200)
    fresh_feed: dict[str, QuantityInput] = Field(min_length=1)
    reaction: ProcessReaction
    single_pass_conversion: EvidenceQuantity
    product_recoveries: dict[str, EvidenceQuantity]
    recycle: RecyclePolicy
    target_species: str = Field(pattern=COMPONENT_ID.pattern)
    utility_per_extent: EvidenceQuantity | None = None
    max_internal_molar_flow: QuantityInput | None = None
    convergence_tolerance: float = Field(default=1e-10, gt=0, le=1e-4)
    max_iterations: int = Field(default=1000, ge=10, le=100000)

    @model_validator(mode="after")
    def validate_cross_references(self) -> "ProcessFlowsheetRequest":
        ids = [item.id for item in self.components]
        if len(ids) != len(set(ids)):
            raise ValueError("Process component ids must be unique")
        known = set(ids)
        if set(self.fresh_feed) - known:
            raise ValueError("Fresh feed references unknown components")
        if set(self.reaction.stoichiometry) - known:
            raise ValueError("Reaction references unknown components")
        if set(self.product_recoveries) != known:
            raise ValueError("A sourced product recovery is required for every component")
        if self.target_species not in known:
            raise ValueError("target_species is not declared")
        lookup = {item.id: item for item in self.components}
        element_delta: dict[str, float] = {}
        charge_delta = 0.0
        for key, coefficient in self.reaction.stoichiometry.items():
            item = lookup[key]
            for element, count in item.elements.items():
                element_delta[element] = element_delta.get(element, 0.0) + coefficient * count
            charge_delta += coefficient * item.charge
        residual = {key: value for key, value in element_delta.items() if abs(value) > 1e-10}
        if residual or abs(charge_delta) > 1e-10:
            raise ValueError(f"Process reaction is not conserved: elements={residual}, charge={charge_delta:g}")
        return self

