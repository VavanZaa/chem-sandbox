from __future__ import annotations

import re
from typing import Literal

from pydantic import Field, model_validator

from app.operations.schema import QuantityInput, StrictModel
from app.thermochemistry.schema import EvidenceQuantity


SPECIES_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.+-]{0,63}$")


class CatalyticSpecies(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    name: str = Field(min_length=1, max_length=120)
    phase: Literal["gas", "surface"]
    elements: dict[str, int] = Field(default_factory=dict)
    charge: int = Field(default=0, ge=-20, le=20)
    site_occupancy: int = Field(default=0, ge=0, le=8)
    vacant_site: bool = False

    @model_validator(mode="after")
    def validate_species(self) -> "CatalyticSpecies":
        if any(not re.fullmatch(r"[A-Z][a-z]?", key) or value <= 0 for key, value in self.elements.items()):
            raise ValueError("Element keys must be symbols and counts must be positive integers")
        if self.phase == "gas" and (self.site_occupancy or self.vacant_site):
            raise ValueError("Gas species cannot occupy or represent surface sites")
        if self.phase == "surface" and self.site_occupancy < 1:
            raise ValueError("Surface species require a positive site_occupancy")
        if self.vacant_site and (self.elements or self.charge):
            raise ValueError("The vacant-site pseudo-species cannot carry atoms or charge")
        if not self.vacant_site and not self.elements:
            raise ValueError("Non-vacant species require an elemental composition")
        return self


class SurfaceElementaryReaction(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    process: Literal["adsorption", "desorption", "surface_reaction", "diffusion"]
    stoichiometry: dict[str, float] = Field(min_length=2, max_length=40)
    forward_activation_gibbs: EvidenceQuantity
    reverse_activation_gibbs: EvidenceQuantity
    reaction_gibbs: EvidenceQuantity
    transmission_coefficient: float = Field(default=1.0, gt=0, le=1)

    @model_validator(mode="after")
    def validate_coefficients(self) -> "SurfaceElementaryReaction":
        values = list(self.stoichiometry.values())
        if any(abs(value) <= 1e-15 for value in values):
            raise ValueError("Reaction stoichiometric coefficients cannot be zero")
        if not any(value < 0 for value in values) or not any(value > 0 for value in values):
            raise ValueError("An elementary reaction requires reactants and products")
        return self


class FirstOrderDeactivation(StrictModel):
    model: Literal["none", "first_order"] = "none"
    rate_constant: EvidenceQuantity | None = None

    @model_validator(mode="after")
    def validate_model(self) -> "FirstOrderDeactivation":
        if self.model == "first_order" and self.rate_constant is None:
            raise ValueError("First-order deactivation requires a sourced rate constant")
        if self.model == "none" and self.rate_constant is not None:
            raise ValueError("A deactivation rate is only valid for first_order")
        return self


class TransportEffectiveness(StrictModel):
    value: float = Field(gt=0, le=1)
    source: str = Field(min_length=3, max_length=500)
    basis: Literal["computed", "experimental", "literature", "hypothetical"]


class CatalysisSimulationRequest(StrictModel):
    id: str = Field(min_length=1, max_length=80)
    temperature: QuantityInput
    pressure: QuantityInput
    duration: QuantityInput
    points: int = Field(default=101, ge=2, le=1001)
    active_site_amount: QuantityInput
    species: list[CatalyticSpecies] = Field(min_length=3, max_length=200)
    reactions: list[SurfaceElementaryReaction] = Field(min_length=1, max_length=500)
    initial_surface_amounts_per_site: dict[str, float] = Field(default_factory=dict)
    gas_activities: dict[str, float] = Field(default_factory=dict)
    deactivation: FirstOrderDeactivation = Field(default_factory=FirstOrderDeactivation)
    transport_effectiveness: TransportEffectiveness | None = None
    target_gas_species: list[str] = Field(default_factory=list, max_length=20)
    detailed_balance_tolerance: QuantityInput = Field(
        default_factory=lambda: QuantityInput(value=1.0, unit="kJ/mol")
    )

    @model_validator(mode="after")
    def validate_network(self) -> "CatalysisSimulationRequest":
        ids = [item.id for item in self.species]
        if len(ids) != len(set(ids)):
            raise ValueError("Catalytic species ids must be unique")
        lookup = {item.id: item for item in self.species}
        vacant = [item.id for item in self.species if item.vacant_site]
        if len(vacant) != 1:
            raise ValueError("Exactly one explicit vacant-site pseudo-species is required")
        surface_ids = {item.id for item in self.species if item.phase == "surface"}
        gas_ids = set(ids) - surface_ids
        if set(self.initial_surface_amounts_per_site) - surface_ids:
            raise ValueError("Initial surface amounts reference a non-surface species")
        if set(self.gas_activities) - gas_ids:
            raise ValueError("Gas activities reference a non-gas species")
        if set(self.target_gas_species) - gas_ids:
            raise ValueError("Target gas species must be declared gas species")
        if any(value < 0 for value in self.initial_surface_amounts_per_site.values()):
            raise ValueError("Initial surface amounts cannot be negative")
        if any(value < 0 for value in self.gas_activities.values()):
            raise ValueError("Gas activities cannot be negative")
        occupied = sum(lookup[key].site_occupancy * value for key, value in self.initial_surface_amounts_per_site.items())
        if occupied > 1 + 1e-10:
            raise ValueError("Initial surface site balance exceeds one total site")
        known = set(ids)
        reaction_ids = [reaction.id for reaction in self.reactions]
        if len(reaction_ids) != len(set(reaction_ids)):
            raise ValueError("Catalytic reaction ids must be unique")
        for reaction in self.reactions:
            unknown = set(reaction.stoichiometry) - known
            if unknown:
                raise ValueError(f"Reaction {reaction.id} references unknown species: {sorted(unknown)}")
            changed_surface = [
                key for key in reaction.stoichiometry
                if key in surface_ids and not lookup[key].vacant_site
            ]
            if not changed_surface:
                raise ValueError(
                    f"Reaction {reaction.id} has no bound intermediate; catalyst-free endpoint cancellation is forbidden"
                )
            element_delta: dict[str, float] = {}
            charge_delta = 0.0
            site_delta = 0.0
            for key, coefficient in reaction.stoichiometry.items():
                item = lookup[key]
                for element, count in item.elements.items():
                    element_delta[element] = element_delta.get(element, 0.0) + coefficient * count
                charge_delta += coefficient * item.charge
                if item.phase == "surface":
                    site_delta += coefficient * item.site_occupancy
            residual = {key: value for key, value in element_delta.items() if abs(value) > 1e-10}
            if residual or abs(charge_delta) > 1e-10 or abs(site_delta) > 1e-10:
                raise ValueError(
                    f"Reaction {reaction.id} is not element/charge/site balanced: "
                    f"elements={residual}, charge={charge_delta:g}, sites={site_delta:g}"
                )
        return self

