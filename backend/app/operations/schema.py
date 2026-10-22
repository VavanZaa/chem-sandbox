from __future__ import annotations

import math
import re
from typing import Literal

import pint
from pydantic import BaseModel, ConfigDict, Field, model_validator


UNIT_REGISTRY = pint.UnitRegistry(autoconvert_offset_to_baseunit=True)
SPECIES_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.+-]{0,63}$")


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)


class QuantityInput(StrictModel):
    value: float
    unit: str = Field(min_length=1, max_length=40)

    def as_unit(self, target: str) -> float:
        try:
            converted = UNIT_REGISTRY.Quantity(self.value, self.unit).to(target)
        except (pint.DimensionalityError, pint.UndefinedUnitError) as error:
            raise ValueError(f"{self.value} {self.unit} is not compatible with {target}") from error
        magnitude = float(converted.magnitude)
        if not math.isfinite(magnitude):
            raise ValueError("Converted quantity must be finite")
        return magnitude


class MechanismSpecies(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    name: str = Field(min_length=1, max_length=120)
    elements: dict[str, int] = Field(min_length=1)
    charge: int = Field(default=0, ge=-20, le=20)
    surface_sites: float = Field(default=0, ge=0)

    @model_validator(mode="after")
    def validate_elements(self) -> "MechanismSpecies":
        if any(not re.fullmatch(r"[A-Z][a-z]?", symbol) or count <= 0 for symbol, count in self.elements.items()):
            raise ValueError("Element keys must be symbols and counts must be positive integers")
        return self


class ElementaryReaction(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    stoichiometry: dict[str, float] = Field(min_length=2)
    reversible: bool = True
    evidence_source: str = Field(min_length=3, max_length=500)

    @model_validator(mode="after")
    def validate_coefficients(self) -> "ElementaryReaction":
        if any(abs(value) <= 1e-15 for value in self.stoichiometry.values()):
            raise ValueError("Reaction stoichiometric coefficients cannot be zero")
        if not any(value < 0 for value in self.stoichiometry.values()) or not any(
            value > 0 for value in self.stoichiometry.values()
        ):
            raise ValueError("A reaction requires both reactants and products")
        return self


class MechanismDefinition(StrictModel):
    id: str = Field(min_length=1, max_length=80)
    species: list[MechanismSpecies] = Field(min_length=1)
    reactions: list[ElementaryReaction] = Field(min_length=1)

    @model_validator(mode="after")
    def validate_identifiers(self) -> "MechanismDefinition":
        species_ids = [item.id for item in self.species]
        reaction_ids = [item.id for item in self.reactions]
        if len(species_ids) != len(set(species_ids)):
            raise ValueError("Mechanism species ids must be unique")
        if len(reaction_ids) != len(set(reaction_ids)):
            raise ValueError("Reaction ids must be unique")
        known = set(species_ids)
        for reaction in self.reactions:
            unknown = set(reaction.stoichiometry) - known
            if unknown:
                raise ValueError(f"Reaction {reaction.id} references unknown species: {sorted(unknown)}")
        return self


class PhaseDefinition(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    kind: Literal["gas", "ideal_solution", "aqueous_electrolyte", "liquid", "solid", "surface"]
    species_ids: list[str] = Field(min_length=1)
    activity_model: str = Field(min_length=2, max_length=120)


class InitialAmount(StrictModel):
    species_id: str = Field(pattern=SPECIES_ID.pattern)
    phase_id: str = Field(pattern=SPECIES_ID.pattern)
    amount: QuantityInput


class SpeciesMolarFlow(StrictModel):
    species_id: str = Field(pattern=SPECIES_ID.pattern)
    flow: QuantityInput


class FeedStream(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    phase_id: str = Field(pattern=SPECIES_ID.pattern)
    species_flows: list[SpeciesMolarFlow] = Field(min_length=1)
    temperature: QuantityInput
    pressure: QuantityInput
    start_time: QuantityInput = Field(default_factory=lambda: QuantityInput(value=0, unit="s"))
    end_time: QuantityInput | None = None


class HeatTransfer(StrictModel):
    mode: Literal["adiabatic", "isothermal", "overall_UA"] = "isothermal"
    setpoint_temperature: QuantityInput | None = None
    ambient_temperature: QuantityInput | None = None
    overall_UA: QuantityInput | None = None

    @model_validator(mode="after")
    def required_heat_data(self) -> "HeatTransfer":
        if self.mode == "isothermal" and self.setpoint_temperature is None:
            raise ValueError("Isothermal operation requires setpoint_temperature")
        if self.mode == "overall_UA" and (self.ambient_temperature is None or self.overall_UA is None):
            raise ValueError("overall_UA mode requires ambient_temperature and overall_UA")
        return self


class MixingModel(StrictModel):
    model: Literal["ideal", "stirred", "axial_dispersion"] = "ideal"
    speed: QuantityInput | None = None
    dispersion_coefficient: QuantityInput | None = None


class GasTransferModel(StrictModel):
    model: Literal["none", "equilibrium", "kla"] = "none"
    kla: QuantityInput | None = None
    equilibrium_model: str | None = Field(default=None, max_length=120)

    @model_validator(mode="after")
    def required_transfer_data(self) -> "GasTransferModel":
        if self.model == "kla" and (self.kla is None or not self.equilibrium_model):
            raise ValueError("kLa gas transfer requires kla and equilibrium_model")
        return self


class CatalystDefinition(StrictModel):
    mode: Literal["none", "homogeneous", "surface"] = "none"
    species_id: str | None = None
    loading: QuantityInput | None = None
    active_site_amount: QuantityInput | None = None

    @model_validator(mode="after")
    def required_catalyst_data(self) -> "CatalystDefinition":
        if self.mode != "none" and not self.species_id:
            raise ValueError("Catalytic operation requires an explicit catalyst species_id")
        if self.mode == "homogeneous" and self.loading is None:
            raise ValueError("Homogeneous catalyst requires loading")
        if self.mode == "surface" and self.active_site_amount is None:
            raise ValueError("Surface catalyst requires active_site_amount")
        return self


class ReactorOperation(StrictModel):
    id: str = Field(min_length=1, max_length=80)
    reactor_type: Literal["batch", "semibatch", "cstr", "pfr"]
    volume: QuantityInput
    temperature: QuantityInput
    pressure: QuantityInput
    duration: QuantityInput | None = None
    residence_time: QuantityInput | None = None
    phases: list[PhaseDefinition] = Field(min_length=1)
    initial_inventory: list[InitialAmount] = Field(default_factory=list)
    feeds: list[FeedStream] = Field(default_factory=list)
    heat_transfer: HeatTransfer = Field(
        default_factory=lambda: HeatTransfer(
            mode="isothermal", setpoint_temperature=QuantityInput(value=298.15, unit="K")
        )
    )
    mixing: MixingModel = Field(default_factory=MixingModel)
    gas_transfer: GasTransferModel = Field(default_factory=GasTransferModel)
    catalyst: CatalystDefinition = Field(default_factory=CatalystDefinition)

    @model_validator(mode="after")
    def validate_reactor_requirements(self) -> "ReactorOperation":
        if self.reactor_type in {"batch", "semibatch"} and self.duration is None:
            raise ValueError(f"{self.reactor_type} operation requires duration")
        if self.reactor_type in {"cstr", "pfr"} and self.residence_time is None:
            raise ValueError(f"{self.reactor_type} operation requires residence_time")
        if self.reactor_type in {"cstr", "pfr"} and not self.feeds:
            raise ValueError(f"{self.reactor_type} operation requires at least one feed")
        if self.reactor_type == "batch" and self.feeds:
            raise ValueError("Batch operation cannot contain feeds; use semibatch")
        if self.reactor_type == "batch" and not self.initial_inventory:
            raise ValueError("Batch operation requires an initial inventory")
        phase_ids = [phase.id for phase in self.phases]
        if len(phase_ids) != len(set(phase_ids)):
            raise ValueError("Phase ids must be unique")
        valid_phases = set(phase_ids)
        if any(item.phase_id not in valid_phases for item in self.initial_inventory):
            raise ValueError("Initial inventory references an unknown phase")
        if any(feed.phase_id not in valid_phases for feed in self.feeds):
            raise ValueError("Feed references an unknown phase")
        for phase in self.phases:
            if len(phase.species_ids) != len(set(phase.species_ids)):
                raise ValueError(f"Phase {phase.id} contains duplicate species ids")
        return self


class OperationValidationRequest(StrictModel):
    operation: ReactorOperation
    mechanism: MechanismDefinition

    @model_validator(mode="after")
    def validate_cross_references(self) -> "OperationValidationRequest":
        known_species = {item.id for item in self.mechanism.species}
        phase_species: dict[str, set[str]] = {}
        for phase in self.operation.phases:
            unknown = set(phase.species_ids) - known_species
            if unknown:
                raise ValueError(f"Phase {phase.id} references unknown species: {sorted(unknown)}")
            phase_species[phase.id] = set(phase.species_ids)
        for amount in self.operation.initial_inventory:
            if amount.species_id not in phase_species[amount.phase_id]:
                raise ValueError(f"Species {amount.species_id} is not declared in phase {amount.phase_id}")
        for feed in self.operation.feeds:
            for flow in feed.species_flows:
                if flow.species_id not in phase_species[feed.phase_id]:
                    raise ValueError(f"Feed species {flow.species_id} is not declared in phase {feed.phase_id}")
        if self.operation.catalyst.species_id and self.operation.catalyst.species_id not in known_species:
            raise ValueError("Catalyst species_id is not present in the mechanism")
        return self
