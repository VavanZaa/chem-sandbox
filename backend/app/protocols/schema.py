from __future__ import annotations

import re
from typing import Literal

from pydantic import Field, model_validator

from app.operations.schema import QuantityInput, StrictModel


SPECIES_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.+-]{0,63}$")


class ProtocolSpecies(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    name: str = Field(min_length=1, max_length=120)
    smiles: str | None = Field(default=None, min_length=1, max_length=2000)
    composition: dict[str, int] | None = Field(default=None, min_length=1)
    charge: int | None = Field(default=None, ge=-100, le=100)
    identity_note: str | None = Field(default=None, max_length=500)

    @model_validator(mode="after")
    def identity(self) -> "ProtocolSpecies":
        if (self.smiles is None) == (self.composition is None):
            raise ValueError("Each protocol species requires exactly one of smiles or composition")
        if self.composition is not None and any(
            not re.fullmatch(r"[A-Z][a-z]?", element) or count <= 0
            for element, count in self.composition.items()
        ):
            raise ValueError("Composition keys must be element symbols with positive counts")
        if self.composition is not None and self.charge is None:
            raise ValueError("Composition-only species require an explicit charge")
        return self


class ProtocolRateLaw(StrictModel):
    kind: Literal["mass_action", "eyring_mass_action", "michaelis_menten"]
    source: str = Field(min_length=3, max_length=1000)
    basis: Literal["hypothetical", "external_calculation", "experimental", "literature"]
    rate_constant: QuantityInput | None = None
    activation_free_energy: QuantityInput | None = None
    reference_temperature: QuantityInput | None = None
    activation_entropy: QuantityInput | None = None
    substrate_id: str | None = Field(default=None, pattern=SPECIES_ID.pattern)
    catalyst_id: str | None = Field(default=None, pattern=SPECIES_ID.pattern)
    kcat: QuantityInput | None = None
    km: QuantityInput | None = None

    @model_validator(mode="after")
    def required_parameters(self) -> "ProtocolRateLaw":
        if self.kind == "mass_action":
            if self.rate_constant is None:
                raise ValueError("mass_action requires a rate_constant")
            forbidden = [self.activation_free_energy, self.reference_temperature, self.activation_entropy,
                         self.substrate_id, self.catalyst_id, self.kcat, self.km]
        elif self.kind == "eyring_mass_action":
            if self.activation_free_energy is None:
                raise ValueError("eyring_mass_action requires activation_free_energy")
            forbidden = [self.rate_constant, self.substrate_id, self.catalyst_id, self.kcat, self.km]
        else:
            if not self.substrate_id or not self.catalyst_id or self.kcat is None or self.km is None:
                raise ValueError("michaelis_menten requires substrate_id, catalyst_id, kcat and km")
            forbidden = [self.rate_constant, self.activation_free_energy, self.reference_temperature,
                         self.activation_entropy]
        if any(value is not None for value in forbidden):
            raise ValueError(f"Parameters incompatible with {self.kind} were supplied")
        return self


class ProtocolReaction(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    name: str = Field(min_length=1, max_length=160)
    stoichiometry: dict[str, int] = Field(min_length=2, max_length=50)
    rate_law: ProtocolRateLaw

    @model_validator(mode="after")
    def coefficients(self) -> "ProtocolReaction":
        if any(value == 0 or abs(value) > 20 for value in self.stoichiometry.values()):
            raise ValueError("Reaction coefficients must be nonzero integers with magnitude at most 20")
        if not any(value < 0 for value in self.stoichiometry.values()) or not any(
            value > 0 for value in self.stoichiometry.values()
        ):
            raise ValueError("A protocol reaction requires reactants and products")
        return self


class ProtocolAddition(StrictModel):
    species_id: str = Field(pattern=SPECIES_ID.pattern)
    amount: QuantityInput


class ProtocolStage(StrictModel):
    id: str = Field(pattern=SPECIES_ID.pattern)
    name: str = Field(min_length=1, max_length=160)
    operation: Literal["addition", "reaction", "hold", "sample", "transfer"]
    duration: QuantityInput = Field(default_factory=lambda: QuantityInput(value=0, unit="s"))
    additions: list[ProtocolAddition] = Field(default_factory=list, max_length=100)
    active_reactions: list[str] = Field(default_factory=list, max_length=200)
    temperature: QuantityInput | None = None
    volume: QuantityInput | None = None
    withdraw_fraction: float = Field(default=0, ge=0, lt=1)
    recovery_fraction: float = Field(default=1, gt=0, le=1)
    note: str | None = Field(default=None, max_length=1000)

    @model_validator(mode="after")
    def operation_fields(self) -> "ProtocolStage":
        duration_s = self.duration.as_unit("s")
        if self.operation == "reaction" and (duration_s <= 0 or not self.active_reactions):
            raise ValueError("A reaction stage requires positive duration and active_reactions")
        if self.operation == "hold" and (duration_s <= 0 or self.active_reactions):
            raise ValueError("A hold stage requires positive duration and no active reactions")
        if self.operation == "addition" and (not self.additions or duration_s != 0 or self.active_reactions):
            raise ValueError("An addition stage requires additions, zero duration and no active reactions")
        if self.operation == "sample" and (self.withdraw_fraction <= 0 or duration_s != 0 or self.active_reactions):
            raise ValueError("A sample stage requires a positive withdraw_fraction and zero duration")
        if self.operation == "transfer" and (self.recovery_fraction >= 1 or duration_s != 0 or self.active_reactions):
            raise ValueError("A transfer stage requires recovery_fraction below one and zero duration")
        if self.operation not in {"sample", "transfer"} and (
            self.withdraw_fraction or self.recovery_fraction != 1
        ):
            raise ValueError("Withdrawal/recovery fields are restricted to sample or transfer stages")
        return self


class LabProtocolRequest(StrictModel):
    id: str = Field(min_length=1, max_length=100)
    species: list[ProtocolSpecies] = Field(min_length=2, max_length=200)
    reactions: list[ProtocolReaction] = Field(min_length=1, max_length=500)
    initial_amounts: dict[str, QuantityInput] = Field(min_length=1)
    initial_volume: QuantityInput
    initial_temperature: QuantityInput
    stages: list[ProtocolStage] = Field(min_length=1, max_length=200)
    points_per_reaction_stage: int = Field(default=21, ge=2, le=501)
    include_sbml: bool = False
    cross_check_roadrunner: bool = True
    acknowledge_evidence_limits: Literal[True]

    @model_validator(mode="after")
    def unique_ids(self) -> "LabProtocolRequest":
        species_ids = [item.id for item in self.species]
        reaction_ids = [item.id for item in self.reactions]
        stage_ids = [item.id for item in self.stages]
        if len(species_ids) != len(set(species_ids)):
            raise ValueError("Protocol species ids must be unique")
        if len(reaction_ids) != len(set(reaction_ids)):
            raise ValueError("Protocol reaction ids must be unique")
        if len(stage_ids) != len(set(stage_ids)):
            raise ValueError("Protocol stage ids must be unique")
        known_species, known_reactions = set(species_ids), set(reaction_ids)
        if set(self.initial_amounts) - known_species:
            raise ValueError("Initial amounts reference unknown species")
        for reaction in self.reactions:
            if set(reaction.stoichiometry) - known_species:
                raise ValueError(f"Reaction {reaction.id} references unknown species")
            rate = reaction.rate_law
            if rate.substrate_id and rate.substrate_id not in known_species:
                raise ValueError(f"Reaction {reaction.id} has an unknown Michaelis-Menten substrate")
            if rate.catalyst_id and rate.catalyst_id not in known_species:
                raise ValueError(f"Reaction {reaction.id} has an unknown catalyst")
        for stage in self.stages:
            if set(stage.active_reactions) - known_reactions:
                raise ValueError(f"Stage {stage.id} activates unknown reactions")
            if any(addition.species_id not in known_species for addition in stage.additions):
                raise ValueError(f"Stage {stage.id} adds an unknown species")
        return self


class SBMLInspectionRequest(StrictModel):
    sbml: str = Field(min_length=20, max_length=5_000_000)
