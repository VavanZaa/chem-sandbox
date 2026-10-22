from __future__ import annotations

import re
from typing import Literal

from pydantic import Field, model_validator

from app.operations.schema import QuantityInput, StrictModel


class EvidenceQuantity(QuantityInput):
    source: str = Field(min_length=3, max_length=500)
    method: str = Field(min_length=2, max_length=200)
    basis: Literal["computed", "experimental", "literature", "hypothetical"]
    uncertainty: QuantityInput | None = None


class FrequencyEvidence(StrictModel):
    values_cm_inv: list[float]
    source: str = Field(min_length=3, max_length=500)
    method: str = Field(min_length=2, max_length=200)
    basis: Literal["computed", "experimental", "literature", "hypothetical"]


class MolecularIdentity(StrictModel):
    formula: str = Field(min_length=1, max_length=120)
    elements: dict[str, int] = Field(min_length=1)
    charge: int = Field(default=0, ge=-20, le=20)
    multiplicity: int = Field(default=1, ge=1, le=20)
    phase: Literal["gas", "solution"]
    standard_state: str = Field(min_length=2, max_length=120)
    electronic_method: str = Field(min_length=2, max_length=200)
    solvation_model: str | None = Field(default=None, max_length=200)

    @model_validator(mode="after")
    def validate_identity(self) -> "MolecularIdentity":
        if any(not re.fullmatch(r"[A-Z][a-z]?", element) or count <= 0 for element, count in self.elements.items()):
            raise ValueError("Identity element keys must be symbols with positive counts")
        if self.phase == "solution" and not self.solvation_model:
            raise ValueError("Solution-phase thermochemistry requires a declared solvation_model")
        return self


class StationaryPointEvidence(StrictModel):
    kind: Literal["minimum", "transition_state"]
    atom_count: int = Field(ge=2, le=10000)
    linear: bool = False
    geometry_optimized: bool
    optimizer_converged: bool
    max_force: EvidenceQuantity
    frequencies: FrequencyEvidence
    reaction_coordinate_verified: bool = False
    coordinate_verification_source: str | None = Field(default=None, max_length=500)


class ThermochemicalConformer(StrictModel):
    id: str = Field(min_length=1, max_length=80)
    degeneracy: int = Field(default=1, ge=1, le=1000)
    identity: MolecularIdentity
    stationary_point: StationaryPointEvidence
    electronic_energy: EvidenceQuantity
    zero_point_energy: EvidenceQuantity | None = None
    thermal_enthalpy_correction: EvidenceQuantity | None = None
    entropy: EvidenceQuantity | None = None
    solvation_free_energy: EvidenceQuantity | None = None
    standard_state_correction: EvidenceQuantity | None = None


class ThermochemicalEnsembleRequest(StrictModel):
    id: str = Field(min_length=1, max_length=80)
    temperature: QuantityInput
    conformers: list[ThermochemicalConformer] = Field(min_length=1)
    conformer_search_complete: bool
    conformer_search_protocol: str = Field(min_length=3, max_length=1000)

    @model_validator(mode="after")
    def unique_conformers(self) -> "ThermochemicalEnsembleRequest":
        ids = [item.id for item in self.conformers]
        if len(ids) != len(set(ids)):
            raise ValueError("Conformer ids must be unique")
        return self


class ActivationFreeEnergyRequest(StrictModel):
    reactant: ThermochemicalEnsembleRequest
    transition_state: ThermochemicalEnsembleRequest


class StoichiometricEnsemble(StrictModel):
    coefficient: float = Field(gt=0, le=1000)
    ensemble: ThermochemicalEnsembleRequest


class ReactionFreeEnergyRequest(StrictModel):
    reactants: list[StoichiometricEnsemble] = Field(min_length=1)
    products: list[StoichiometricEnsemble] = Field(min_length=1)


class HarmonicCorrectionRequest(StrictModel):
    temperature: QuantityInput
    frequencies: FrequencyEvidence
    stationary_point_kind: Literal["minimum", "transition_state"]
    low_frequency_treatment: Literal["strict_harmonic", "raise_to_cutoff"] = "raise_to_cutoff"
    low_frequency_cutoff_cm_inv: float = Field(default=50.0, gt=0, le=500)
