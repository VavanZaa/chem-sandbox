from __future__ import annotations

from enum import Enum
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class SubstanceId(str, Enum):
    sodium_bicarbonate = "sodium_bicarbonate"
    acetic_acid = "acetic_acid"
    water = "water"


SUBSTANCE_LIBRARY = {
    SubstanceId.sodium_bicarbonate: {
        "name": "Baking soda",
        "formula": "NaHCO₃",
        "default_amount": 0.05952,
        "unit": "mol",
        "allowed_units": ["mol", "g"],
        "molar_mass_g_mol": 84.0066,
        "accent": "mint",
    },
    SubstanceId.acetic_acid: {
        "name": "Acetic acid",
        "formula": "CH₃COOH",
        "default_amount": 0.05,
        "unit": "mol",
        "allowed_units": ["mol", "g"],
        "molar_mass_g_mol": 60.052,
        "accent": "amber",
    },
    SubstanceId.water: {
        "name": "Water",
        "formula": "H₂O",
        "default_amount": 2.767,
        "unit": "mol",
        "allowed_units": ["mol", "mL", "g"],
        "molar_mass_g_mol": 18.01528,
        "density_g_mL": 0.997,
        "accent": "blue",
    },
}


class SubstanceInput(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    id: SubstanceId
    amount: float = Field(gt=0, le=10_000)
    unit: Literal["g", "mol", "mL"]


class EnvironmentInput(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    temperature_K: float = Field(default=298.15, ge=250, le=450)
    pressure_atm: float = Field(default=1.0, ge=0.1, le=10)
    solvent: Literal["water", "vacuum", "other"] = "water"
    mixing: Literal["gentle", "standard", "vigorous"] = "standard"


class MolecularInput(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    name: str = Field(min_length=1, max_length=80)
    smiles: str = Field(min_length=1, max_length=500)
    amount_mol: float = Field(default=.05, gt=0, le=10000)
    copies: int = Field(default=1, ge=1, le=3)


class ComputeInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    engine: Literal["all_compatible", "mace_mp", "ani2x", "aimnet2", "aimnet2_rxn"] = "mace_mp"
    device: Literal["auto", "cuda", "cpu"] = "auto"


class SearchInput(BaseModel):
    model_config = ConfigDict(extra="forbid")
    effort: Literal["low", "medium", "high"] = "medium"
    runs: int = Field(default=8, ge=1, le=50)
    seed: int = Field(default=42, ge=0, le=2**31 - 1)
    max_candidates_per_state: int = Field(default=8, ge=1, le=24)
    optimization_steps: int = Field(default=10, ge=0, le=60)
    candidate_optimization_steps: int = Field(default=6, ge=1, le=60)
    max_depth: int = Field(default=3, ge=1, le=12)
    max_evaluations_per_run: int = Field(default=24, ge=1, le=160)
    max_states_per_run: int = Field(default=16, ge=1, le=80)
    frontier_width: int = Field(default=6, ge=1, le=24)


class KnowledgePolicy(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)
    mode: Literal["blind_physics_only"] = "blind_physics_only"
    reference_access: Literal["post_run_evaluation_only"] = "post_run_evaluation_only"
    reaction_family_priors: Literal["generic_non_targeted_operators_only"] = "generic_non_targeted_operators_only"


class ExperimentDefinition(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    experiment_id: str = Field(default="bicarbonate_volcano_v02", min_length=3, max_length=80)
    mode: Literal["bicarbonate", "general"] = "bicarbonate"
    substances: list[SubstanceInput] = Field(default_factory=list)
    molecules: list[MolecularInput] = Field(default_factory=list, max_length=6)
    compute: ComputeInput = Field(default_factory=ComputeInput)
    environment: EnvironmentInput = Field(default_factory=EnvironmentInput)
    search: SearchInput = Field(default_factory=SearchInput)
    knowledge_policy: KnowledgePolicy = Field(default_factory=KnowledgePolicy)

    @model_validator(mode="after")
    def validate_v01_chemistry(self) -> "ExperimentDefinition":
        if self.mode == "general":
            if not self.molecules:
                raise ValueError("General mode requires at least one SMILES molecule")
            from app.chemistry.general_builder import validate_molecules
            validate_molecules(self.molecules, self.compute.engine)
            return self
        if self.compute.engine not in {"mace_mp", "all_compatible"}:
            raise ValueError("Only MACE-MP supports every element in the sodium bicarbonate demo; choose MACE or All compatible")
        ids = [item.id for item in self.substances]
        if len(ids) != len(set(ids)):
            raise ValueError("Each supported substance may appear only once")
        required = set(SUBSTANCE_LIBRARY)
        if set(ids) != required:
            raise ValueError(
                "Bicarbonate mode requires baking soda, acetic acid, and water; choose general mode for SMILES molecules"
            )
        for item in self.substances:
            allowed_units = SUBSTANCE_LIBRARY[item.id]["allowed_units"]
            if item.unit not in allowed_units:
                raise ValueError(
                    f"{item.id.value} supports {', '.join(str(unit) for unit in allowed_units)}"
                )
        return self


def default_experiment() -> ExperimentDefinition:
    return ExperimentDefinition(
        substances=[
            SubstanceInput(id=key, amount=value["default_amount"], unit=value["unit"])
            for key, value in SUBSTANCE_LIBRARY.items()
        ]
    )
