"""Conditional first-order state-network kinetics, not bulk reaction yields."""
import math
from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator
from scipy.linalg import expm


class RateEvidence(BaseModel):
    model_config = ConfigDict(allow_inf_nan=False, extra="forbid")
    edge_id: str
    rate_s_inv: float | None = Field(default=None, ge=0, le=1e14)
    activation_free_energy_ev: float | None = Field(default=None, ge=0, le=20)
    temperature_K: float = Field(ge=250, le=450)
    source: str = Field(min_length=3, max_length=1000)
    basis: Literal["hypothetical", "experimental", "external_calculation"] = "hypothetical"

    @model_validator(mode="after")
    def one_parameter(self):
        if (self.rate_s_inv is None) == (self.activation_free_energy_ev is None):
            raise ValueError("Supply exactly one first-order rate or activation FREE-energy barrier per directed edge")
        if not self.source.strip():
            raise ValueError("A source or explicit hypothetical assumption is required")
        return self


class KineticsRequest(BaseModel):
    model_config = ConfigDict(allow_inf_nan=False, extra="forbid")
    initial_state_id: str
    temperature_K: float = Field(default=298.15, ge=250, le=450)
    duration_s: float = Field(default=10, gt=0, le=1e8)
    points: int = Field(default=51, ge=2, le=201)
    rates: list[RateEvidence] = Field(max_length=1000)
    acknowledge_conditional_model: Literal[True]


def simulate_kinetics(graph: dict, request: KineticsRequest) -> dict:
    ids = [node["id"] for node in graph["nodes"]]
    if request.initial_state_id not in ids:
        raise ValueError("Initial state must be in the searched network")
    if len(ids) > 160:
        raise ValueError("This bounded kinetics solver supports at most 160 states")
    edges = {edge["id"]: edge for edge in graph["edges"]}
    if not edges:
        raise ValueError("No candidate pathways exist; kinetics are unavailable")
    evidence = {item.edge_id: item for item in request.rates}
    if len(evidence) != len(request.rates):
        raise ValueError("Duplicate edge evidence is not allowed")
    if set(evidence) != set(edges):
        missing = sorted(set(edges) - set(evidence))
        extra = sorted(set(evidence) - set(edges))
        raise ValueError(f"Supply evidence for every directed pathway; missing={missing}, unknown={extra}")
    index = {state: i for i, state in enumerate(ids)}
    generator = np.zeros((len(ids), len(ids)))
    used_rates = []
    for edge_id, item in evidence.items():
        if abs(item.temperature_K - request.temperature_K) > 1e-6:
            raise ValueError("All evidence temperatures must match the simulation; no temperature extrapolation is assumed")
        # Eyring with transmission coefficient 1 and fixed first-order state convention.
        rate = item.rate_s_inv
        if rate is None:
            rate = (1.380649e-23 * request.temperature_K / 6.62607015e-34) * math.exp(
                -item.activation_free_energy_ev / (8.617333262145e-5 * request.temperature_K))
        source, target = index[edges[edge_id]["source"]], index[edges[edge_id]["target"]]
        generator[target, source] += rate
        generator[source, source] -= rate
        used_rates.append({**item.model_dump(), "computed_rate_s_inv": rate})
    if np.linalg.norm(generator, 1) * request.duration_s > 1e10:
        raise ValueError("Rate × duration exceeds the numerical safety limit; shorten duration or refine rates")
    initial = np.zeros(len(ids))
    initial[index[request.initial_state_id]] = 1
    samples = []
    for time_s in np.linspace(0, request.duration_s, request.points):
        population = expm(generator * time_s) @ initial
        if not np.isfinite(population).all() or population.min() < -1e-8 or abs(population.sum() - 1) > 1e-6:
            raise ValueError("Kinetic integration failed probability conservation; revise rate/time scales")
        population = np.maximum(population, 0)
        population /= population.sum()
        samples.append({"time_s": float(time_s), "fractions": dict(zip(ids, population.tolist()))})
    return {"status": "conditional", "temperature_K": request.temperature_K,
            "initial_state_id": request.initial_state_id, "rates": used_rates, "samples": samples,
            "actual_yield": None,
            "warnings": ["These are conditional state populations, not measured conversion, selectivity or isolated yield.",
                         "All supplied evidence is user-declared, not independently validated by ChemSandbox.",
                         "Assumes constant effective first-order rates, one initial state and a closed, complete network. Missing pathways can change the answer.",
                         "Reverse pathways require separate evidence; reversibility and detailed balance are not inferred.",
                         "For barriers: k = (kB T / h) exp(−ΔG‡ / kB T), transmission coefficient 1. Raw MACE ΔE is never used as ΔG‡.",
                         "Bulk multi-reactant yields additionally need concentration-dependent mass action, solvent/pH, gas exchange, and calibrated rates."],
            "references": ["https://goldbook.iupac.org/terms/view/T06470", "https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.expm.html"]}
