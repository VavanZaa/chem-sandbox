"""Condition-dependent, mass-balanced batch kinetics with declared evidence.

No probability is inferred from engine votes. Model and barrier-sensitivity
spreads are reported separately from experimental/calibrated uncertainty.
"""
from collections import Counter
import warnings
import threading
from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator
from scipy.integrate import solve_ivp
from scipy.optimize import linprog, OptimizeWarning
from app.chemistry.general_builder import parse_molecule

# HiGHS' native worker pool is unnecessary for these tiny LPs and can conflict
# with request-worker lifecycle on Windows. Keep this calculation single-threaded.
_ceiling_lock = threading.Lock()


class SpeciesAmount(BaseModel):
    model_config = ConfigDict(allow_inf_nan=False, extra="forbid")
    id: str = Field(min_length=1, max_length=40, pattern=r"^[A-Za-z][A-Za-z0-9_]*$")
    name: str = Field(min_length=1, max_length=80)
    smiles: str = Field(min_length=1, max_length=500)
    initial_mol: float = Field(ge=0, le=10000)


class KineticChannel(BaseModel):
    model_config = ConfigDict(allow_inf_nan=False, extra="forbid")
    id: str = Field(min_length=1, max_length=80)
    reactants: dict[str, int]
    products: dict[str, int]
    activation_free_energy_kj_mol: float = Field(ge=0, le=500)
    reference_temperature_K: float = Field(default=298.15, ge=200, le=1500)
    activation_entropy_j_mol_K: float = Field(default=0, ge=-1000, le=1000)
    source: str = Field(min_length=3, max_length=1000)

    @model_validator(mode="after")
    def coefficients(self):
        if not self.reactants or not self.products or not self.source.strip():
            raise ValueError("A channel requires reactants, products and a source")
        if any(not 1 <= value <= 10 for value in [*self.reactants.values(), *self.products.values()]):
            raise ValueError("Stoichiometric coefficients must be integers from 1 to 10")
        return self


class KineticModel(BaseModel):
    id: str = Field(min_length=1, max_length=80)
    basis: Literal["hypothetical", "external_calculation", "experimental"] = "hypothetical"
    channels: list[KineticChannel] = Field(min_length=1, max_length=32)


class EfficiencyRequest(BaseModel):
    model_config = ConfigDict(allow_inf_nan=False, extra="forbid")
    species: list[SpeciesAmount] = Field(min_length=2, max_length=24)
    models: list[KineticModel] = Field(min_length=1, max_length=4)
    reference_reactant: str
    target_product: str
    product_pool: list[str] = Field(min_length=1, max_length=24)
    temperature_K: float = Field(default=298.15, ge=200, le=1500)
    duration_s: float = Field(default=60, gt=0, le=1e7)
    phase: Literal["solution", "ideal_gas"] = "solution"
    volume_L: float = Field(default=.05, gt=1e-8, le=10000)
    pressure_atm: float = Field(default=1, gt=.001, le=1000)
    solvent: str = Field(default="water", min_length=1, max_length=80)
    sensitivity_kj_mol: float = Field(default=5, ge=0, le=30)
    acknowledge_assumptions: Literal[True]


def simulate_efficiency(request: EfficiencyRequest) -> dict:
    ids = [item.id for item in request.species]
    if len(set(ids)) != len(ids) or len({model.id for model in request.models}) != len(request.models):
        raise ValueError("Species IDs and model IDs must be unique")
    if request.reference_reactant not in ids or request.target_product not in ids:
        raise ValueError("Select a known reference reactant and target product")
    if request.reference_reactant == request.target_product:
        raise ValueError("Reference reactant and target product must differ")
    if len(set(request.product_pool)) != len(request.product_pool) or not set(request.product_pool) <= set(ids) or request.target_product not in request.product_pool:
        raise ValueError("Product pool must contain unique known IDs including the target")
    index = {name: i for i, name in enumerate(ids)}
    initial = np.array([item.initial_mol for item in request.species])
    ref, target = index[request.reference_reactant], index[request.target_product]
    if initial[ref] <= 0:
        raise ValueError("Reference reactant must have a positive initial amount")
    inventories = []
    for item in request.species:
        molecule = parse_molecule(item.smiles)
        inventory = Counter(atom.GetSymbol() for atom in molecule.GetAtoms())
        inventory["charge"] = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
        inventories.append(inventory)
    elements = sorted({key for inventory in inventories for key in inventory})
    atomic_matrix = np.array([[inventory.get(key, 0) for inventory in inventories] for key in elements])
    results = []
    topology = None
    for model in request.models:
        if len({channel.id for channel in model.channels}) != len(model.channels):
            raise ValueError("Channel IDs must be unique in each model")
        signature = sorted((tuple(sorted(channel.reactants.items())), tuple(sorted(channel.products.items()))) for channel in model.channels)
        if topology is not None and signature != topology:
            raise ValueError("Model comparison requires identical reaction channels; missing channels cannot be treated as zero probability")
        topology = signature
        matrix = np.zeros((len(ids), len(model.channels)))
        orders = np.zeros_like(matrix)
        for j, channel in enumerate(model.channels):
            if not (set(channel.reactants) | set(channel.products)) <= set(ids):
                raise ValueError(f"Unknown species in {channel.id}")
            for species, coefficient in channel.reactants.items():
                matrix[index[species], j] -= coefficient
                orders[index[species], j] = coefficient
            for species, coefficient in channel.products.items():
                matrix[index[species], j] += coefficient
            if np.any(atomic_matrix @ matrix[:, j]):
                raise ValueError(f"Channel {channel.id} does not conserve atoms and charge")
            if not np.any(matrix[:, j]):
                raise ValueError(f"Channel {channel.id} changes no species")
        with _ceiling_lock, warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Unrecognized options detected.*threads.*", category=OptimizeWarning)
            ceiling = linprog(-matrix[target], A_ub=-matrix, b_ub=initial, bounds=(0, None), method="highs", options={"threads": 1})
        if not ceiling.success:
            raise ValueError("Could not establish a finite stoichiometric product ceiling")
        theoretical = max(0.0, -float(ceiling.fun))
        base_barriers = np.array([channel.activation_free_energy_kj_mol -
                                 (request.temperature_K - channel.reference_temperature_K) * channel.activation_entropy_j_mol_K / 1000
                                 for channel in model.channels])
        if np.any(base_barriers < 0):
            raise ValueError("Temperature extrapolation gave a negative activation free energy; supply evidence valid at these conditions")
        # Unit-activity transition-state convention: 1 mol/L solution or 1 atm gas.
        standard_M = 1.0 if request.phase == "solution" else 1 / (.082057366080960 * request.temperature_K)
        prefactor = 1.380649e-23 * request.temperature_K / 6.62607015e-34
        scenarios = [("nominal", np.zeros(len(base_barriers)))]
        if request.sensitivity_kj_mol:
            direction = np.where(matrix[target] > 0, -1.0, 1.0)
            scenarios += [("all_barriers_higher", np.ones(len(base_barriers))), ("all_barriers_lower", -np.ones(len(base_barriers))),
                          ("target_favored", direction), ("target_disfavored", -direction)]
        scenario_results = []
        for label, shifts in scenarios:
            barriers = np.maximum(0, base_barriers + shifts * request.sensitivity_kj_mol)
            rates = prefactor * np.exp(-barriers * 1000 / (8.31446261815324 * request.temperature_K))
            def volume(moles):
                if request.phase == "solution":
                    return request.volume_L
                return max(moles.sum() * .082057366080960 * request.temperature_K / request.pressure_atm, 1e-15)
            def rhs(_time, moles):
                positive = np.maximum(moles, 0)
                vessel = volume(positive)
                activity = positive / vessel / standard_M
                flux = rates * standard_M * np.prod(activity[:, None] ** orders, axis=0)
                return vessel * matrix @ flux
            solution = solve_ivp(rhs, [0, request.duration_s], initial, method="BDF",
                                 t_eval=np.linspace(0, request.duration_s, 51), rtol=1e-7,
                                 atol=max(initial.max() * 1e-10, 1e-14))
            if not solution.success or not np.isfinite(solution.y).all():
                raise ValueError(f"Kinetic integration failed for {model.id}/{label}")
            tolerance = max(initial.max() * 1e-6, 1e-12)
            if solution.y.min() < -tolerance or np.max(np.abs(atomic_matrix @ (solution.y - initial[:, None]))) > tolerance * 20:
                raise ValueError("Kinetic integration failed atom/charge conservation")
            final = np.maximum(solution.y[:, -1], 0)
            formed = max(0.0, final[target] - initial[target])
            product_total = sum(max(0.0, final[index[name]] - initial[index[name]]) for name in request.product_pool)
            metrics = {"scenario": label, "conversion": float((initial[ref] - final[ref]) / initial[ref]),
                       "target_yield": float(formed / theoretical) if theoretical > 1e-14 else None,
                       "target_formed_mol": float(formed),
                       "target_product_share": float(formed / product_total) if product_total > 1e-14 else None,
                       "final_moles": dict(zip(ids, final.tolist()))}
            scenario_results.append(metrics)
            if label == "nominal":
                samples = [{"time_s": float(t), "moles": dict(zip(ids, np.maximum(y, 0).tolist()))}
                           for t, y in zip(solution.t, solution.y.T)]
        nominal = scenario_results[0]
        yields = [item["target_yield"] for item in scenario_results if item["target_yield"] is not None]
        results.append({"model_id": model.id, "basis": model.basis, "nominal": nominal, "scenarios": scenario_results,
                        "yield_sensitivity_range": [min(yields), max(yields)] if yields else None,
                        "stoichiometric_ceiling_mol": theoretical, "samples": samples,
                        "barriers_at_temperature_kj_mol": dict(zip([channel.id for channel in model.channels], base_barriers.tolist()))})
    return {"status": "conditional_prediction", "conditions": request.model_dump(exclude={"models", "species"}),
            "models": results, "evidence": request.model_dump(), "calibrated_product_probabilities": None,
            "warnings": ["Yield is net target formation divided by a stoichiometric network upper bound, not an imposed 100% outcome.",
                         "The listed-product share is a molar composition of the selected product pool, not a probability that an engine is correct.",
                         "The model/sensitivity ranges are scenario spreads, not confidence intervals. Agreement can share systematic bias.",
                         "Each directed channel is treated as elementary mass action. Reverse and competing channels must be explicitly included with appropriate standard-state barriers.",
                         "Solution: fixed volume and 1 M standard state; gas: ideal, isothermal constant-pressure vessel and 1 atm standard state.",
                         "Solvent/pH/catalysts/mixing effects must already be encoded in supplied barriers; merely naming a solvent does not calculate its effect.",
                         "Temperature extrapolation assumes constant activation entropy; zero entropy assumes constant activation free energy. Validate this approximation.",
                         "Reverse rates and thermodynamic cycle consistency are user-supplied, not inferred from reaction energies. Sensitivity shifts are clipped at zero barrier.",
                         "No work-up, gas collection or isolation loss is modeled. Near-complete conversion is possible and is not artificially capped below 100%."],
            "references": ["https://goldbook.iupac.org/terms/view/T06470"]}
