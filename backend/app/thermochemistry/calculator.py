from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.special import logsumexp

from .schema import (
    ActivationFreeEnergyRequest,
    EvidenceQuantity,
    MolecularIdentity,
    ReactionFreeEnergyRequest,
    ThermochemicalConformer,
    ThermochemicalEnsembleRequest,
)


R_KJ_MOL_K = 0.00831446261815324
IMAGINARY_THRESHOLD_CM_INV = -20.0
MAX_CERTIFIED_FORCE_EV_A = 0.05


def _evidence_value(value: EvidenceQuantity | None, target: str, label: str, problems: list[str]) -> float | None:
    if value is None:
        problems.append(f"missing {label}")
        return None
    if value.basis == "hypothetical":
        problems.append(f"{label} is hypothetical")
    try:
        return value.as_unit(target)
    except ValueError as error:
        problems.append(str(error))
        return None


def _identity_key(identity: MolecularIdentity) -> tuple[Any, ...]:
    return (
        tuple(sorted(identity.elements.items())),
        identity.charge,
        identity.multiplicity,
        identity.phase,
        identity.standard_state,
        identity.electronic_method,
        identity.solvation_model,
    )


def certify_conformer(conformer: ThermochemicalConformer, temperature_K: float) -> dict[str, Any]:
    problems: list[str] = []
    point = conformer.stationary_point
    if conformer.electronic_energy.method != conformer.identity.electronic_method:
        problems.append("electronic-energy method does not match molecular identity method")
    if not point.geometry_optimized:
        problems.append("geometry is not marked optimized")
    if not point.optimizer_converged:
        problems.append("geometry optimizer did not converge")
    max_force = _evidence_value(point.max_force, "eV/angstrom", "maximum force", problems)
    if max_force is not None and max_force > MAX_CERTIFIED_FORCE_EV_A:
        problems.append(f"maximum force {max_force:.6g} eV/angstrom exceeds {MAX_CERTIFIED_FORCE_EV_A}")
    expected_modes = 3 * point.atom_count - (5 if point.linear else 6)
    frequencies = point.frequencies.values_cm_inv
    if len(frequencies) != expected_modes:
        problems.append(f"expected {expected_modes} vibrational modes, received {len(frequencies)}")
    if point.frequencies.basis == "hypothetical":
        problems.append("frequencies are hypothetical")
    imaginary = [value for value in frequencies if value < IMAGINARY_THRESHOLD_CM_INV]
    if point.kind == "minimum" and imaginary:
        problems.append(f"minimum has {len(imaginary)} significant imaginary mode(s)")
    if point.kind == "transition_state":
        if len(imaginary) != 1:
            problems.append(f"transition state requires exactly one significant imaginary mode; found {len(imaginary)}")
        if not point.reaction_coordinate_verified or not point.coordinate_verification_source:
            problems.append("transition-state reaction coordinate is not verified")

    electronic = _evidence_value(conformer.electronic_energy, "kJ/mol", "electronic energy", problems)
    zpe = _evidence_value(conformer.zero_point_energy, "kJ/mol", "zero-point energy", problems)
    thermal_h = _evidence_value(
        conformer.thermal_enthalpy_correction, "kJ/mol", "thermal enthalpy correction", problems
    )
    entropy = _evidence_value(conformer.entropy, "J/mol/K", "entropy", problems)
    standard = _evidence_value(
        conformer.standard_state_correction, "kJ/mol", "standard-state correction", problems
    )
    if conformer.identity.phase == "solution":
        solvation = _evidence_value(
            conformer.solvation_free_energy, "kJ/mol", "solvation free energy", problems
        )
    elif conformer.solvation_free_energy is None:
        solvation = 0.0
    else:
        solvation = _evidence_value(
            conformer.solvation_free_energy, "kJ/mol", "solvation free energy", problems
        )
    values = (electronic, zpe, thermal_h, entropy, solvation, standard)
    gibbs = None
    if all(value is not None for value in values):
        gibbs = float(electronic + zpe + thermal_h - temperature_K * entropy / 1000.0 + solvation + standard)
    return {
        "conformer_id": conformer.id,
        "certified": not problems,
        "stationary_point_kind": point.kind,
        "expected_vibrational_modes": expected_modes,
        "significant_imaginary_modes_cm_inv": imaginary,
        "gibbs_free_energy_kj_mol": gibbs if not problems else None,
        "uncertified_calculated_value_kj_mol": gibbs if problems else None,
        "problems": problems,
    }


def evaluate_ensemble(request: ThermochemicalEnsembleRequest) -> dict[str, Any]:
    temperature_K = request.temperature.as_unit("K")
    if temperature_K <= 0:
        raise ValueError("Temperature must be greater than 0 K")
    identity_keys = {_identity_key(item.identity) for item in request.conformers}
    shared_problems = []
    if len(identity_keys) != 1:
        shared_problems.append("conformers do not share composition, charge, spin, phase, standard state, and method")
    kinds = {item.stationary_point.kind for item in request.conformers}
    if len(kinds) != 1:
        shared_problems.append("an ensemble cannot mix minima and transition states")
    records = [certify_conformer(item, temperature_K) for item in request.conformers]
    all_conformers_certified = all(item["certified"] for item in records) and not shared_problems
    submitted_gibbs = None
    populations: list[dict[str, Any]] = []
    if all_conformers_certified:
        energies = np.asarray([item["gibbs_free_energy_kj_mol"] for item in records], dtype=float)
        degeneracies = np.asarray([item.degeneracy for item in request.conformers], dtype=float)
        log_weights = np.log(degeneracies) - energies / (R_KJ_MOL_K * temperature_K)
        log_partition = float(logsumexp(log_weights))
        submitted_gibbs = float(-R_KJ_MOL_K * temperature_K * log_partition)
        probabilities = np.exp(log_weights - log_partition)
        populations = [
            {"conformer_id": conformer.id, "population": float(probability)}
            for conformer, probability in zip(request.conformers, probabilities, strict=True)
        ]
    certified = all_conformers_certified and request.conformer_search_complete
    certification_problems = list(shared_problems)
    if not request.conformer_search_complete:
        certification_problems.append("conformer search is not declared complete")
    return {
        "ensemble_id": request.id,
        "temperature_K": temperature_K,
        "certified": certified,
        "stationary_point_kind": next(iter(kinds)) if len(kinds) == 1 else None,
        "ensemble_gibbs_free_energy_kj_mol": submitted_gibbs if certified else None,
        "conditional_submitted_ensemble_gibbs_kj_mol": submitted_gibbs,
        "conformer_populations": populations,
        "conformers": records,
        "certification_problems": certification_problems,
        "conformer_search_protocol": request.conformer_search_protocol,
        "claim": (
            "Certified Gibbs free energy for the declared conformer ensemble."
            if certified
            else "No certified ensemble Gibbs free energy is claimed."
        ),
    }


def _single_identity(request: ThermochemicalEnsembleRequest) -> MolecularIdentity:
    first = request.conformers[0].identity
    if any(_identity_key(item.identity) != _identity_key(first) for item in request.conformers[1:]):
        raise ValueError(f"Ensemble {request.id} does not have one consistent molecular identity")
    return first


def evaluate_activation_free_energy(request: ActivationFreeEnergyRequest) -> dict[str, Any]:
    reactant = evaluate_ensemble(request.reactant)
    transition_state = evaluate_ensemble(request.transition_state)
    problems: list[str] = []
    if abs(reactant["temperature_K"] - transition_state["temperature_K"]) > 1e-9:
        problems.append("reactant and transition state temperatures differ")
    reactant_identity = _single_identity(request.reactant)
    ts_identity = _single_identity(request.transition_state)
    if _identity_key(reactant_identity) != _identity_key(ts_identity):
        problems.append("reactant and transition state identity/method/standard state differ")
    if reactant["stationary_point_kind"] != "minimum":
        problems.append("reactant ensemble is not composed of minima")
    if transition_state["stationary_point_kind"] != "transition_state":
        problems.append("transition-state ensemble is not composed of transition states")
    if not reactant["certified"]:
        problems.append("reactant ensemble is not certified")
    if not transition_state["certified"]:
        problems.append("transition-state ensemble is not certified")
    barrier = None
    if not problems:
        barrier = transition_state["ensemble_gibbs_free_energy_kj_mol"] - reactant["ensemble_gibbs_free_energy_kj_mol"]
        if barrier < 0:
            problems.append("certified inputs produce a negative activation free energy; inspect reference states")
            barrier = None
    return {
        "certified": not problems,
        "temperature_K": reactant["temperature_K"],
        "activation_gibbs_free_energy_kj_mol": barrier,
        "reactant": reactant,
        "transition_state": transition_state,
        "problems": problems,
        "claim": "Certified activation Gibbs free energy." if not problems else "No ΔG‡ claim is made.",
    }


def evaluate_reaction_free_energy(request: ReactionFreeEnergyRequest) -> dict[str, Any]:
    reactants = [(term, evaluate_ensemble(term.ensemble)) for term in request.reactants]
    products = [(term, evaluate_ensemble(term.ensemble)) for term in request.products]
    all_items = reactants + products
    problems: list[str] = []
    temperatures = {round(result["temperature_K"], 10) for _, result in all_items}
    if len(temperatures) != 1:
        problems.append("all thermochemical ensembles must use the same temperature")
    if any(not result["certified"] for _, result in all_items):
        problems.append("every reactant and product ensemble must be certified")

    element_delta: dict[str, float] = {}
    charge_delta = 0.0
    methods: set[str] = set()
    for sign, terms in ((-1.0, reactants), (1.0, products)):
        for term, _ in terms:
            identity = _single_identity(term.ensemble)
            methods.add(identity.electronic_method)
            for element, count in identity.elements.items():
                element_delta[element] = element_delta.get(element, 0.0) + sign * term.coefficient * count
            charge_delta += sign * term.coefficient * identity.charge
    residual = {key: value for key, value in element_delta.items() if abs(value) > 1e-10}
    if residual or abs(charge_delta) > 1e-10:
        problems.append(f"reaction is not atom/charge balanced: elements={residual}, charge={charge_delta}")
    if len(methods) != 1:
        problems.append("reactant and product electronic methods differ")
    delta = None
    if not problems:
        product_g = sum(term.coefficient * result["ensemble_gibbs_free_energy_kj_mol"] for term, result in products)
        reactant_g = sum(term.coefficient * result["ensemble_gibbs_free_energy_kj_mol"] for term, result in reactants)
        delta = float(product_g - reactant_g)
    return {
        "certified": not problems,
        "temperature_K": next(iter(temperatures)) if len(temperatures) == 1 else None,
        "reaction_gibbs_free_energy_kj_mol": delta,
        "reactants": [result for _, result in reactants],
        "products": [result for _, result in products],
        "balance": {"element_residual": residual, "charge_residual": charge_delta},
        "problems": problems,
        "claim": "Certified reaction Gibbs free energy." if not problems else "No ΔG reaction claim is made.",
    }
