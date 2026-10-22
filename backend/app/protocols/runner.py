"""Piecewise laboratory-protocol simulator with explicit material accounting."""
from __future__ import annotations

from collections import Counter
import math
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from scipy.integrate import solve_ivp

from app.chemistry.general_builder import parse_molecule

from .schema import LabProtocolRequest, ProtocolRateLaw, ProtocolReaction


R_GAS_J_MOL_K = 8.31446261815324
K_B_J_K = 1.380649e-23
H_J_S = 6.62607015e-34


def _species_identity(item) -> tuple[Counter[str], dict[str, Any]]:
    if item.smiles is None:
        inventory = Counter(item.composition)
        inventory["charge"] = int(item.charge or 0)
        return inventory, {
            "id": item.id,
            "name": item.name,
            "formula": "".join(f"{key}{value if value != 1 else ''}" for key, value in item.composition.items()),
            "isomeric_smiles": None,
            "stereo_status": "declared_composition_only",
            "stereochemistry": None,
            "identity_note": item.identity_note,
        }
    molecule = parse_molecule(item.smiles)
    formal_charge = Chem.GetFormalCharge(molecule)
    if item.charge is not None and item.charge != formal_charge:
        raise ValueError(f"Species {item.id} charge disagrees with its SMILES formal charge")
    inventory = Counter(atom.GetSymbol() for atom in molecule.GetAtoms())
    inventory["charge"] = formal_charge
    parent = Chem.RemoveHs(Chem.Mol(molecule))
    Chem.AssignStereochemistry(parent, cleanIt=True, force=True)
    potential = list(Chem.FindPotentialStereo(parent))
    centers = Chem.FindMolChiralCenters(
        parent, force=True, includeUnassigned=True, useLegacyImplementation=False
    )
    unspecified = any(str(value.specified) == "Unspecified" for value in potential)
    status = "achiral" if not potential else "unassigned" if unspecified else "specified"
    return inventory, {
        "id": item.id,
        "name": item.name,
        "formula": rdMolDescriptors.CalcMolFormula(parent),
        "isomeric_smiles": Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True),
        "stereo_status": status,
        "stereochemistry": [f"atom {index + 1}: {label}" for index, label in centers],
        "identity_note": item.identity_note,
    }


def _prepare(request: LabProtocolRequest) -> dict[str, Any]:
    ids = [item.id for item in request.species]
    index = {key: value for value, key in enumerate(ids)}
    identities: list[dict[str, Any]] = []
    inventories: list[Counter[str]] = []
    for item in request.species:
        inventory, identity = _species_identity(item)
        inventories.append(inventory)
        identities.append(identity)
    conserved = sorted({element for inventory in inventories for element in inventory})
    atomic = np.asarray(
        [[inventory.get(element, 0) for inventory in inventories] for element in conserved], dtype=float
    )
    reaction_index = {reaction.id: value for value, reaction in enumerate(request.reactions)}
    stoichiometry = np.zeros((len(ids), len(request.reactions)), dtype=float)
    reaction_metadata: list[dict[str, Any]] = []
    for column, reaction in enumerate(request.reactions):
        for species_id, coefficient in reaction.stoichiometry.items():
            stoichiometry[index[species_id], column] = coefficient
        residual = atomic @ stoichiometry[:, column]
        if np.max(np.abs(residual), initial=0) > 1e-12:
            detail = {element: float(value) for element, value in zip(conserved, residual) if abs(value) > 1e-12}
            raise ValueError(f"Reaction {reaction.id} does not conserve atoms/charge: {detail}")
        reactant_orders = {
            index[key]: -coefficient for key, coefficient in reaction.stoichiometry.items() if coefficient < 0
        }
        molecularity = sum(reactant_orders.values())
        rate = reaction.rate_law
        if rate.kind == "mass_action":
            target_unit = "1 / second" if molecularity == 1 else (
                f"liter ** {molecularity - 1} / mole ** {molecularity - 1} / second"
            )
            constant = rate.rate_constant.as_unit(target_unit)
            if constant <= 0:
                raise ValueError(f"Reaction {reaction.id} rate constant must be positive")
            parameters = {"rate_constant": constant, "rate_constant_unit": target_unit}
        elif rate.kind == "eyring_mass_action":
            parameters = {
                "activation_free_energy_kj_mol": rate.activation_free_energy.as_unit("kJ/mol"),
                "reference_temperature_K": (
                    rate.reference_temperature.as_unit("K") if rate.reference_temperature else 298.15
                ),
                "activation_entropy_j_mol_K": (
                    rate.activation_entropy.as_unit("J/mol/K") if rate.activation_entropy else 0.0
                ),
            }
        else:
            substrate = index[rate.substrate_id]
            catalyst = index[rate.catalyst_id]
            if substrate not in reactant_orders:
                raise ValueError(f"Michaelis-Menten substrate must be a reactant in {reaction.id}")
            if abs(stoichiometry[catalyst, column]) > 1e-15:
                raise ValueError(f"Michaelis-Menten catalyst must not be consumed in {reaction.id}")
            parameters = {
                "substrate_index": substrate,
                "catalyst_index": catalyst,
                "kcat_s_inv": rate.kcat.as_unit("1/s"),
                "km_mol_L": rate.km.as_unit("mol/L"),
            }
            if parameters["kcat_s_inv"] <= 0 or parameters["km_mol_L"] <= 0:
                raise ValueError(f"Reaction {reaction.id} kcat and Km must be positive")
        reaction_metadata.append({
            "id": reaction.id,
            "name": reaction.name,
            "rate_kind": rate.kind,
            "basis": rate.basis,
            "source": rate.source,
            "reactant_orders": reactant_orders,
            "molecularity": molecularity,
            **parameters,
        })
    return {
        "ids": ids,
        "index": index,
        "identities": identities,
        "conserved": conserved,
        "atomic": atomic,
        "reaction_index": reaction_index,
        "stoichiometry": stoichiometry,
        "reaction_metadata": reaction_metadata,
    }


def _reaction_flux(
    metadata: dict[str, Any], amounts: np.ndarray, volume_L: float, temperature_K: float
) -> float:
    concentrations = np.maximum(amounts, 0) / volume_L
    kind = metadata["rate_kind"]
    if kind == "mass_action":
        activity = math.prod(
            concentrations[index] ** order for index, order in metadata["reactant_orders"].items()
        )
        return metadata["rate_constant"] * activity * volume_L
    if kind == "eyring_mass_action":
        delta_g = metadata["activation_free_energy_kj_mol"] - (
            temperature_K - metadata["reference_temperature_K"]
        ) * metadata["activation_entropy_j_mol_K"] / 1000
        if delta_g < 0:
            raise ValueError(f"Temperature extrapolation gives a negative barrier for {metadata['id']}")
        rate = K_B_J_K * temperature_K / H_J_S * math.exp(
            -delta_g * 1000 / (R_GAS_J_MOL_K * temperature_K)
        )
        standard = 1.0  # mol/L
        activity = math.prod(
            (concentrations[index] / standard) ** order
            for index, order in metadata["reactant_orders"].items()
        )
        return rate * standard * activity * volume_L
    substrate = concentrations[metadata["substrate_index"]]
    return (
        metadata["kcat_s_inv"]
        * max(amounts[metadata["catalyst_index"]], 0)
        * substrate
        / (metadata["km_mol_L"] + substrate)
    )


def _clean_amounts(values: np.ndarray, scale: float) -> np.ndarray:
    tolerance = max(scale * 1e-8, 1e-13)
    if np.min(values, initial=0) < -tolerance:
        raise ValueError("Protocol integration produced a materially negative species amount")
    return np.maximum(values, 0)


def simulate_protocol(request: LabProtocolRequest) -> dict[str, Any]:
    prepared = _prepare(request)
    ids, index = prepared["ids"], prepared["index"]
    initial = np.zeros(len(ids), dtype=float)
    for species_id, amount in request.initial_amounts.items():
        value = amount.as_unit("mol")
        if value < 0:
            raise ValueError("Protocol initial amounts cannot be negative")
        initial[index[species_id]] = value
    if initial.sum() <= 0:
        raise ValueError("Protocol initial inventory must contain a positive amount")
    volume = request.initial_volume.as_unit("L")
    temperature = request.initial_temperature.as_unit("K")
    if volume <= 0 or temperature <= 0:
        raise ValueError("Initial volume and temperature must be positive")

    amounts = initial.copy()
    additions_total = np.zeros_like(amounts)
    removed_total = np.zeros_like(amounts)
    time_s = 0.0
    scale = max(initial.max(), 1e-12)
    timeline = [{"time_s": 0.0, "stage_id": "initial", "amounts_mol": dict(zip(ids, amounts.tolist()))}]
    stages: list[dict[str, Any]] = []

    for stage in request.stages:
        before = amounts.copy()
        volume_before = volume
        if stage.temperature is not None:
            temperature = stage.temperature.as_unit("K")
            if temperature <= 0:
                raise ValueError(f"Stage {stage.id} temperature must be positive")
        if stage.volume is not None:
            volume = stage.volume.as_unit("L")
            if volume <= 0:
                raise ValueError(f"Stage {stage.id} volume must be positive")
        stage_additions = np.zeros_like(amounts)
        for addition in stage.additions:
            value = addition.amount.as_unit("mol")
            if value < 0:
                raise ValueError("Protocol additions cannot be negative")
            stage_additions[index[addition.species_id]] += value
        amounts += stage_additions
        additions_total += stage_additions

        stage_removed = np.zeros_like(amounts)
        if stage.operation == "sample":
            stage_removed = amounts * stage.withdraw_fraction
            amounts -= stage_removed
            volume *= 1 - stage.withdraw_fraction
        elif stage.operation == "transfer":
            stage_removed = amounts * (1 - stage.recovery_fraction)
            amounts *= stage.recovery_fraction
            volume *= stage.recovery_fraction
        removed_total += stage_removed
        after_actions = amounts.copy()
        duration_s = stage.duration.as_unit("s")
        extents: dict[str, float] = {}
        samples: list[dict[str, Any]] = []

        if stage.operation == "reaction":
            active_columns = [prepared["reaction_index"][key] for key in stage.active_reactions]
            matrix = prepared["stoichiometry"][:, active_columns]
            metadata = [prepared["reaction_metadata"][column] for column in active_columns]

            def rhs(_time, state):
                current = np.maximum(state[:len(ids)], 0)
                fluxes = np.asarray([
                    _reaction_flux(reaction, current, volume, temperature) for reaction in metadata
                ])
                return np.concatenate((matrix @ fluxes, fluxes))

            state0 = np.concatenate((amounts, np.zeros(len(active_columns))))
            solution = solve_ivp(
                rhs,
                (0, duration_s),
                state0,
                method="BDF",
                t_eval=np.linspace(0, duration_s, request.points_per_reaction_stage),
                rtol=1e-8,
                atol=max(scale * 1e-11, 1e-15),
            )
            if not solution.success or not np.isfinite(solution.y).all():
                raise ValueError(f"Kinetic integration failed in stage {stage.id}: {solution.message}")
            amounts = _clean_amounts(solution.y[:len(ids), -1], scale)
            extents = dict(zip(stage.active_reactions, solution.y[len(ids):, -1].tolist()))
            samples = [
                {"time_s": time_s + float(relative_time), "stage_id": stage.id,
                 "amounts_mol": dict(zip(ids, _clean_amounts(vector[:len(ids)], scale).tolist()))}
                for relative_time, vector in zip(solution.t, solution.y.T)
            ]
            timeline.extend(samples[1:])
        elif duration_s > 0:
            timeline.append({"time_s": time_s + duration_s, "stage_id": stage.id,
                             "amounts_mol": dict(zip(ids, amounts.tolist()))})
        else:
            timeline.append({"time_s": time_s, "stage_id": stage.id,
                             "amounts_mol": dict(zip(ids, amounts.tolist()))})
        time_s += duration_s
        stages.append({
            "id": stage.id,
            "name": stage.name,
            "operation": stage.operation,
            "start_time_s": time_s - duration_s,
            "end_time_s": time_s,
            "temperature_K": temperature,
            "volume_before_L": volume_before,
            "volume_after_L": volume,
            "active_reactions": stage.active_reactions,
            "inventory_before_mol": dict(zip(ids, before.tolist())),
            "additions_mol": dict(zip(ids, stage_additions.tolist())),
            "removed_mol": dict(zip(ids, stage_removed.tolist())),
            "inventory_after_actions_mol": dict(zip(ids, after_actions.tolist())),
            "inventory_final_mol": dict(zip(ids, amounts.tolist())),
            "reaction_extents_mol": extents,
            "note": stage.note,
        })

    expected_atoms = prepared["atomic"] @ (initial + additions_total - removed_total)
    final_atoms = prepared["atomic"] @ amounts
    residual = final_atoms - expected_atoms
    tolerance = max(float(np.max(np.abs(expected_atoms), initial=0)) * 1e-7, 1e-10)
    if np.max(np.abs(residual), initial=0) > tolerance:
        raise ValueError("Whole-protocol material ledger failed atom/charge closure")
    result: dict[str, Any] = {
        "status": "conditional_protocol_simulation",
        "protocol_id": request.id,
        "total_duration_s": time_s,
        "species_identity": prepared["identities"],
        "reactions": prepared["reaction_metadata"],
        "stages": stages,
        "timeline": timeline,
        "final_amounts_mol": dict(zip(ids, amounts.tolist())),
        "material_ledger": {
            "initial_mol": dict(zip(ids, initial.tolist())),
            "added_mol": dict(zip(ids, additions_total.tolist())),
            "removed_mol": dict(zip(ids, removed_total.tolist())),
            "final_mol": dict(zip(ids, amounts.tolist())),
            "conserved_quantities": prepared["conserved"],
            "residual": dict(zip(prepared["conserved"], residual.tolist())),
            "closed": True,
        },
        "scientific_boundary": {
            "recipe_used_for_blind_discovery": False,
            "rate_evidence_is_user_supplied": True,
            "stereochemical_identity_preserved": True,
            "enantioselectivity_inferred": False,
        },
        "warnings": [
            "This is a user-declared staged kinetic protocol, not blind product discovery or an automatically validated laboratory recipe.",
            "A stage activates only its listed reactions. Missing reactions, inhibition, enzyme deactivation, transport, pH and phase behavior remain absent unless explicitly modeled.",
            "R/S and E/Z identities are retained, but different stereoisomer rates require separately sourced stereospecific kinetic evidence and an explicit chiral catalyst/environment.",
            "Additions, samples and transfer losses are explicit ledger events; they are not chemical reactions.",
        ],
    }
    from .sbml_adapter import build_sbml, cross_check_roadrunner
    sbml = build_sbml(request, prepared)
    result["sbml"] = {
        "validation": sbml["validation"],
        "sha256": sbml["sha256"],
        "network_only": True,
        "document": sbml["document"] if request.include_sbml else None,
    }
    result["roadrunner_cross_check"] = (
        cross_check_roadrunner(request, prepared, sbml["document"], result["final_amounts_mol"])
        if request.cross_check_roadrunner else {"status": "not_requested"}
    )
    return result
