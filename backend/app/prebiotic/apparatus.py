"""Closed recirculating gas-liquid apparatus for prebiotic chemistry studies.

This layer handles transport, declared kinetics, interventions, and ledgers. It
does not contain a Miller-Urey product list or a hidden reaction recipe.
"""
from __future__ import annotations

from collections import Counter
import gc
import math
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from scipy.constants import R
from scipy.integrate import solve_ivp

from app.chemistry.general_builder import parse_molecule

from .aqueous_speciation import (
    augment_prepared_with_aqueous_ions,
    equilibrate_aqueous_ion_inventory,
)
from .blind_aqueous import discover_aqueous_candidates
from .boltzmann import collision_set_targets
from .coupled_plasma import solve_apparatus_plasma_snapshot
from .formaldehyde_hydration import (
    METHANEDIOL_AQ,
    augment_prepared_with_methanediol,
    formaldehyde_hydration_enthalpy_J_mol,
    formaldehyde_hydration_parameters,
    source_bounded_formaldehyde_contact,
    validate_diagnostic_formaldehyde_range,
)
from .heavy_gas import (
    PAVAN_HCN_DIAGNOSTIC_REACTIONS,
    augment_prepared_with_gri,
    pavan_internal_conserved_mol,
    run_gri_heavy_species_cycle,
    run_pavan_charged_species_cycle,
)
from .mechanism_overlap import audit_pavan_gri_overlap
from .phase_partition import (
    capture_gas_water_by_flow,
    closed_flask_coupled_henry_coefficients,
    coupled_henry_transport_coefficients,
    equilibrate_gas_water,
    enforce_water_vapor_capacity,
    transfer_gas_water_finite_rate,
)
from .plasma import calculate_rate_coefficient
from .schema import COMPARTMENTS, PrebioticApparatusRequest


def _prepare(request: PrebioticApparatusRequest) -> dict[str, Any]:
    ids = [item.id for item in request.species]
    index = {key: position for position, key in enumerate(ids)}
    inventories: list[Counter[str]] = []
    identities: list[dict[str, Any]] = []
    transfer: list[dict[str, Any]] = []
    hypothetical_sources: list[str] = []
    for species in request.species:
        inventory = Counter(species.composition)
        inventory["charge"] = species.charge
        canonical_smiles = None
        formula = (
            "e-"
            if species.phase_behavior == "electron"
            else "".join(
                f"{element}{count if count != 1 else ''}"
                for element, count in species.composition.items()
            )
        )
        stereo_status = "declared_composition_only"
        stereochemistry: list[str] = []
        if species.smiles is not None:
            molecule = parse_molecule(species.smiles)
            smiles_inventory = Counter(atom.GetSymbol() for atom in molecule.GetAtoms())
            smiles_charge = Chem.GetFormalCharge(molecule)
            if smiles_inventory != Counter(species.composition) or smiles_charge != species.charge:
                raise ValueError(
                    f"Species {species.id} SMILES disagrees with its declared composition or charge"
                )
            parent = Chem.RemoveHs(Chem.Mol(molecule))
            Chem.AssignStereochemistry(parent, cleanIt=True, force=True)
            potential = list(Chem.FindPotentialStereo(parent))
            centers = Chem.FindMolChiralCenters(
                parent, force=True, includeUnassigned=True, useLegacyImplementation=False
            )
            unspecified = any(str(value.specified) == "Unspecified" for value in potential)
            stereo_status = "achiral" if not potential else "unassigned" if unspecified else "specified"
            stereochemistry = [f"atom {position + 1}: {label}" for position, label in centers]
            canonical_smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
            formula = rdMolDescriptors.CalcMolFormula(parent)
        inventories.append(inventory)
        rates: dict[str, float] = {}
        evidence: dict[str, Any] = {}
        for name, quantity in (
            ("evaporation", species.evaporation_rate),
            ("condensation", species.condensation_rate),
            ("wall_loss", species.wall_loss_rate),
        ):
            rates[name] = quantity.as_unit("1/s") if quantity is not None else 0.0
            evidence[name] = quantity.model_dump(mode="json") if quantity is not None else None
            if quantity is not None and quantity.basis == "hypothetical":
                hypothetical_sources.append(f"{species.id}:{name}")
        transfer.append({"rates_s_inv": rates, "evidence": evidence})
        identities.append(
            {
                "id": species.id,
                "name": species.name,
                "composition": species.composition,
                "formula": formula,
                "charge": species.charge,
                "isomeric_smiles": canonical_smiles,
                "stereo_status": stereo_status,
                "stereochemistry": stereochemistry,
                "phase_behavior": species.phase_behavior,
                "enantiomer_of": species.enantiomer_of,
                "identity_note": species.identity_note,
            }
        )
    conserved = sorted({key for inventory in inventories for key in inventory})
    atomic = np.asarray(
        [[inventory.get(key, 0) for inventory in inventories] for key in conserved], dtype=float
    )
    stoichiometry = np.zeros((len(ids), len(request.reactions)), dtype=float)
    reactions: list[dict[str, Any]] = []
    for column, reaction in enumerate(request.reactions):
        for species_id, coefficient in reaction.stoichiometry.items():
            stoichiometry[index[species_id], column] = coefficient
        residual = atomic @ stoichiometry[:, column]
        if np.max(np.abs(residual), initial=0) > 1e-12:
            detail = {
                key: float(value)
                for key, value in zip(conserved, residual, strict=True)
                if abs(value) > 1e-12
            }
            raise ValueError(f"Reaction {reaction.id} does not conserve atoms/charge: {detail}")
        collision_result = None
        if reaction.rate_law.kind == "electron_collision":
            target_id = reaction.rate_law.target_species_id
            electron_id = reaction.rate_law.electron_species_id
            reactant_orders = {index[target_id]: 1, index[electron_id]: 1}
            molecularity = 2
            target_unit = "liter / mole / second"
            if request.spark.mode == "cross_section_eedf":
                collision_result = calculate_rate_coefficient(
                    reaction.rate_law.collision_dataset, request.spark.eedf
                )
                rate = collision_result["rate_coefficient_L_mol_s"]
            else:
                rate = 0.0
        else:
            reactant_orders = {
                index[species_id]: -coefficient
                for species_id, coefficient in reaction.stoichiometry.items()
                if coefficient < 0
            }
            molecularity = sum(reactant_orders.values())
            target_unit = (
                "1/s"
                if molecularity == 1
                else f"liter ** {molecularity - 1} / mole ** {molecularity - 1} / second"
            )
            rate = reaction.rate_law.rate_constant.as_unit(target_unit)
            if rate <= 0:
                raise ValueError(f"Reaction {reaction.id} rate constant must be positive")
        energy_cost = (
            reaction.rate_law.energy_cost.as_unit("J/mol")
            if reaction.rate_law.energy_cost is not None
            else None
        )
        if (
            reaction.rate_law.rate_constant is not None
            and reaction.rate_law.rate_constant.basis == "hypothetical"
        ):
            hypothetical_sources.append(f"{reaction.id}:rate_constant")
        if reaction.rate_law.energy_cost is not None and reaction.rate_law.energy_cost.basis == "hypothetical":
            hypothetical_sources.append(f"{reaction.id}:energy_cost")
        reactions.append(
            {
                "id": reaction.id,
                "name": reaction.name,
                "phase": reaction.phase,
                "kind": reaction.rate_law.kind,
                "reactant_orders": reactant_orders,
                "molecularity": molecularity,
                "rate_constant": rate,
                "rate_constant_unit": target_unit,
                "rate_evidence": (
                    reaction.rate_law.rate_constant.model_dump(mode="json")
                    if reaction.rate_law.rate_constant is not None
                    else {
                        "collision_dataset": reaction.rate_law.collision_dataset.model_dump(
                            mode="json", exclude_none=True
                        ),
                        "eedf": (
                            request.spark.eedf.model_dump(mode="json", exclude_none=True)
                            if request.spark.eedf is not None
                            else None
                        ),
                    }
                ),
                "collision_rate": collision_result,
                "target_species_id": reaction.rate_law.target_species_id,
                "electron_species_id": reaction.rate_law.electron_species_id,
                "energy_cost_J_mol": energy_cost,
            }
        )
    return {
        "ids": ids,
        "index": index,
        "identities": identities,
        "transfer": transfer,
        "conserved": conserved,
        "atomic": atomic,
        "stoichiometry": stoichiometry,
        "reactions": reactions,
        "hypothetical_sources": sorted(set(hypothetical_sources)),
    }


def _pass_fraction(flow_L_s: float, duration_s: float, gas_volume_L: float) -> float:
    if flow_L_s <= 0 or duration_s <= 0:
        return 0.0
    return 1.0 - math.exp(-flow_L_s * duration_s / gas_volume_L)


def _move_fraction(source: np.ndarray, destination: np.ndarray, fractions: np.ndarray) -> np.ndarray:
    moved = np.maximum(source, 0) * np.clip(fractions, 0, 1)
    source -= moved
    destination += moved
    return moved


def _apply_coupled_henry_transfer(
    compartments: dict[str, np.ndarray],
    index: dict[str, int],
    transport_report: dict[str, Any] | None,
    *,
    max_negative_roundoff_mol: float = 1e-11,
) -> None:
    if not transport_report:
        return
    if max_negative_roundoff_mol < 0:
        raise ValueError("Coupled transport roundoff limit cannot be negative")
    clipped: dict[str, float] = {}
    for species, delta in transport_report.get("boiler_delta_mol", {}).items():
        position = index[species]
        compartments["boiler_liquid"][position] += float(delta)
        if compartments["boiler_liquid"][position] < -max_negative_roundoff_mol:
            raise ValueError("Coupled boiler transfer produced negative inventory")
        if compartments["boiler_liquid"][position] < 0:
            clipped[species] = -float(compartments["boiler_liquid"][position])
        compartments["boiler_liquid"][position] = max(
            compartments["boiler_liquid"][position], 0.0
        )
    transport_report["boiler_negative_roundoff_clipped_mol"] = clipped
    for species, captured in transport_report.get("captured_mol", {}).items():
        compartments["condensate_trap"][index[species]] += float(captured)


def _integrate_reactions(
    amounts: np.ndarray,
    volume_L: float,
    duration_s: float,
    reaction_columns: list[int],
    prepared: dict[str, Any],
    multipliers: dict[str, float],
    scale: float,
    energy_budget_J: float | None = None,
) -> tuple[np.ndarray, dict[str, float], bool]:
    if duration_s <= 0 or not reaction_columns:
        return amounts, {}, False
    matrix = prepared["stoichiometry"][:, reaction_columns]
    metadata = [prepared["reactions"][column] for column in reaction_columns]

    def rhs(_time, state):
        current = np.maximum(state[: len(amounts)], 0)
        concentrations = current / volume_L
        fluxes = []
        for reaction in metadata:
            if any(current[position] <= 0 for position in reaction["reactant_orders"]):
                fluxes.append(0.0)
                continue
            activity = math.prod(
                concentrations[position] ** order
                for position, order in reaction["reactant_orders"].items()
            )
            fluxes.append(
                reaction["rate_constant"]
                * activity
                * volume_L
                * multipliers.get(reaction["kind"], 1.0)
            )
        vector = np.asarray(fluxes, dtype=float)
        return np.concatenate((matrix @ vector, vector))

    state0 = np.concatenate((amounts, np.zeros(len(reaction_columns))))
    energy_costs = np.asarray(
        [float(reaction["energy_cost_J_mol"] or 0.0) for reaction in metadata], dtype=float
    )
    can_enforce_energy = (
        energy_budget_J is not None
        and energy_budget_J > 0
        and np.any(energy_costs > 0)
        and all(
            reaction["energy_cost_J_mol"] is not None
            for reaction in metadata
            if reaction["kind"] in {"electron_impact_reduced", "electron_collision"}
        )
    )

    def energy_remaining(_time, state):
        return float(energy_budget_J - np.dot(np.maximum(state[len(amounts) :], 0), energy_costs))

    energy_remaining.terminal = True
    energy_remaining.direction = -1
    solution = solve_ivp(
        rhs,
        (0, duration_s),
        state0,
        method="BDF",
        t_eval=None if can_enforce_energy else [duration_s],
        rtol=1e-8,
        atol=max(scale * 1e-11, 1e-15),
        events=energy_remaining if can_enforce_energy else None,
    )
    if not solution.success or not np.isfinite(solution.y).all():
        raise ValueError(f"Apparatus kinetic integration failed: {solution.message}")
    final_state = solution.y[:, -1]
    energy_limited = bool(can_enforce_energy and solution.t[-1] < duration_s * (1 - 1e-10))
    if energy_limited:
        remaining_s = duration_s - float(solution.t[-1])
        original_multipliers = multipliers
        multipliers = {
            **multipliers,
            "electron_impact_reduced": 0.0,
            "electron_collision": 0.0,
        }
        continuation = solve_ivp(
            rhs,
            (0, remaining_s),
            final_state,
            method="BDF",
            t_eval=[remaining_s],
            rtol=1e-8,
            atol=max(scale * 1e-11, 1e-15),
        )
        multipliers = original_multipliers
        if not continuation.success or not np.isfinite(continuation.y).all():
            raise ValueError(f"Energy-limited apparatus integration failed: {continuation.message}")
        final_state = continuation.y[:, -1]
    final = final_state[: len(amounts)]
    tolerance = max(scale * 1e-8, 1e-13)
    if np.min(final, initial=0) < -tolerance:
        raise ValueError("Apparatus integration produced a materially negative amount")
    final = np.maximum(final, 0)
    extents = {
        reaction["id"]: float(value)
        for reaction, value in zip(metadata, final_state[len(amounts) :], strict=True)
    }
    return final, extents, energy_limited


def _compartment_dict(ids: list[str], amounts: np.ndarray) -> dict[str, float]:
    return {species_id: float(value) for species_id, value in zip(ids, amounts, strict=True)}


def _sum_record_maps(records: list[dict[str, Any]], field: str) -> dict[str, float]:
    keys = {key for record in records for key in record.get(field, {})}
    return {
        key: sum(float(record.get(field, {}).get(key, 0.0)) for record in records)
        for key in sorted(keys)
    }


def _aggregate_integration_steps(
    records: list[dict[str, Any]], substeps_per_cycle: int
) -> list[dict[str, Any]]:
    """Collapse split operator steps back to declared apparatus cycles."""

    if substeps_per_cycle <= 1:
        return records
    aggregated = []
    additive_maps = (
        "additions_mol",
        "sampled_mol",
        "evaporated_mol",
        "condensed_mol",
        "returned_mol",
        "wall_deposition_mol",
        "reaction_extents_mol",
        "local_ion_clipping_projection_mol",
        "local_unprojected_ion_roundoff_clip_mol",
    )
    for start in range(0, len(records), substeps_per_cycle):
        group = records[start : start + substeps_per_cycle]
        if len(group) != substeps_per_cycle:
            raise ValueError("Continuous-overlap integration ended with a partial cycle")
        combined = dict(group[-1])
        combined["cycle"] = len(aggregated) + 1
        combined.pop("substep", None)
        combined["integration_substeps"] = substeps_per_cycle
        combined["start_time_s"] = float(group[0]["start_time_s"])
        combined["end_time_s"] = float(group[-1]["end_time_s"])
        combined["intervention_notes"] = [
            note for record in group for note in record.get("intervention_notes", [])
        ]
        for field in additive_maps:
            if any(field in record for record in group):
                combined[field] = _sum_record_maps(group, field)
        for field in ("deposited_electrical_energy_J", "deposited_energy_used_J"):
            if any(field in record for record in group):
                combined[field] = sum(float(record.get(field, 0.0)) for record in group)
        if any("energy_limited" in record for record in group):
            combined["energy_limited"] = any(
                bool(record.get("energy_limited", False)) for record in group
            )
        if any("solver_status" in record for record in group):
            statuses = [str(record.get("solver_status", "not_run")) for record in group]
            combined["integration_solver_statuses"] = statuses
            combined["integration_solver_errors"] = [
                record.get("solver_error") for record in group
            ]
            combined["solver_failed_in_any_substep"] = any(
                status.endswith("_failed") or status == "not_run"
                for status in statuses
            )
        if any(record.get("local_atom_residual_mol") is not None for record in group):
            worst = max(
                (record for record in group if record.get("local_atom_residual_mol") is not None),
                key=lambda record: float(record["local_atom_residual_mol"]),
            )
            combined["worst_local_atom_residual_substep"] = worst.get("substep")
            combined["worst_local_atom_residual_mol"] = worst["local_atom_residual_mol"]
            combined["worst_local_raw_atom_residual_mol"] = worst.get(
                "local_raw_atom_residual_mol"
            )
            combined["worst_local_clipped_negative_species_mol"] = worst.get(
                "local_clipped_negative_species_mol"
            )
            combined["worst_local_atom_residual_by_element_mol"] = worst.get(
                "local_atom_residual_by_element_mol"
            )
        aggregated.append(combined)
    return aggregated


def simulate_prebiotic_apparatus(
    request: PrebioticApparatusRequest,
    *,
    diagnostic_pathway_species_each_substep: tuple[str, ...] = (),
    diagnostic_gri_reaction_multipliers: tuple[tuple[int, float], ...] = (),
    diagnostic_literature_h2cn_h_to_hcn: bool = False,
    diagnostic_hebrard_h2cn_n_branches: bool = False,
    diagnostic_hebrard_h2cn_h_isomer_branches: bool = False,
    diagnostic_hebrard_hnc_h_isomerization: bool = False,
    diagnostic_charged_hcn_final_rate_each_substep: bool = False,
    diagnostic_charged_hcn_integrated_each_substep: bool = False,
    diagnostic_pavan_reaction_multipliers: tuple[tuple[int, float], ...] = (),
) -> dict[str, Any]:
    if len(set(diagnostic_pathway_species_each_substep)) != len(
        diagnostic_pathway_species_each_substep
    ):
        raise ValueError("Diagnostic pathway species must be unique")
    pavan_multiplier_indices = [
        index for index, _ in diagnostic_pavan_reaction_multipliers
    ]
    if len(pavan_multiplier_indices) != len(set(pavan_multiplier_indices)):
        raise ValueError("Charged diagnostic reaction multipliers must be unique")
    if any(index not in PAVAN_HCN_DIAGNOSTIC_REACTIONS for index in pavan_multiplier_indices):
        raise ValueError("Unsupported charged diagnostic reaction index")
    if any(
        not math.isfinite(factor) or factor < 0
        for _, factor in diagnostic_pavan_reaction_multipliers
    ):
        raise ValueError("Charged diagnostic reaction factors must be finite and nonnegative")
    if diagnostic_hebrard_h2cn_h_isomer_branches and request.spark.mode != "boltzmann_eedf":
        raise ValueError("HNC branch diagnostic requires the Boltzmann/GRI apparatus")
    if diagnostic_hebrard_hnc_h_isomerization and not diagnostic_hebrard_h2cn_h_isomer_branches:
        raise ValueError("HNC + H isomerization requires the opt-in HNC species branch")
    prepared = _prepare(request)
    if request.spark.mode == "boltzmann_eedf":
        augment_prepared_with_gri(
            prepared,
            diagnostic_hnc_temperature_K=(
                request.apparatus.discharge_temperature.as_unit("K")
                if diagnostic_hebrard_h2cn_h_isomer_branches else None
            ),
        )
        if request.spark.boltzmann.electron_number_density.basis == "hypothetical":
            prepared["hypothetical_sources"].append("spark:electron_number_density")
            prepared["hypothetical_sources"] = sorted(set(prepared["hypothetical_sources"]))
    if request.apparatus.aqueous_acid_base_model == "hcn_nh3_water_equilibrium":
        augment_prepared_with_aqueous_ions(prepared)
    if request.apparatus.closed_flask_formaldehyde_hydration != "off":
        augment_prepared_with_methanediol(prepared)
    ids, index = prepared["ids"], prepared["index"]
    compartments = {
        name: np.zeros(len(ids), dtype=float) for name in COMPARTMENTS
    }
    initial = np.zeros(len(ids), dtype=float)
    for item in request.initial_inventory:
        amount = item.amount.as_unit("mol")
        if amount < 0:
            raise ValueError("Initial inventory amounts cannot be negative")
        compartments[item.compartment][index[item.species_id]] += amount
        initial[index[item.species_id]] += amount
    if initial.sum() <= 0:
        raise ValueError("Initial inventory must contain a positive amount")
    if request.apparatus.operation_mode == "closed_stirred_flask":
        if any(np.any(compartments[place] > 0) for place in ("condensate_trap", "collected_condensate")):
            raise ValueError("Closed stirred flask initial inventory cannot use condenser compartments")
        if request.spark.mode == "boltzmann_eedf":
            collision_targets = set(collision_set_targets(request.spark.boltzmann.collision_set))
            unsupported_feed = sorted({
                item.species_id
                for item in request.initial_inventory
                if item.compartment == "gas_loop"
                and item.amount.as_unit("mol") > 0
                and item.species_id not in collision_targets
                and prepared["identities"][index[item.species_id]]["charge"] == 0
            })
            if unsupported_feed:
                raise ValueError(
                    "Closed-flask Boltzmann spark has no collision-set coverage for initial "
                    f"neutral gas feed: {', '.join(unsupported_feed)}"
                )
    acid_base_updates = 0
    acid_base_maximum_residual_mol = 0.0
    hydration_minimum_pH_ideal: float | None = None
    hydration_extrapolation_used = False

    def equilibrate_liquids(branch: dict[str, np.ndarray]) -> dict[str, dict[str, Any]]:
        nonlocal acid_base_updates, acid_base_maximum_residual_mol
        reports: dict[str, dict[str, Any]] = {}
        if request.apparatus.aqueous_acid_base_model == "off":
            return reports
        for place, temperature in (
            ("boiler_liquid", request.apparatus.boiler_temperature.as_unit("K")),
            ("condensate_trap", request.apparatus.condenser_temperature.as_unit("K")),
            ("collected_condensate", request.apparatus.condenser_temperature.as_unit("K")),
        ):
            branch[place], report = equilibrate_aqueous_ion_inventory(
                branch[place], prepared, temperature
            )
            reports[place] = report
            if report["status"] == "liquid_ideal_acid_base_equilibrium_estimated":
                acid_base_updates += 1
                acid_base_maximum_residual_mol = max(
                    acid_base_maximum_residual_mol,
                    float(report["ledger_maximum_residual_mol"]),
                )
        return reports

    equilibrate_liquids(compartments)
    charged_compartments = (
        {name: values.copy() for name, values in compartments.items()}
        if request.spark.mode == "boltzmann_eedf"
        else None
    )
    charged_internal_gas_state: dict[str, float] = {}
    charged_internal_sampled_state: dict[str, float] = {}
    charged_sampled_total = np.zeros(len(ids), dtype=float)
    scale = max(float(initial.max()), 1e-12)
    apparatus = request.apparatus
    hydration_diagnostic_extrapolation = (
        apparatus.closed_flask_formaldehyde_hydration
        == "winkelman_2002_pH4_to_5_extrapolation_diagnostic"
    )
    liquid_volume_L = apparatus.boiler_liquid_volume.as_unit("L")
    gas_volume_L = apparatus.gas_loop_volume.as_unit("L")
    discharge_zone_volume_L = (
        apparatus.discharge_zone_volume.as_unit("L")
        if apparatus.discharge_zone_volume is not None
        else gas_volume_L
    )
    discharge_zone_fraction = min(discharge_zone_volume_L / gas_volume_L, 1.0)
    trap_volume_L = apparatus.condensate_trap_volume.as_unit("L")
    contact_water_holdup_L = (
        apparatus.condenser_contact_water_holdup.as_unit("L")
        if apparatus.condenser_contact_water_holdup is not None else None
    )
    pressure_Pa = apparatus.pressure.as_unit("Pa")
    flow_L_s = apparatus.gas_circulation_rate.as_unit("L/s")
    finite_closed_contact = (
        apparatus.operation_mode == "closed_stirred_flask"
        and apparatus.closed_flask_gas_water_transfer_model == "finite_rate_contact"
    )
    coupled_closed_contact = (
        finite_closed_contact
        and apparatus.closed_flask_finite_contact_solver == "coupled_bdf"
    )
    closed_flask_kla_per_s = (
        apparatus.closed_flask_liquid_side_kla.as_unit("1/s")
        if apparatus.closed_flask_liquid_side_kla is not None else None
    )
    if (
        apparatus.closed_flask_liquid_side_kla is not None
        and apparatus.closed_flask_liquid_side_kla.basis == "hypothetical"
    ):
        prepared["hypothetical_sources"].append("apparatus:closed_flask_liquid_side_kla")
        prepared["hypothetical_sources"] = sorted(set(prepared["hypothetical_sources"]))
    operation_substeps = (
        apparatus.integration_substeps_per_cycle
        if apparatus.operation_mode in {"continuous_overlap", "closed_stirred_flask"}
        else 1
    )
    declared_boil_s = apparatus.boil_duration.as_unit("s")
    declared_discharge_s = apparatus.discharge_duration.as_unit("s")
    declared_condense_s = apparatus.condense_duration.as_unit("s")
    boil_s = declared_boil_s / operation_substeps
    discharge_s = declared_discharge_s / operation_substeps
    # The recirculating loop uses mean residence in the localized discharge.
    # A closed mixed flask retains its gas for the entire interval; only the
    # plasma source inventory is scaled by the active-zone volume fraction.
    heavy_chemistry_residence_fraction = (
        1.0 if apparatus.operation_mode == "closed_stirred_flask"
        else discharge_zone_fraction
    )
    effective_heavy_chemistry_residence_s = discharge_s * heavy_chemistry_residence_fraction
    declared_effective_heavy_chemistry_residence_s = (
        declared_discharge_s * heavy_chemistry_residence_fraction
    )
    condense_s = declared_condense_s / operation_substeps
    if (
        apparatus.discharge_zone_volume is not None
        and apparatus.discharge_zone_volume.basis == "hypothetical"
    ):
        prepared["hypothetical_sources"].append("apparatus:discharge_zone_volume")
        prepared["hypothetical_sources"] = sorted(set(prepared["hypothetical_sources"]))
    if (
        apparatus.condenser_contact_water_holdup is not None
        and apparatus.condenser_contact_water_holdup.basis == "hypothetical"
    ):
        prepared["hypothetical_sources"].append("apparatus:condenser_contact_water_holdup")
        prepared["hypothetical_sources"] = sorted(set(prepared["hypothetical_sources"]))
    cycle_s = boil_s + discharge_s + condense_s
    declared_cycle_s = declared_boil_s + declared_discharge_s + declared_condense_s
    integration_steps = apparatus.cycles * operation_substeps
    boil_pass = _pass_fraction(flow_L_s, boil_s, gas_volume_L)
    discharge_pass = _pass_fraction(flow_L_s, discharge_s, gas_volume_L)
    condense_pass = _pass_fraction(flow_L_s, condense_s, gas_volume_L)
    # A closed stirred flask has no circulation pass, but gas-phase chemistry
    # still evolves throughout the flask. Electron-impact channels are instead
    # limited to the declared discharge-zone fraction of a well-mixed gas.
    gas_chemistry_fraction = (
        1.0 if apparatus.operation_mode == "closed_stirred_flask" else discharge_pass
    )
    electron_exposure_fraction = (
        discharge_zone_fraction
        if apparatus.operation_mode == "closed_stirred_flask"
        else discharge_pass
    )
    rates = [item["rates_s_inv"] for item in prepared["transfer"]]
    evaporation = np.asarray([item["evaporation"] for item in rates])
    condensation = np.asarray([item["condensation"] for item in rates])
    wall_loss = np.asarray([item["wall_loss"] for item in rates])
    interventions: dict[int, list[Any]] = {}
    for intervention in request.interventions:
        interventions.setdefault(intervention.cycle, []).append(intervention)
    additions_total = np.zeros(len(ids), dtype=float)
    sampled_total = np.zeros(len(ids), dtype=float)
    extents_total = {reaction.id: 0.0 for reaction in request.reactions}
    electron_extent_total = {
        reaction.id: 0.0
        for reaction in request.reactions
        if reaction.rate_law.kind in {"electron_impact_reduced", "electron_collision"}
    }
    gas_columns = [
        column for column, reaction in enumerate(prepared["reactions"]) if reaction["phase"] == "gas"
    ]
    liquid_columns = [
        column for column, reaction in enumerate(prepared["reactions"]) if reaction["phase"] == "liquid"
    ]
    duty_fraction = 0.0
    electrical_estimate_per_cycle_J = 0.0
    deposited_per_cycle_J = 0.0
    deposited_energy_basis = "spark_off"
    if request.spark.mode in {"reduced_order", "cross_section_eedf", "boltzmann_eedf"}:
        duty_fraction = (
            request.spark.pulse_width.as_unit("s") * request.spark.frequency.as_unit("Hz")
        )
        electrical_estimate_per_cycle_J = (
            request.spark.voltage.as_unit("V")
            * request.spark.current.as_unit("A")
            * discharge_s
            * duty_fraction
            * float(request.spark.coupling_efficiency)
        )
        if request.spark.calorimetric_deposited_energy_per_cycle is not None:
            deposited_per_cycle_J = (
                request.spark.calorimetric_deposited_energy_per_cycle.as_unit("J")
                / operation_substeps
            )
            deposited_energy_basis = "calorimetric_measurement"
            if request.spark.calorimetric_deposited_energy_per_cycle.basis == "hypothetical":
                prepared["hypothetical_sources"].append(
                    "spark:calorimetric_deposited_energy_per_cycle"
                )
                prepared["hypothetical_sources"] = sorted(
                    set(prepared["hypothetical_sources"])
                )
        else:
            deposited_per_cycle_J = electrical_estimate_per_cycle_J
            deposited_energy_basis = "declared_electrical_estimate"
    cycles: list[dict[str, Any]] = []
    plasma_cycles: list[dict[str, Any]] = []
    heavy_species_cycles: list[dict[str, Any]] = []
    partition_cycles: list[dict[str, Any]] = []
    charged_branch_cycles: list[dict[str, Any]] = []
    charged_hcn_final_rate_snapshots: list[dict[str, Any]] = []
    charged_hcn_integrated_snapshots: list[dict[str, Any]] = []
    charged_gas_pressures_Pa: list[float] = []
    any_energy_limited = False

    for integration_step in range(1, integration_steps + 1):
        cycle = (integration_step - 1) // operation_substeps + 1
        substep = (integration_step - 1) % operation_substeps + 1
        cycle_additions = np.zeros(len(ids), dtype=float)
        cycle_sampled = np.zeros(len(ids), dtype=float)
        intervention_notes: list[str] = []
        for intervention in interventions.get(cycle, []) if substep == 1 else []:
            for addition in intervention.additions:
                value = addition.amount.as_unit("mol")
                if value < 0:
                    raise ValueError("Intervention additions cannot be negative")
                position = index[addition.species_id]
                compartments[addition.compartment][position] += value
                if charged_compartments is not None:
                    charged_compartments[addition.compartment][position] += value
                cycle_additions[position] += value
            if intervention.sample_compartment is not None:
                removed = (
                    compartments[intervention.sample_compartment] * intervention.sample_fraction
                )
                compartments[intervention.sample_compartment] -= removed
                cycle_sampled += removed
                if charged_compartments is not None:
                    charged_removed = (
                        charged_compartments[intervention.sample_compartment]
                        * intervention.sample_fraction
                    )
                    charged_compartments[intervention.sample_compartment] -= charged_removed
                    charged_sampled_total += charged_removed
                    if intervention.sample_compartment == "gas_loop":
                        for species, amount in list(charged_internal_gas_state.items()):
                            removed_internal = amount * intervention.sample_fraction
                            charged_internal_gas_state[species] = amount - removed_internal
                            charged_internal_sampled_state[species] = (
                                charged_internal_sampled_state.get(species, 0.0)
                                + removed_internal
                            )
            if intervention.note:
                intervention_notes.append(intervention.note)
        additions_total += cycle_additions
        sampled_total += cycle_sampled
        primary_aqueous_reports = equilibrate_liquids(compartments)
        if charged_compartments is not None:
            charged_aqueous_reports = equilibrate_liquids(charged_compartments)
        else:
            charged_aqueous_reports = {}

        evaporated = _move_fraction(
            compartments["boiler_liquid"],
            compartments["gas_loop"],
            (1 - np.exp(-evaporation * boil_s)) * boil_pass,
        )
        boiler_partition = None
        water_vapor_capacity = None
        condenser_water_vapor_capacity = None
        if request.spark.mode == "boltzmann_eedf":
            (
                compartments["gas_loop"],
                compartments["boiler_liquid"],
                water_vapor_capacity,
            ) = enforce_water_vapor_capacity(
                compartments["gas_loop"],
                compartments["boiler_liquid"],
                prepared,
                apparatus.boiler_temperature.as_unit("K"),
                gas_volume_L,
            )
            if apparatus.condenser_phase_transfer_model != "coupled_flow_transport" and not finite_closed_contact:
                (
                    compartments["gas_loop"],
                    compartments["boiler_liquid"],
                    boiler_partition,
                ) = equilibrate_gas_water(
                    compartments["gas_loop"],
                    compartments["boiler_liquid"],
                    prepared,
                    apparatus.boiler_temperature.as_unit("K"),
                    float(water_vapor_capacity["ideal_total_pressure_Pa"]),
                    gas_volume_L,
                )
        charged_boiler_partition = None
        charged_water_vapor_capacity = None
        if charged_compartments is not None:
            _move_fraction(
                charged_compartments["boiler_liquid"],
                charged_compartments["gas_loop"],
                (1 - np.exp(-evaporation * boil_s)) * boil_pass,
            )
            (
                charged_compartments["gas_loop"],
                charged_compartments["boiler_liquid"],
                charged_water_vapor_capacity,
            ) = enforce_water_vapor_capacity(
                charged_compartments["gas_loop"],
                charged_compartments["boiler_liquid"],
                prepared,
                apparatus.boiler_temperature.as_unit("K"),
                gas_volume_L,
            )
            if apparatus.condenser_phase_transfer_model != "coupled_flow_transport" and not finite_closed_contact:
                (
                    charged_compartments["gas_loop"],
                    charged_compartments["boiler_liquid"],
                    charged_boiler_partition,
                ) = equilibrate_gas_water(
                    charged_compartments["gas_loop"],
                    charged_compartments["boiler_liquid"],
                    prepared,
                    apparatus.boiler_temperature.as_unit("K"),
                    float(charged_water_vapor_capacity["ideal_total_pressure_Pa"]),
                    gas_volume_L,
                )
        if apparatus.closed_flask_formaldehyde_hydration != "off":
            # The water-vapor cap above may change the liquid-water amount.
            # Refresh pH at the actual contact boundary before selecting the
            # source-bounded hydration rate and free-molecule Henry ratio.
            primary_aqueous_reports = equilibrate_liquids(compartments)
            if charged_compartments is not None:
                charged_aqueous_reports = equilibrate_liquids(charged_compartments)
        gas_before = compartments["gas_loop"].copy()
        coupled_primary_transport = None
        coupled_charged_transport = None
        primary_formaldehyde_hydration = None
        charged_formaldehyde_hydration = None
        if (
            apparatus.condenser_phase_transfer_model == "coupled_flow_transport"
            and request.spark.mode == "boltzmann_eedf"
        ):
            water_position = index.get("H2O")
            if water_position is not None:
                def coefficients_for(branch: dict[str, np.ndarray]) -> dict[str, dict[str, float]]:
                    return coupled_henry_transport_coefficients(
                        prepared,
                        gas_volume_L,
                        flow_L_s,
                        apparatus.boiler_temperature.as_unit("K"),
                        apparatus.condenser_temperature.as_unit("K"),
                        float(branch["boiler_liquid"][water_position]),
                        float(branch["boiler_liquid"][water_position]
                              + branch["condensate_trap"][water_position]),
                        float(contact_water_holdup_L),
                        boil_s,
                        condense_s,
                        effective_heavy_chemistry_residence_s,
                    )
                coupled_primary_transport = coefficients_for(compartments)
                if charged_compartments is not None:
                    coupled_charged_transport = coefficients_for(charged_compartments)
        elif coupled_closed_contact:
            water_position = index.get("H2O")
            if water_position is not None:
                def closed_coefficients(
                    branch: dict[str, np.ndarray], aqueous_reports: dict[str, dict[str, Any]],
                ) -> tuple[dict[str, dict[str, float]], dict[str, Any] | None]:
                    nonlocal hydration_minimum_pH_ideal, hydration_extrapolation_used
                    water_amount = float(branch["boiler_liquid"][water_position])
                    coefficients = closed_flask_coupled_henry_coefficients(
                        prepared, gas_volume_L,
                        apparatus.boiler_temperature.as_unit("K"),
                        water_amount,
                        float(closed_flask_kla_per_s),
                    )
                    hydration = None
                    if apparatus.closed_flask_formaldehyde_hydration != "off":
                        acid_report = aqueous_reports.get("boiler_liquid", {})
                        if acid_report.get("status") != "liquid_ideal_acid_base_equilibrium_estimated":
                            raise ValueError("Formaldehyde hydration requires calculable boiler-liquid pH")
                        current_pH = float(acid_report["pH_ideal_aqueous_model"])
                        current_temperature = apparatus.boiler_temperature.as_unit("K")
                        if hydration_diagnostic_extrapolation:
                            validate_diagnostic_formaldehyde_range(current_temperature, current_pH)
                        hydration = source_bounded_formaldehyde_contact(
                            temperature_K=current_temperature,
                            pH=current_pH,
                            gas_volume_L=gas_volume_L,
                            liquid_water_mol=water_amount,
                            liquid_side_kla_per_s=float(closed_flask_kla_per_s),
                            allow_extrapolation=hydration_diagnostic_extrapolation,
                        )
                        hydration_minimum_pH_ideal = (
                            current_pH if hydration_minimum_pH_ideal is None
                            else min(hydration_minimum_pH_ideal, current_pH)
                        )
                        hydration_extrapolation_used |= bool(hydration["explicit_extrapolation"])
                        if "CH2O" not in coefficients:
                            raise ValueError("Formaldehyde hydration requires a CH2O gas-transfer species")
                        ratio = float(hydration["free_liquid_to_gas_ratio"])
                        coefficients["CH2O"]["boiler_liquid_to_gas_equilibrium_ratio"] = ratio
                        coefficients["CH2O"]["boiler_relaxation_rate_per_chemistry_s"] = (
                            float(closed_flask_kla_per_s) * (1.0 + ratio)
                        )
                    return coefficients, hydration
                coupled_primary_transport, primary_formaldehyde_hydration = closed_coefficients(
                    compartments, primary_aqueous_reports,
                )
                if charged_compartments is not None:
                    coupled_charged_transport, charged_formaldehyde_hydration = closed_coefficients(
                        charged_compartments, charged_aqueous_reports,
                    )
        plasma_snapshot = None
        heavy_species_snapshot = None
        charged_plasma_snapshot = None
        charged_shadow_snapshot = None
        if request.spark.mode == "boltzmann_eedf":
            discharge_zone_gas = gas_before * discharge_zone_fraction
            plasma_snapshot = solve_apparatus_plasma_snapshot(
                request,
                prepared["identities"],
                discharge_zone_gas,
                discharge_zone_volume_L,
                discharge_s,
                duty_fraction,
                deposited_per_cycle_J,
            )
            plasma_snapshot.update(
                {
                    "discharge_zone_volume_L": discharge_zone_volume_L,
                    "discharge_zone_fraction": discharge_zone_fraction,
                    "whole_loop_mean_residence_time_s": effective_heavy_chemistry_residence_s,
                    "residence_time_method": (
                        "closed mixed flask: whole gas remains for the full integration interval"
                        if apparatus.operation_mode == "closed_stirred_flask"
                        else "well-mixed loop occupancy: discharge stage duration multiplied "
                        "by discharge-zone volume / gas-loop volume"
                    ),
                }
            )
            plasma_cycles.append(
                {"cycle": cycle, "substep": substep, **plasma_snapshot}
            )
            if charged_compartments is not None:
                charged_gas_before = charged_compartments["gas_loop"].copy()
                charged_plasma_snapshot = solve_apparatus_plasma_snapshot(
                    request,
                    prepared["identities"],
                    charged_gas_before * discharge_zone_fraction,
                    discharge_zone_volume_L,
                    discharge_s,
                    duty_fraction,
                    deposited_per_cycle_J,
                )
                charged_plasma_snapshot.update(
                    {
                        "discharge_zone_volume_L": discharge_zone_volume_L,
                        "discharge_zone_fraction": discharge_zone_fraction,
                        "whole_loop_mean_residence_time_s": effective_heavy_chemistry_residence_s,
                        "residence_time_method": (
                            "closed mixed flask: whole charged branch remains for the full interval"
                            if apparatus.operation_mode == "closed_stirred_flask"
                            else "well-mixed propagated charged branch: discharge stage duration "
                            "multiplied by discharge-zone volume / gas-loop volume"
                        ),
                    }
                )
                charged_initial_state = {
                    ids[position]: float(amount)
                    for position, amount in enumerate(charged_gas_before)
                    if amount > 0
                }
                charged_initial_state.update(charged_internal_gas_state)
                charged_gas_pressures_Pa.append(
                    sum(max(float(amount), 0.0) for amount in charged_initial_state.values())
                    * R
                    * apparatus.discharge_temperature.as_unit("K")
                    / (gas_volume_L * 1e-3)
                )
                try:
                    supplied_elements = prepared["atomic"] @ (initial + additions_total)
                    declared_absent_elements = {
                        element
                        for element, amount in zip(prepared["conserved"], supplied_elements, strict=True)
                        if element not in {"charge", "E"} and amount == 0.0
                    }
                    charged_shadow_snapshot = run_pavan_charged_species_cycle(
                        charged_gas_before,
                        prepared,
                        charged_plasma_snapshot,
                        apparatus.discharge_temperature.as_unit("K"),
                        gas_volume_L,
                        effective_heavy_chemistry_residence_s,
                        deposited_per_cycle_J,
                        initial_state_mol=charged_initial_state,
                        coupled_henry_transport=coupled_charged_transport,
                        boiler_species_mol={
                            species: float(charged_compartments["boiler_liquid"][index[species]])
                            for species in (coupled_charged_transport or {})
                        } | ({
                            species: float(charged_compartments["boiler_liquid"][index[species]])
                            for species in (METHANEDIOL_AQ, "H2O")
                        } if charged_formaldehyde_hydration is not None else {}),
                        formaldehyde_hydration=charged_formaldehyde_hydration,
                        declared_absent_elements=declared_absent_elements,
                        integrated_pathway_species=(
                            ("HCN", "HCNN", "CH", "CH2", "CH2(S)", "CH2OH", "CH3O", "CH2O", "HCO", "CO")
                            if diagnostic_charged_hcn_integrated_each_substep else ()
                        ),
                        diagnostic_reaction_multipliers=diagnostic_pavan_reaction_multipliers,
                    )
                except Exception as error:  # pragma: no cover - data-dependent guard
                    if charged_formaldehyde_hydration is not None:
                        raise ValueError(
                            "Closed-flask charged hydration branch failed at "
                            f"cycle {cycle}, substep {substep}: {error}"
                        ) from error
                    charged_shadow_snapshot = {
                        "status": "propagated_pavan_charged_shadow_failed",
                        "error": f"{type(error).__name__}: {error}"[:500],
                        "formed_mol": {},
                        "net_change_mol": {},
                        "final_state_mol": charged_initial_state,
                        "warnings": [
                            "The propagated charged branch failed this cycle and retained its pre-discharge state; the primary inventory was not changed."
                        ],
                        "scientific_firewall": {
                            "primary_inventory_mutated_by_shadow": False,
                            "quantitative_yield_claim_allowed": False,
                        },
                    }
                charged_shadow_snapshot["branch_plasma_snapshot"] = {
                    "status": charged_plasma_snapshot.get("status"),
                    "mole_fractions": charged_plasma_snapshot.get("mole_fractions", {}),
                    "covered_gas_fraction": charged_plasma_snapshot.get(
                        "covered_gas_fraction", 0.0
                    ),
                    "mean_electron_energy_eV": charged_plasma_snapshot.get(
                        "mean_electron_energy_eV"
                    ),
                    "reactive_channel_count": charged_plasma_snapshot.get(
                        "reactive_channel_count", 0
                    ),
                    "balanced_reactive_channel_count": charged_plasma_snapshot.get(
                        "balanced_reactive_channel_count", 0
                    ),
                }
                charged_shadow_snapshot["propagation_cycle"] = cycle
                charged_shadow_snapshot["propagation_substep"] = substep
                charged_shadow_snapshot["carried_from_previous_cycle"] = (
                    integration_step > 1
                )
                final_charged_state = charged_shadow_snapshot.get(
                    "final_state_mol", charged_initial_state
                )
                for species in charged_shadow_snapshot.get("state_species_names", []):
                    if species in index:
                        charged_compartments["gas_loop"][index[species]] = 0.0
                for species, amount in final_charged_state.items():
                    if species in index:
                        charged_compartments["gas_loop"][index[species]] = float(amount)
                charged_internal_gas_state = {
                    species: float(amount)
                    for species, amount in final_charged_state.items()
                    if species not in index and float(amount) > 0
                }
                _apply_coupled_henry_transfer(
                    charged_compartments,
                    index,
                    charged_shadow_snapshot.get("coupled_henry_transport"),
                    max_negative_roundoff_mol=(1e-18 if coupled_closed_contact else 1e-11),
                )
            gas_before_heavy = gas_before.copy()
            boiler_before_heavy = compartments["boiler_liquid"].copy()
            gas_before, heavy_species_snapshot = run_gri_heavy_species_cycle(
                gas_before,
                prepared,
                plasma_snapshot,
                apparatus.discharge_temperature.as_unit("K"),
                gas_volume_L,
                effective_heavy_chemistry_residence_s,
                deposited_per_cycle_J,
                propagated_charged_shadow=charged_shadow_snapshot,
                source_forcing=apparatus.primary_source_forcing,
                include_diagnostic_shadows=substep == operation_substeps,
                coupled_henry_transport=coupled_primary_transport,
                boiler_species_mol={
                    species: float(compartments["boiler_liquid"][index[species]])
                    for species in (coupled_primary_transport or {})
                } | ({
                    species: float(compartments["boiler_liquid"][index[species]])
                    for species in (METHANEDIOL_AQ, "H2O")
                } if primary_formaldehyde_hydration is not None else {}),
                formaldehyde_hydration=primary_formaldehyde_hydration,
                integrated_pathway_species=tuple(dict.fromkeys(
                    (
                        ("N2", "N", "NO", "CH3", "H2CN", "HCN", "CH2O")
                        if integration_step == 1 else ()
                    ) + diagnostic_pathway_species_each_substep
                )),
                diagnostic_reaction_multipliers=diagnostic_gri_reaction_multipliers,
                diagnostic_literature_h2cn_h_to_hcn=diagnostic_literature_h2cn_h_to_hcn,
                diagnostic_hebrard_h2cn_n_branches=diagnostic_hebrard_h2cn_n_branches,
                diagnostic_hebrard_h2cn_h_isomer_branches=(
                    diagnostic_hebrard_h2cn_h_isomer_branches
                ),
                diagnostic_hebrard_hnc_h_isomerization=(
                    diagnostic_hebrard_hnc_h_isomerization
                ),
            )
            _apply_coupled_henry_transfer(
                compartments,
                index,
                heavy_species_snapshot.get("coupled_henry_transport"),
                max_negative_roundoff_mol=(1e-18 if coupled_closed_contact else 1e-11),
            )
            if apparatus.operation_mode == "closed_stirred_flask":
                # A coupled contact transfers equal amounts from gas to water;
                # check their combined local ledger, not gas alone.
                local_residual = prepared["atomic"] @ (
                    gas_before + compartments["boiler_liquid"]
                    - gas_before_heavy - boiler_before_heavy
                )
                heavy_species_snapshot["local_primary_conserved_residual_mol"] = {
                    key: float(value)
                    for key, value in zip(prepared["conserved"], local_residual, strict=True)
                }
            if (
                heavy_species_snapshot.get("coupled_henry_transport") is not None
                and not coupled_closed_contact
            ):
                boiler_partition = {
                    "status": "coupled_flow_transport_with_gri_bdf",
                    "transferred_to_water_mol": heavy_species_snapshot[
                        "coupled_henry_transport"
                    ]["boiler_delta_mol"],
                }
            heavy_species_cycles.append(
                {"cycle": cycle, "substep": substep, **heavy_species_snapshot}
            )
        gas_multipliers = {
            "mass_action": gas_chemistry_fraction,
            "electron_impact_reduced": (
                electron_exposure_fraction * duty_fraction
                if request.spark.mode == "reduced_order" else 0.0
            ),
            "electron_collision": (
                electron_exposure_fraction * duty_fraction
                if request.spark.mode == "cross_section_eedf"
                else 0.0
            ),
        }
        gas_after, gas_extents, cycle_energy_limited = _integrate_reactions(
            gas_before,
            gas_volume_L,
            discharge_s,
            gas_columns,
            prepared,
            gas_multipliers,
            scale,
            deposited_per_cycle_J if request.spark.mode != "off" else None,
        )
        any_energy_limited = any_energy_limited or cycle_energy_limited
        compartments["gas_loop"] = gas_after
        for reaction_id, extent in gas_extents.items():
            extents_total[reaction_id] += extent
            if reaction_id in electron_extent_total:
                electron_extent_total[reaction_id] += extent

        if apparatus.operation_mode == "closed_stirred_flask":
            pre_spark_partition = boiler_partition
            (
                compartments["gas_loop"],
                compartments["boiler_liquid"],
                water_vapor_capacity,
            ) = enforce_water_vapor_capacity(
                compartments["gas_loop"], compartments["boiler_liquid"], prepared,
                apparatus.boiler_temperature.as_unit("K"), gas_volume_L,
            )
            if coupled_closed_contact:
                post_spark_partition = {
                    "status": "closed_stirred_flask_finite_contact_coupled_with_gri_bdf",
                    "transferred_to_water_mol": (
                        (heavy_species_snapshot or {}).get("coupled_henry_transport") or {}
                    ).get("boiler_delta_mol", {}),
                }
            elif finite_closed_contact:
                (
                    compartments["gas_loop"],
                    compartments["boiler_liquid"],
                    post_spark_partition,
                ) = transfer_gas_water_finite_rate(
                    compartments["gas_loop"], compartments["boiler_liquid"], prepared,
                    apparatus.boiler_temperature.as_unit("K"),
                    float(water_vapor_capacity["ideal_total_pressure_Pa"]), gas_volume_L,
                    discharge_s, float(closed_flask_kla_per_s),
                )
            else:
                (
                    compartments["gas_loop"],
                    compartments["boiler_liquid"],
                    post_spark_partition,
                ) = equilibrate_gas_water(
                    compartments["gas_loop"], compartments["boiler_liquid"], prepared,
                    apparatus.boiler_temperature.as_unit("K"),
                    float(water_vapor_capacity["ideal_total_pressure_Pa"]), gas_volume_L,
                )
            boiler_partition = {
                "status": (
                    "closed_stirred_flask_finite_contact_coupled_bdf"
                    if coupled_closed_contact else
                    "closed_stirred_flask_finite_rate_contact"
                    if finite_closed_contact else
                    "closed_stirred_flask_step_equilibrium_upper_bound"
                ),
                "before_spark": pre_spark_partition,
                "after_spark": post_spark_partition,
            }
            if charged_compartments is not None:
                (
                    charged_compartments["gas_loop"],
                    charged_compartments["boiler_liquid"],
                    charged_water_vapor_capacity,
                ) = enforce_water_vapor_capacity(
                    charged_compartments["gas_loop"], charged_compartments["boiler_liquid"],
                    prepared, apparatus.boiler_temperature.as_unit("K"), gas_volume_L,
                )
                if coupled_closed_contact:
                    charged_boiler_partition = {
                        "status": "closed_stirred_flask_finite_contact_coupled_with_pavan_bdf",
                        "transferred_to_water_mol": (
                            (charged_shadow_snapshot or {}).get("coupled_henry_transport") or {}
                        ).get("boiler_delta_mol", {}),
                    }
                elif finite_closed_contact:
                    (
                        charged_compartments["gas_loop"],
                        charged_compartments["boiler_liquid"],
                        charged_boiler_partition,
                    ) = transfer_gas_water_finite_rate(
                        charged_compartments["gas_loop"], charged_compartments["boiler_liquid"],
                        prepared, apparatus.boiler_temperature.as_unit("K"),
                        float(charged_water_vapor_capacity["ideal_total_pressure_Pa"]), gas_volume_L,
                        discharge_s, float(closed_flask_kla_per_s),
                    )
                else:
                    (
                        charged_compartments["gas_loop"],
                        charged_compartments["boiler_liquid"],
                        charged_boiler_partition,
                    ) = equilibrate_gas_water(
                        charged_compartments["gas_loop"], charged_compartments["boiler_liquid"],
                        prepared, apparatus.boiler_temperature.as_unit("K"),
                        float(charged_water_vapor_capacity["ideal_total_pressure_Pa"]), gas_volume_L,
                    )
            partition_cycles.append({
                "cycle": cycle, "substep": substep, "boiler": boiler_partition,
                "water_vapor_capacity": water_vapor_capacity,
                "condenser_water_vapor_capacity": None, "condenser": None,
            })

        condensed = _move_fraction(
            compartments["gas_loop"],
            compartments["condensate_trap"],
            (1 - np.exp(-condensation * condense_s)) * condense_pass,
        )
        condenser_partition = None
        if apparatus.operation_mode != "closed_stirred_flask" and (
            request.spark.mode == "boltzmann_eedf"
            or apparatus.condenser_phase_transfer_model in {
                "flow_limited_capture", "coupled_flow_transport"
            }
        ):
            (
                compartments["gas_loop"],
                compartments["condensate_trap"],
                condenser_water_vapor_capacity,
            ) = enforce_water_vapor_capacity(
                compartments["gas_loop"],
                compartments["condensate_trap"],
                prepared,
                apparatus.condenser_temperature.as_unit("K"),
                gas_volume_L,
            )
            if (
                apparatus.condenser_phase_transfer_model == "coupled_flow_transport"
                and request.spark.mode == "boltzmann_eedf"
            ):
                condenser_partition = {
                    "status": "coupled_with_gri_bdf",
                    "transferred_to_water_mol": (
                        (heavy_species_snapshot or {}).get("coupled_henry_transport") or {}
                    ).get("captured_mol", {}),
                    "limitations": [
                        "Mean-field simultaneous exchange is not spatially resolved and its contact hold-up is hypothetical."
                    ],
                }
            elif apparatus.condenser_phase_transfer_model in {
                "flow_limited_capture", "coupled_flow_transport"
            }:
                water_position = index.get("H2O")
                available_water_mol = (
                    float(compartments["boiler_liquid"][water_position]
                          + compartments["condensate_trap"][water_position])
                    if water_position is not None else 0.0
                )
                (
                    compartments["gas_loop"],
                    compartments["condensate_trap"],
                    condenser_partition,
                ) = capture_gas_water_by_flow(
                    compartments["gas_loop"],
                    compartments["condensate_trap"],
                    prepared,
                    apparatus.condenser_temperature.as_unit("K"),
                    float(condenser_water_vapor_capacity["ideal_total_pressure_Pa"]),
                    gas_volume_L,
                    flow_L_s,
                    condense_s,
                    float(contact_water_holdup_L),
                    available_water_mol,
                )
            else:
                (
                    compartments["gas_loop"],
                    compartments["condensate_trap"],
                    condenser_partition,
                ) = equilibrate_gas_water(
                    compartments["gas_loop"],
                    compartments["condensate_trap"],
                    prepared,
                    apparatus.condenser_temperature.as_unit("K"),
                    float(condenser_water_vapor_capacity["ideal_total_pressure_Pa"]),
                    gas_volume_L,
                )
            partition_cycles.append(
                {
                    "cycle": cycle,
                    "substep": substep,
                    "boiler": boiler_partition,
                    "water_vapor_capacity": water_vapor_capacity,
                    "condenser_water_vapor_capacity": condenser_water_vapor_capacity,
                    "condenser": condenser_partition,
                }
            )
        charged_condenser_partition = None
        charged_condenser_water_vapor_capacity = None
        charged_condensed = np.zeros(len(ids), dtype=float)
        if charged_compartments is not None and apparatus.operation_mode != "closed_stirred_flask":
            charged_condensed = _move_fraction(
                charged_compartments["gas_loop"],
                charged_compartments["condensate_trap"],
                (1 - np.exp(-condensation * condense_s)) * condense_pass,
            )
            (
                charged_compartments["gas_loop"],
                charged_compartments["condensate_trap"],
                charged_condenser_water_vapor_capacity,
            ) = enforce_water_vapor_capacity(
                charged_compartments["gas_loop"],
                charged_compartments["condensate_trap"],
                prepared,
                apparatus.condenser_temperature.as_unit("K"),
                gas_volume_L,
            )
            if (
                apparatus.condenser_phase_transfer_model == "coupled_flow_transport"
                and request.spark.mode == "boltzmann_eedf"
            ):
                charged_condenser_partition = {
                    "status": "coupled_with_pavan_bdf",
                    "transferred_to_water_mol": (
                        (charged_shadow_snapshot or {}).get("coupled_henry_transport") or {}
                    ).get("captured_mol", {}),
                }
            elif apparatus.condenser_phase_transfer_model in {
                "flow_limited_capture", "coupled_flow_transport"
            }:
                water_position = index.get("H2O")
                charged_available_water_mol = (
                    float(charged_compartments["boiler_liquid"][water_position]
                          + charged_compartments["condensate_trap"][water_position])
                    if water_position is not None else 0.0
                )
                (
                    charged_compartments["gas_loop"],
                    charged_compartments["condensate_trap"],
                    charged_condenser_partition,
                ) = capture_gas_water_by_flow(
                    charged_compartments["gas_loop"],
                    charged_compartments["condensate_trap"],
                    prepared,
                    apparatus.condenser_temperature.as_unit("K"),
                    float(charged_condenser_water_vapor_capacity["ideal_total_pressure_Pa"]),
                    gas_volume_L,
                    flow_L_s,
                    condense_s,
                    float(contact_water_holdup_L),
                    charged_available_water_mol,
                )
            else:
                (
                    charged_compartments["gas_loop"],
                    charged_compartments["condensate_trap"],
                    charged_condenser_partition,
                ) = equilibrate_gas_water(
                    charged_compartments["gas_loop"],
                    charged_compartments["condensate_trap"],
                    prepared,
                    apparatus.condenser_temperature.as_unit("K"),
                    float(charged_condenser_water_vapor_capacity["ideal_total_pressure_Pa"]),
                    gas_volume_L,
                )
        return_fractions = np.full(
            len(ids), apparatus.condensate_nonwater_return_fraction, dtype=float
        )
        if "H2O" in index:
            return_fractions[index["H2O"]] = apparatus.condensate_return_fraction
        returned = compartments["condensate_trap"] * return_fractions
        collected = compartments["condensate_trap"] - returned
        compartments["condensate_trap"] -= returned + collected
        compartments["boiler_liquid"] += returned
        compartments["collected_condensate"] += collected
        charged_returned = np.zeros(len(ids), dtype=float)
        charged_collected = np.zeros(len(ids), dtype=float)
        if charged_compartments is not None:
            charged_returned = charged_compartments["condensate_trap"] * return_fractions
            charged_collected = charged_compartments["condensate_trap"] - charged_returned
            charged_compartments["condensate_trap"] -= charged_returned + charged_collected
            charged_compartments["boiler_liquid"] += charged_returned
            charged_compartments["collected_condensate"] += charged_collected

        liquid_after, liquid_extents, _ = _integrate_reactions(
            compartments["boiler_liquid"].copy(),
            liquid_volume_L,
            cycle_s,
            liquid_columns,
            prepared,
            {"mass_action": 1.0, "electron_impact_reduced": 0.0, "electron_collision": 0.0},
            scale,
        )
        compartments["boiler_liquid"] = liquid_after
        primary_end_aqueous_reports = equilibrate_liquids(compartments)
        if charged_compartments is not None:
            charged_end_aqueous_reports = equilibrate_liquids(charged_compartments)
        else:
            charged_end_aqueous_reports = {}
        if apparatus.closed_flask_formaldehyde_hydration != "off":
            for branch_name, reports in (
                ("primary", primary_end_aqueous_reports),
                ("charged", charged_end_aqueous_reports),
            ):
                if branch_name == "charged" and charged_compartments is None:
                    continue
                report = reports.get("boiler_liquid", {})
                if report.get("status") != "liquid_ideal_acid_base_equilibrium_estimated":
                    raise ValueError("Source-bounded hydration needs calculable aqueous pH")
                pH = float(report["pH_ideal_aqueous_model"])
                hydration_minimum_pH_ideal = (
                    pH if hydration_minimum_pH_ideal is None
                    else min(hydration_minimum_pH_ideal, pH)
                )
                try:
                    current_temperature = apparatus.boiler_temperature.as_unit("K")
                    if hydration_diagnostic_extrapolation:
                        validate_diagnostic_formaldehyde_range(current_temperature, pH)
                    rate_report = formaldehyde_hydration_parameters(
                        current_temperature, pH,
                        allow_extrapolation=hydration_diagnostic_extrapolation,
                    )
                    hydration_extrapolation_used |= bool(rate_report["explicit_extrapolation"])
                except ValueError as error:
                    raise ValueError(
                        f"{branch_name} formaldehyde hydration left "
                        f"{'diagnostic' if hydration_diagnostic_extrapolation else 'measured'} domain "
                        f"after cycle {cycle}, substep {substep}: pH={pH:.5g}"
                    ) from error
        for reaction_id, extent in liquid_extents.items():
            extents_total[reaction_id] += extent

        wall_from_gas = _move_fraction(
            compartments["gas_loop"],
            compartments["wall"],
            1 - np.exp(-wall_loss * cycle_s),
        )
        wall_from_liquid = _move_fraction(
            compartments["boiler_liquid"],
            compartments["wall"],
            1 - np.exp(-wall_loss * cycle_s),
        )
        if charged_compartments is not None:
            _move_fraction(
                charged_compartments["gas_loop"],
                charged_compartments["wall"],
                1 - np.exp(-wall_loss * cycle_s),
            )
            _move_fraction(
                charged_compartments["boiler_liquid"],
                charged_compartments["wall"],
                1 - np.exp(-wall_loss * cycle_s),
            )
            if diagnostic_charged_hcn_final_rate_each_substep:
                charged_hcn_final_rate_snapshots.append({
                    "cycle": cycle,
                    "substep": substep,
                    "solver_status": (
                        charged_shadow_snapshot.get("status")
                        if charged_shadow_snapshot is not None else "not_run"
                    ),
                    "rate_basis": "instantaneous net rates at charged BDF endpoint; not integrated reaction extents",
                    "hcn": (
                        charged_shadow_snapshot.get(
                            "final_instantaneous_pathway_diagnostics", {}
                        ).get("HCN")
                        if charged_shadow_snapshot is not None else None
                    ),
                })
            if diagnostic_charged_hcn_integrated_each_substep:
                charged_hcn_integrated_snapshots.append({
                    "cycle": cycle,
                    "substep": substep,
                    "solver_status": (
                        charged_shadow_snapshot.get("status")
                        if charged_shadow_snapshot is not None else "not_run"
                    ),
                    "diagnostic": (
                        charged_shadow_snapshot.get("integrated_pathway_diagnostics")
                        if charged_shadow_snapshot is not None else None
                    ),
                    "electron_sources": (
                        charged_shadow_snapshot.get("source_extents")
                        if charged_shadow_snapshot is not None else None
                    ),
                })
            charged_branch_cycles.append(
                {
                    "cycle": cycle,
                    "substep": substep,
                    "start_time_s": (integration_step - 1) * cycle_s,
                    "end_time_s": integration_step * cycle_s,
                    "deposited_electrical_energy_J": electrical_estimate_per_cycle_J,
                    "deposited_energy_used_J": deposited_per_cycle_J,
                    "deposited_energy_basis": deposited_energy_basis,
                    "solver_status": (
                        charged_shadow_snapshot.get("status")
                        if charged_shadow_snapshot is not None
                        else "not_run"
                    ),
                    "solver_error": (
                        charged_shadow_snapshot.get("error")
                        if charged_shadow_snapshot is not None
                        else None
                    ),
                    "final_charge_residual_mol_e": (
                        charged_shadow_snapshot.get("final_charge_residual_mol_e")
                        if charged_shadow_snapshot is not None
                        else None
                    ),
                    "local_atom_residual_mol": (
                        charged_shadow_snapshot.get("max_atom_residual_mol")
                        if charged_shadow_snapshot is not None
                        else None
                    ),
                    "local_raw_atom_residual_mol": (
                        charged_shadow_snapshot.get("max_raw_atom_residual_mol")
                        if charged_shadow_snapshot is not None
                        else None
                    ),
                    "local_atom_residual_by_element_mol": (
                        charged_shadow_snapshot.get("atom_residual_by_element_mol")
                        if charged_shadow_snapshot is not None
                        else None
                    ),
                    "local_clipped_negative_species_mol": (
                        charged_shadow_snapshot.get("clipped_negative_species_mol")
                        if charged_shadow_snapshot is not None
                        else None
                    ),
                    "local_ion_clipping_projection_mol": (
                        charged_shadow_snapshot.get("conservative_ion_clipping_projection_mol", {})
                        if charged_shadow_snapshot is not None
                        else {}
                    ),
                    "local_unprojected_ion_roundoff_clip_mol": (
                        charged_shadow_snapshot.get("unprojected_ion_roundoff_clip_mol", {})
                        if charged_shadow_snapshot is not None
                        else {}
                    ),
                    "mean_electron_energy_eV": (
                        charged_plasma_snapshot.get("mean_electron_energy_eV")
                        if charged_plasma_snapshot is not None
                        else None
                    ),
                    "compartments_mol": {
                        name: _compartment_dict(ids, values)
                        for name, values in charged_compartments.items()
                    },
                }
            )
        cycles.append(
            {
                "cycle": cycle,
                "substep": substep,
                "start_time_s": (integration_step - 1) * cycle_s,
                "end_time_s": integration_step * cycle_s,
                "intervention_notes": intervention_notes,
                "additions_mol": _compartment_dict(ids, cycle_additions),
                "sampled_mol": _compartment_dict(ids, cycle_sampled),
                "evaporated_mol": _compartment_dict(ids, evaporated),
                "condensed_mol": _compartment_dict(ids, condensed),
                "returned_mol": _compartment_dict(ids, returned),
                "collected_mol": _compartment_dict(ids, collected),
                "wall_deposition_mol": _compartment_dict(ids, wall_from_gas + wall_from_liquid),
                "reaction_extents_mol": {**gas_extents, **liquid_extents},
                "deposited_electrical_energy_J": electrical_estimate_per_cycle_J,
                "deposited_energy_used_J": deposited_per_cycle_J,
                "deposited_energy_basis": deposited_energy_basis,
                "energy_limited": cycle_energy_limited,
                "plasma_snapshot": plasma_snapshot,
                "heavy_species_snapshot": heavy_species_snapshot,
                "phase_partition": {
                    "boiler": boiler_partition,
                    "water_vapor_capacity": water_vapor_capacity,
                    "condenser_water_vapor_capacity": condenser_water_vapor_capacity,
                    "condenser": condenser_partition,
                }
                if condenser_partition is not None or apparatus.operation_mode == "closed_stirred_flask"
                else None,
                "compartments_mol": {
                    name: _compartment_dict(ids, values) for name, values in compartments.items()
                },
            }
        )
        # Cantera ReactorNet/Solution graphs contain Python-to-C++ reference
        # cycles. Reclaim completed solver graphs during long resolution scans
        # instead of waiting for the next automatic collection threshold.
        if request.spark.mode == "boltzmann_eedf" and integration_step % 8 == 0:
            gc.collect()

    # Retain integration-step snapshots for a post-discovery HCN occupancy
    # diagnostic. They must never feed a reaction or change the state update.
    hcn_integration_steps = cycles
    charged_hcn_integration_steps = charged_branch_cycles
    cycles = _aggregate_integration_steps(cycles, operation_substeps)
    charged_branch_cycles = _aggregate_integration_steps(
        charged_branch_cycles, operation_substeps
    )

    final = sum(compartments.values(), np.zeros(len(ids), dtype=float))
    expected_atoms = prepared["atomic"] @ (initial + additions_total - sampled_total)
    final_atoms = prepared["atomic"] @ final
    residual = final_atoms - expected_atoms
    tolerances = np.maximum(np.abs(expected_atoms) * 1e-12, 1e-12)
    closed = bool(np.all(np.abs(residual) <= tolerances))
    if not closed:
        failures = {
            key: {"residual_mol": float(error), "tolerance_mol": float(tolerance)}
            for key, error, tolerance in zip(prepared["conserved"], residual, tolerances, strict=True)
            if abs(error) > tolerance
        }
        worst_heavy = max(
            heavy_species_cycles,
            key=lambda item: max(
                (abs(value) for value in item.get("local_primary_conserved_residual_mol", {}).values()),
                default=0.0,
            ),
            default={},
        )
        raise ValueError(
            "Whole-apparatus material ledger failed atom/charge closure: "
            f"{failures}; worst_primary_heavy_substep="
            f"{ {key: worst_heavy.get(key) for key in ('cycle', 'substep', 'local_primary_conserved_residual_mol')} }; "
            f"worst_primary_solver="
            f"{ {key: (worst_heavy.get('solver') or {}).get(key) for key in ('negative_species_clipped_mol', 'raw_gas_maximum_element_residual_mol')} }; "
            f"maximum_aqueous_acid_base_residual_mol={acid_base_maximum_residual_mol:.6g}"
        )

    active_liquid = compartments["boiler_liquid"] + compartments["condensate_trap"]
    observed_liquid = active_liquid + compartments["collected_condensate"]
    aqueous_discovery = discover_aqueous_candidates(observed_liquid, prepared)
    # Keep the known-reference layer out of candidate generation.  The import
    # and evidence read occur only after discovery has been serialized/frozen.
    from .post_discovery import evaluate_frozen_aqueous_discovery

    water_index = prepared["index"].get("H2O")
    water_mol = float(active_liquid[water_index]) if water_index is not None else 0.0
    # The reaction concentration basis is the water actually present, not the
    # nominal capacity of the boiler and trap.  At the final/cold sampling
    # state, 18.01528 mL per mol is an explicit density approximation.
    effective_aqueous_volume_L = (
        water_mol * 0.01801528
        if water_mol > 0
        else liquid_volume_L + trap_volume_L
    )
    post_discovery_validation = evaluate_frozen_aqueous_discovery(
        aqueous_discovery,
        active_liquid,
        prepared,
        effective_aqueous_volume_L,
        integration_steps * cycle_s,
        (
            min(
                apparatus.condenser_temperature.as_unit("K"),
                apparatus.boiler_temperature.as_unit("K"),
            ),
            max(
                apparatus.condenser_temperature.as_unit("K"),
                apparatus.boiler_temperature.as_unit("K"),
            ),
        ),
        cycles,
    )
    from .historical_benchmark import evaluate_historical_benchmark

    final_liquid_mol = _compartment_dict(ids, observed_liquid)
    historical_benchmark = evaluate_historical_benchmark(
        aqueous_discovery,
        post_discovery_validation["discovery_artifact_sha256"],
        _compartment_dict(ids, initial),
        final_liquid_mol,
        gas_volume_L,
        integration_steps * cycle_s,
    )

    energy_cost_complete = all(
        reaction["energy_cost_J_mol"] is not None
        for reaction in prepared["reactions"]
        if reaction["kind"] in {"electron_impact_reduced", "electron_collision"}
    )
    declared_channel_energy_J = sum(
        electron_extent_total[reaction["id"]] * float(reaction["energy_cost_J_mol"] or 0)
        for reaction in prepared["reactions"]
        if reaction["kind"] in {"electron_impact_reduced", "electron_collision"}
    )
    electrical_estimate_J = electrical_estimate_per_cycle_J * integration_steps
    deposited_energy_J = deposited_per_cycle_J * integration_steps
    from .precursor_benchmark import evaluate_precursor_energy_benchmark

    simulated_gas_pressures_Pa = [
        float(item["initial_pressure_Pa_from_inventory"])
        for item in heavy_species_cycles
        if item.get("initial_pressure_Pa_from_inventory") is not None
    ]
    simulated_pressure_range_Pa = (
        [min(simulated_gas_pressures_Pa), max(simulated_gas_pressures_Pa)]
        if simulated_gas_pressures_Pa
        else None
    )
    charged_pressure_range_Pa = (
        [min(charged_gas_pressures_Pa), max(charged_gas_pressures_Pa)]
        if charged_gas_pressures_Pa
        else None
    )
    precursor_energy_benchmark = evaluate_precursor_energy_benchmark(
        _compartment_dict(ids, initial),
        final_liquid_mol,
        deposited_energy_J,
        cycles=cycles,
        energy_basis=deposited_energy_basis,
        declared_reference_pressure_Pa=pressure_Pa,
        simulated_gas_pressure_range_Pa=simulated_pressure_range_Pa,
        calorimetric_evidence=(
            request.spark.calorimetric_deposited_energy_per_cycle.model_dump(
                mode="json", exclude_none=True
            )
            if request.spark.calorimetric_deposited_energy_per_cycle is not None
            else None
        ),
    )
    propagated_charged_report = None
    if charged_compartments is not None:
        charged_final_ground = sum(
            charged_compartments.values(), np.zeros(len(ids), dtype=float)
        )
        charged_expected_vector = prepared["atomic"] @ (
            initial + additions_total - charged_sampled_total
        )
        charged_final_vector = prepared["atomic"] @ charged_final_ground
        internal_final = pavan_internal_conserved_mol(charged_internal_gas_state)
        internal_sampled = pavan_internal_conserved_mol(charged_internal_sampled_state)
        charged_expected_by_key = {
            key: float(value) - float(internal_sampled.get(key, 0.0))
            for key, value in zip(
                prepared["conserved"], charged_expected_vector, strict=True
            )
        }
        charged_final_by_key = {
            key: float(value) + float(internal_final.get(key, 0.0))
            for key, value in zip(prepared["conserved"], charged_final_vector, strict=True)
        }
        charged_residual_by_key = {
            key: charged_final_by_key[key] - charged_expected_by_key[key]
            for key in prepared["conserved"]
        }
        charged_tolerances = {
            key: max(abs(value) * 1e-12, 1e-11)
            for key, value in charged_expected_by_key.items()
        }
        charged_closed = all(
            abs(charged_residual_by_key[key]) <= charged_tolerances[key]
            for key in charged_expected_by_key
        )
        if not charged_closed:
            failures = {
                key: {
                    "residual_mol": charged_residual_by_key[key],
                    "tolerance_mol": charged_tolerances[key],
                }
                for key in charged_expected_by_key
                if abs(charged_residual_by_key[key]) > charged_tolerances[key]
            }
            worst_cycle = max(
                charged_branch_cycles,
                key=lambda item: float(item.get("worst_local_atom_residual_mol", item.get("local_atom_residual_mol") or 0.0)),
                default={},
            )
            raise ValueError(
                "Propagated charge-balanced Pavan branch failed atom/charge closure: "
                f"{failures}; internal_ions_mol="
                f"{ {key: value for key, value in charged_internal_gas_state.items() if key in ('e', 'N2+', 'CH4+', 'H2O+', 'O2+', 'CO2+')} }; "
                f"worst_local_cycle="
                f"{ {key: worst_cycle.get(key) for key in ('cycle', 'worst_local_atom_residual_substep', 'worst_local_atom_residual_mol', 'worst_local_raw_atom_residual_mol', 'worst_local_clipped_negative_species_mol', 'worst_local_atom_residual_by_element_mol')} }"
            )
        charged_liquid = (
            charged_compartments["boiler_liquid"]
            + charged_compartments["condensate_trap"]
            + charged_compartments["collected_condensate"]
        )
        charged_precursor_benchmark = evaluate_precursor_energy_benchmark(
            _compartment_dict(ids, initial),
            _compartment_dict(ids, charged_liquid),
            deposited_energy_J,
            cycles=charged_branch_cycles,
            energy_basis=deposited_energy_basis,
            declared_reference_pressure_Pa=pressure_Pa,
            simulated_gas_pressure_range_Pa=charged_pressure_range_Pa,
            calorimetric_evidence=(
                request.spark.calorimetric_deposited_energy_per_cycle.model_dump(
                    mode="json", exclude_none=True
                )
                if request.spark.calorimetric_deposited_energy_per_cycle is not None
                else None
            ),
        )
        propagated_charged_report = {
            "status": "propagated_charge_balanced_pavan_branch_complete",
            "mechanism_overlap": audit_pavan_gri_overlap(),
            "diagnostic_reaction_multipliers": [
                {
                    "reaction_index": index,
                    "equation": PAVAN_HCN_DIAGNOSTIC_REACTIONS[index],
                    "factor": factor,
                }
                for index, factor in diagnostic_pavan_reaction_multipliers
            ],
            "diagnostic_rate_stress_is_validated_kinetics": False,
            "cycles_completed": len(charged_branch_cycles),
            "failed_substep_count": sum(
                sum(
                    status.endswith("_failed") or status == "not_run"
                    for status in cycle.get(
                        "integration_solver_statuses", [cycle.get("solver_status", "not_run")]
                    )
                )
                for cycle in charged_branch_cycles
            ),
            "failed_substep_examples": [
                {"cycle": cycle["cycle"], "substep": offset + 1, "status": status,
                 "error": (cycle.get("integration_solver_errors") or [cycle.get("solver_error")])[offset]}
                for cycle in charged_branch_cycles
                for offset, status in enumerate(
                    cycle.get("integration_solver_statuses", [cycle.get("solver_status", "not_run")])
                )
                if status.endswith("_failed") or status == "not_run"
            ][:5],
            "conservative_ion_clipping_projection_total_mol": sum(
                sum(float(amount) for amount in cycle.get("local_ion_clipping_projection_mol", {}).values())
                for cycle in charged_branch_cycles
            ),
            "unprojected_ion_roundoff_clip_total_mol": sum(
                sum(float(amount) for amount in cycle.get("local_unprojected_ion_roundoff_clip_mol", {}).values())
                for cycle in charged_branch_cycles
            ),
            "cycle_summaries": charged_branch_cycles,
            "final_compartments_mol": {
                name: _compartment_dict(ids, values)
                for name, values in charged_compartments.items()
            },
            "internal_gas_state_mol": charged_internal_gas_state,
            "material_ledger": {
                "expected_conserved_mol": charged_expected_by_key,
                "final_conserved_mol": charged_final_by_key,
                "residual": charged_residual_by_key,
                "tolerance_mol_by_conserved_quantity": charged_tolerances,
                "maximum_absolute_residual_mol": max(
                    (abs(value) for value in charged_residual_by_key.values()), default=0.0
                ),
                "closed": charged_closed,
            },
            "precursor_energy_benchmark": charged_precursor_benchmark,
            "simulated_gas_pressure_range_Pa": charged_pressure_range_Pa,
            "warnings": [
                "This is an alternative propagated Pavan branch and does not mutate or vote with the primary GRI inventory.",
                "Its own gas composition drives a separate self-consistent EEDF on every cycle.",
                "Only ground-state species shared with the apparatus ledger undergo reviewed gas/water partitioning; Pavan-only ions and excited states remain in the gas state.",
                "The packaged Pavan branch reuses GRI's complete HCN reaction set and cannot supply independent HCN kinetic confidence or a quantitative yield claim.",
                "The installed mechanism includes CH + N2 to HCNN but has no NCN species; low-temperature prompt-nitrogen branching is not validated.",
            ],
            "scientific_firewall": {
                "known_products_accessible_to_solver": True,
                "historical_yields_accessible_to_solver": False,
                "primary_inventory_mutated_by_shadow": False,
                "comparison_is_independent_confidence": False,
                "quantitative_yield_claim_allowed": False,
            },
        }
        if diagnostic_pavan_reaction_multipliers:
            propagated_charged_report["warnings"].append(
                "Charged HCN reaction-rate multipliers are hypothetical mechanism stresses, not fitted rates, measured uncertainty bounds, or validated predictions."
            )
    plasma_source_threshold_J = sum(
        float(cycle.get("threshold_energy_lower_bound_J", 0.0))
        for cycle in plasma_cycles
    )
    heavy_species_source_energy_J = sum(
        float(cycle.get("source_energy_used_J", 0.0)) for cycle in heavy_species_cycles
    )
    gas_energy_audits = [
        cycle["ideal_gas_internal_energy_audit"]
        for cycle in heavy_species_cycles
        if cycle.get("ideal_gas_internal_energy_audit") is not None
    ]
    gas_energy_audit_complete = bool(
        heavy_species_cycles
        and len(gas_energy_audits) == len(heavy_species_cycles)
        and all(item["all_nonzero_gas_species_covered"] for item in gas_energy_audits)
    )
    used_energy_J = declared_channel_energy_J + heavy_species_source_energy_J
    plasma_source_was_energy_limited = bool(
        plasma_cycles
        and plasma_source_threshold_J > deposited_energy_J * (1 + 1e-8)
    )
    if plasma_source_was_energy_limited and used_energy_J <= deposited_energy_J * (1 + 1e-8):
        energy_status = "accounted_sources_energy_limited_raw_opportunities_exceed_deposit"
    elif plasma_source_was_energy_limited:
        energy_status = "inconsistent_accounted_sources_exceed_deposit"
    elif plasma_cycles:
        energy_status = "conditional_source_opportunities_within_deposit"
    elif not electron_extent_total:
        energy_status = "no_electron_impact_channels"
    elif not energy_cost_complete:
        energy_status = "not_assessable_missing_channel_energy"
    elif any_energy_limited:
        energy_status = "energy_limited_to_declared_deposit"
    elif used_energy_J <= deposited_energy_J * (1 + 1e-8):
        energy_status = "within_declared_deposited_energy"
    else:
        energy_status = "inconsistent_channel_energy_exceeds_deposit"

    paired_diagnostics = []
    for species in request.species:
        if species.enantiomer_of is not None and species.id < species.enantiomer_of:
            left = final[index[species.id]]
            right = final[index[species.enantiomer_of]]
            paired_diagnostics.append(
                {
                    "species": species.id,
                    "enantiomer": species.enantiomer_of,
                    "difference_mol": float(left - right),
                    "racemic_within_numerical_tolerance": bool(
                        abs(left - right) <= max(left, right, 1.0) * 1e-9
                    ),
                    "note": "Equality is diagnostic only; it is not forced by the solver.",
                }
            )

    warnings = [
        "The apparatus layer contains no historical Miller-Urey product list and does not seed product discovery.",
        "Quantitative yield claims require sourced cross sections, thermochemistry, multiphase kinetics and experimental calibration.",
    ]
    if apparatus.operation_mode == "closed_stirred_flask":
        if finite_closed_contact:
            if coupled_closed_contact:
                warnings.append(
                    "Experimental coupled-contact path: targeted short-run runtime and material "
                    "checks pass, but full-duration physical validation is pending. "
                    "The closed stirred flask couples declared finite-rate Henry contact to "
                    "the stiff gas-chemistry solve. kLa is not inferred from HCN yield or "
                    "measured by this model. Aqueous acid-base updates still occur at substep "
                    "boundaries; timestep convergence, water-vapor replenishment, spatial "
                    "spark exposure, and thermal-energy closure remain unresolved."
                )
            else:
                warnings.append(
                    "The closed stirred flask uses a declared liquid-side kLa for finite-rate Henry "
                    "contact. Its value is not inferred from HCN yield or measured by this model. "
                    "Gas chemistry and water transfer are split across substeps, so quantitative use "
                    "requires timestep convergence. Water-vapor replenishment, spatial spark exposure, "
                    "and a complete thermal energy balance remain unresolved."
                )
        else:
            warnings.append(
                "The closed stirred flask uses instantaneous Henry partition at each substep as an "
                "upper-contact limit. It does not solve interfacial mass-transfer kinetics, water-vapor "
                "replenishment, spatial spark exposure, or a complete thermal energy balance."
            )
    if request.spark.mode == "boltzmann_eedf":
        warnings.append(
            "Published GRI/Pavan gas mechanisms already contain HCN species and reactions; "
            "unseeded HCN formation here is not de novo product-identity discovery."
        )
        warnings.append(
            "The gas-loop mixture drives a self-consistent EEDF. Energy-limited, explicitly resolved "
            "neutral fragments enter the primary Cantera gas mechanism. A separate propagated Pavan "
            "branch carries supported ions and excited states without mutating the primary inventory; "
            "unsupported fragment identities remain outside both inventory models."
        )
        if apparatus.discharge_zone_volume is None:
            warnings.append(
                "No discharge-zone volume was supplied, so the full gas loop is treated as active. "
                "This is an upper-bound residence-time assumption, not a measured spark geometry."
            )
        if apparatus.operation_mode == "continuous_overlap":
            warnings.append(
                "Continuous boil/discharge/condensation is approximated by sequential operator-split substeps. "
                "The selected substep count is a numerical control, not apparatus evidence; quantitative use "
                "requires time-step convergence."
            )
    else:
        warnings.append(
            "Electron-impact results are conditional on the selected discharge mode; a declared EEDF "
            "is not a self-consistent Boltzmann solution."
        )
    if apparatus.condensate_nonwater_return_fraction < 1:
        warnings.append(
            "The cold trap isolates a declared fraction of every non-water species in collected_condensate; this is a generic separation scenario, not a measured species-specific collection efficiency. Collected material cannot re-enter the gas loop or react in the active boiler sidecar."
        )
    if apparatus.condenser_phase_transfer_model == "flow_limited_capture":
        warnings.append(
            "Flow-limited condenser capture uses circulation turnover and a declared water contact hold-up, assuming each pass approaches Henry equilibrium. This is a mass-transfer upper bound, not a measured condenser efficiency."
        )
    elif apparatus.condenser_phase_transfer_model == "coupled_flow_transport":
        warnings.append(
            "A mean-field upper-bound Henry exchange couples boiler absorption/desorption and condenser capture to the stiff gas-chemistry integration. Spatial flow and interfacial resistance remain unresolved; the declared contact hold-up is hypothetical and timestep convergence must still be checked."
        )
    if apparatus.aqueous_acid_base_model == "hcn_nh3_water_equilibrium":
        warnings.append(
            "Explicit HCN/CN-, NH3/NH4+, CO2/HCO3-/CO3-- and water equilibrium is updated at apparatus substep boundaries. It omits other buffers, ionic activities and aqueous kinetics, and is not coupled inside the stiff gas/transfer solve; quantitative pH-dependent yield remains unvalidated."
        )
        warnings.append(
            "The electrical energy ledger bounds modeled discharge source terms; it does not include acid-base reaction enthalpy or close the full apparatus thermal energy balance."
        )
    effective_formaldehyde_pool_unresolved = (
        "CH2O" in index and request.spark.mode == "boltzmann_eedf"
        and apparatus.closed_flask_formaldehyde_hydration == "off"
    )
    if effective_formaldehyde_pool_unresolved:
        warnings.append(
            "The selected formaldehyde Henry value is an effective dissolved CH2O-plus-methanediol value. This apparatus currently stores that pool as CH2O without explicit hydration/dehydration; its gas-water timing and downstream HCN amount are not quantitatively validated."
        )
    hydration_heat_report = None
    if apparatus.closed_flask_formaldehyde_hydration != "off":
        hydrate_position = index[METHANEDIOL_AQ]
        net_hydration_extent_mol = float(
            final[hydrate_position] + sampled_total[hydrate_position]
            - initial[hydrate_position] - additions_total[hydrate_position]
        )
        hydration_enthalpy_J_mol = formaldehyde_hydration_enthalpy_J_mol()
        hydration_heat_report = {
            "status": "partial_isothermal_liquid_reaction_enthalpy_not_thermal_closure",
            "source_id": "winkelman_2002_formaldehyde_hydration_v1",
            "net_hydration_extent_mol": net_hydration_extent_mol,
            "reaction_enthalpy_J_mol": hydration_enthalpy_J_mol,
            "liquid_reaction_enthalpy_change_J": (
                net_hydration_extent_mol * hydration_enthalpy_J_mol
            ),
            "heat_release_if_isothermal_J": (
                -net_hydration_extent_mol * hydration_enthalpy_J_mol
            ),
            "not_included_in_temperature_feedback": True,
        }
        warnings.append(
            "Opt-in aqueous CH2O/methanediol hydration uses a cross-source derived free-CH2O "
            "Henry constant and the Winkelman 2002 aqueous rate correlation. The shared "
            "gas-water kLa is independently supplied, not measured here or fitted to HCN. "
            "pH is checked at substep boundaries, not "
            "solved continuously inside the stiff ODE. A partial hydration enthalpy term is "
            "reported without temperature feedback; apparatus heat flows and losses remain "
            "unclosed, so this is not a calibrated HCN yield."
        )
        if hydration_diagnostic_extrapolation:
            warnings.append(
                "Diagnostic hydration continues the 2002 temperature law below its measured "
                "pH-5 boundary, bounded at pH 4. This explicitly extrapolated sensitivity "
                "is not source-backed lower-pH kinetics and is not an experimental HCN yield."
            )
    if prepared["hypothetical_sources"]:
        warnings.append("Hypothetical transport or kinetic evidence is present; results are scenario output only.")
    if energy_status == "inconsistent_channel_energy_exceeds_deposit":
        warnings.append("Declared electron-impact reaction energy exceeds the electrical energy ledger.")
    if simulated_pressure_range_Pa and any(
        abs(value - pressure_Pa) > pressure_Pa * 0.1
        for value in simulated_pressure_range_Pa
    ):
        warnings.append(
            "The fixed-volume gas inventory does not match the declared reference pressure; "
            "kinetics use the inventory-derived pressure, and no hidden vent or pressure controller is applied."
        )
    result = {
        "status": "conditional_apparatus_simulation",
        "apparatus_id": request.id,
        "cycles_completed": apparatus.cycles,
        "total_duration_s": integration_steps * cycle_s,
        "species_identity": prepared["identities"],
        "reactions": prepared["reactions"],
        "cycles": cycles,
        "final_compartments_mol": {
            name: _compartment_dict(ids, values) for name, values in compartments.items()
        },
        "final_total_mol": _compartment_dict(ids, final),
        "reaction_extents_mol": extents_total,
        "material_ledger": {
            "initial_mol": _compartment_dict(ids, initial),
            "added_mol": _compartment_dict(ids, additions_total),
            "sampled_mol": _compartment_dict(ids, sampled_total),
            "final_mol": _compartment_dict(ids, final),
            "conserved_quantities": prepared["conserved"],
            "residual": {
                key: float(value) for key, value in zip(prepared["conserved"], residual, strict=True)
            },
            "tolerance_mol_by_conserved_quantity": {
                key: float(value) for key, value in zip(prepared["conserved"], tolerances, strict=True)
            },
            "maximum_absolute_residual_mol": float(np.max(np.abs(residual), initial=0)),
            "closed": closed,
        },
        "energy_ledger": {
            "scope": "electrical discharge source budget; not a full thermal/chemical enthalpy balance",
            "complete_thermal_energy_balance_closed": False,
            "formaldehyde_hydration_partial_heat": hydration_heat_report,
            "gas_stage_ideal_gas_internal_energy_change_J": (
                sum(float(item["change_J"]) for item in gas_energy_audits)
                if gas_energy_audit_complete else None
            ),
            "gas_stage_thermochemistry_coverage_complete": gas_energy_audit_complete,
            "gas_stage_energy_scope": (
                "GRI ideal-gas inventory change across neutral-discharge stages; "
                "may include coupled gas-water transfer; excludes liquid/ion enthalpy, "
                "wall heat, radiation, plasma internal energy and external heating. "
                "It is not a closed apparatus energy balance."
            ),
            "spark_mode": request.spark.mode,
            "duty_fraction": duty_fraction,
            "deposited_electrical_energy_J": electrical_estimate_J,
            "calorimetric_deposited_energy_J": (
                deposited_energy_J
                if deposited_energy_basis == "calorimetric_measurement"
                else None
            ),
            "deposited_energy_used_J": deposited_energy_J,
            "deposited_energy_basis": deposited_energy_basis,
            "declared_channel_energy_J": declared_channel_energy_J,
            "heavy_species_source_energy_J": heavy_species_source_energy_J,
            "total_accounted_channel_energy_J": used_energy_J,
            "plasma_source_threshold_energy_lower_bound_J": plasma_source_threshold_J,
            "status": energy_status,
            "reaction_was_energy_limited": any_energy_limited,
            "plasma_source_was_energy_limited": plasma_source_was_energy_limited,
        },
        "stereochemistry_diagnostics": paired_diagnostics,
        "aqueous_discovery": aqueous_discovery,
        "post_discovery_validation": post_discovery_validation,
        "historical_benchmark": historical_benchmark,
        "precursor_energy_benchmark": precursor_energy_benchmark,
        "numerical_controls": {
            "boil_pass_fraction": boil_pass,
            "discharge_pass_fraction": discharge_pass,
            "condense_pass_fraction": condense_pass,
            "gas_mass_action_exposure_fraction": gas_chemistry_fraction,
            "declared_electron_channel_exposure_fraction": electron_exposure_fraction,
            "heavy_gas_chemistry_residence_fraction": heavy_chemistry_residence_fraction,
            "discharge_zone_volume_L": discharge_zone_volume_L,
            "discharge_zone_fraction": discharge_zone_fraction,
            "effective_heavy_chemistry_residence_time_s_per_cycle": (
                declared_effective_heavy_chemistry_residence_s
            ),
            "effective_heavy_chemistry_residence_time_s_per_integration_step": (
                effective_heavy_chemistry_residence_s
            ),
            "operation_mode": apparatus.operation_mode,
            "integration_substeps_per_cycle": operation_substeps,
            "primary_source_forcing": apparatus.primary_source_forcing,
            "condenser_phase_transfer_model": apparatus.condenser_phase_transfer_model,
            "closed_flask_gas_water_transfer_model": apparatus.closed_flask_gas_water_transfer_model,
            "closed_flask_finite_contact_solver": apparatus.closed_flask_finite_contact_solver,
            "closed_flask_formaldehyde_hydration": apparatus.closed_flask_formaldehyde_hydration,
            "formaldehyde_hydration_minimum_ideal_pH": hydration_minimum_pH_ideal,
            "formaldehyde_hydration_pH_extrapolation_used": hydration_extrapolation_used,
            "closed_flask_liquid_side_kla_per_s": closed_flask_kla_per_s,
            "aqueous_acid_base_model": apparatus.aqueous_acid_base_model,
            "aqueous_acid_base_equilibrium_updates": acid_base_updates,
            "aqueous_acid_base_maximum_ledger_residual_mol": acid_base_maximum_residual_mol,
            "condenser_contact_water_holdup_L": contact_water_holdup_L,
            "operator_splitting_convergence_assessed_in_this_run": False,
            "discharge_zone_basis": (
                apparatus.discharge_zone_volume.basis
                if apparatus.discharge_zone_volume is not None
                else "unspecified_full_loop_upper_bound"
            ),
        },
        "evidence": {
            "hypothetical_sources": prepared["hypothetical_sources"],
            "all_rates_sourced": not prepared["hypothetical_sources"],
        },
        "plasma_model": {
            "mode": request.spark.mode,
            "eedf": (
                request.spark.eedf.model_dump(mode="json", exclude_none=True)
                if request.spark.eedf is not None
                else None
            ),
            "cross_section_channels": [
                {
                    "reaction_id": reaction["id"],
                    **(reaction["collision_rate"] or {"status": "inactive_spark_off"}),
                }
                for reaction in prepared["reactions"]
                if reaction["kind"] == "electron_collision"
            ],
            "self_consistent_eedf_solved": bool(plasma_cycles),
            "apparatus_mixture_coupled": request.spark.mode == "boltzmann_eedf",
            "cycle_summaries": plasma_cycles,
            "heavy_species_inventory_coupled": bool(heavy_species_cycles),
            "heavy_species_engine": "Cantera gri30.yaml" if heavy_species_cycles else None,
            "heavy_species_cycle_summaries": heavy_species_cycles,
            "phase_partition_cycle_summaries": partition_cycles,
            "propagated_charged_branch_coupled": propagated_charged_report is not None,
            "propagated_charged_branch": propagated_charged_report,
            "diagnostic_charged_hcn_final_rate_snapshots": (
                charged_hcn_final_rate_snapshots
                if diagnostic_charged_hcn_final_rate_each_substep else None
            ),
            "diagnostic_charged_hcn_integrated_snapshots": (
                charged_hcn_integrated_snapshots
                if diagnostic_charged_hcn_integrated_each_substep else None
            ),
        },
        "scientific_firewall": {
            "knowledge_policy": request.knowledge_policy,
            "mechanism_origin": request.mechanism_origin,
            "discovery_manifest_sha256": request.discovery_manifest_sha256,
            "known_products_accessible_to_discovery": bool(heavy_species_cycles),
            "reaction_recipe_used_for_discovery": bool(
                heavy_species_cycles
                or (request.reactions and request.mechanism_origin == "user_declared")
            ),
            "preexisting_gri_hcn_reactions_accessible_to_solver": bool(heavy_species_cycles),
            "preexisting_pavan_hcn_reactions_accessible_to_alternative_branch": (
                propagated_charged_report is not None
            ),
            "aqueous_reaction_templates_used_for_discovery": False,
            "blind_discovery_eligible": (
                request.mechanism_origin == "blind_discovery_artifact"
                and not heavy_species_cycles
            ),
            "quantitative_yield_claim_allowed": False,
            "effective_formaldehyde_pool_without_explicit_hydration": (
                effective_formaldehyde_pool_unresolved
            ),
            "formaldehyde_hydration_pH_extrapolation_used": hydration_extrapolation_used,
        },
        "operating_conditions": {
            "aqueous_acid_base_model": apparatus.aqueous_acid_base_model,
            "pressure_Pa": apparatus.pressure.as_unit("Pa"),
            "declared_reference_pressure_Pa": apparatus.pressure.as_unit("Pa"),
            "simulated_gas_pressure_range_Pa": simulated_pressure_range_Pa,
            "pressure_basis": "ideal-gas pressure from fixed gas volume and cycle inventory; declared pressure is a reference only",
            "boiler_temperature_K": apparatus.boiler_temperature.as_unit("K"),
            "discharge_temperature_K": apparatus.discharge_temperature.as_unit("K"),
            "condenser_temperature_K": apparatus.condenser_temperature.as_unit("K"),
            "boiler_liquid_volume_L": liquid_volume_L,
            "gas_loop_volume_L": gas_volume_L,
            "discharge_zone_volume_L": discharge_zone_volume_L,
            "discharge_zone_fraction": discharge_zone_fraction,
            "effective_heavy_chemistry_residence_time_s_per_cycle": (
                declared_effective_heavy_chemistry_residence_s
            ),
            "effective_heavy_chemistry_residence_time_s_per_integration_step": (
                effective_heavy_chemistry_residence_s
            ),
            "operation_mode": apparatus.operation_mode,
            "integration_substeps_per_cycle": operation_substeps,
            "primary_source_forcing": apparatus.primary_source_forcing,
            "operator_splitting_convergence_assessed_in_this_run": False,
            "discharge_zone_evidence": (
                apparatus.discharge_zone_volume.model_dump(mode="json", exclude_none=True)
                if apparatus.discharge_zone_volume is not None
                else None
            ),
            "condensate_trap_volume_L": trap_volume_L,
            "effective_aqueous_volume_L": effective_aqueous_volume_L,
            "effective_aqueous_volume_basis": "final H2O amount at 18.01528 mL/mol",
            "gas_circulation_rate_L_s": flow_L_s,
            "condensate_water_return_fraction": apparatus.condensate_return_fraction,
            "condensate_nonwater_return_fraction": apparatus.condensate_nonwater_return_fraction,
            "condenser_phase_transfer_model": apparatus.condenser_phase_transfer_model,
            "closed_flask_gas_water_transfer_model": apparatus.closed_flask_gas_water_transfer_model,
            "closed_flask_finite_contact_solver": apparatus.closed_flask_finite_contact_solver,
            "closed_flask_formaldehyde_hydration": apparatus.closed_flask_formaldehyde_hydration,
            "formaldehyde_hydration_minimum_ideal_pH": hydration_minimum_pH_ideal,
            "formaldehyde_hydration_pH_extrapolation_used": hydration_extrapolation_used,
            "closed_flask_liquid_side_kla_evidence": (
                apparatus.closed_flask_liquid_side_kla.model_dump(mode="json", exclude_none=True)
                if apparatus.closed_flask_liquid_side_kla is not None else None
            ),
            "condenser_contact_water_holdup_L": contact_water_holdup_L,
            "condenser_contact_water_holdup_evidence": (
                apparatus.condenser_contact_water_holdup.model_dump(mode="json", exclude_none=True)
                if apparatus.condenser_contact_water_holdup is not None else None
            ),
        },
        "warnings": warnings,
    }
    # HCN is inspected only after the inventory and blind aqueous discovery
    # have been completed; it cannot influence the simulated pathway search.
    from .hcn_formation import evaluate_aqueous_hcn_exposure, evaluate_hcn_formation
    from .aqueous_speciation import evaluate_aqueous_cyanide_speciation

    result["hcn_formation"] = evaluate_hcn_formation(result)
    from .hcn_ceiling import hcn_atom_ceiling
    result["hcn_formation"]["stoichiometric_atom_ceiling"] = hcn_atom_ceiling(request)
    result["hcn_formation"]["aqueous_time_exposure"] = evaluate_aqueous_hcn_exposure(
        hcn_integration_steps, result["operating_conditions"]
    )
    from .aqueous_hcn_hydrolysis import evaluate_aqueous_hcn_hydrolysis_screen
    result["hcn_formation"]["aqueous_hydrolysis_screen"] = (
        evaluate_aqueous_hcn_hydrolysis_screen(
            result["hcn_formation"]["aqueous_time_exposure"]
        )
    )
    if charged_hcn_integration_steps:
        result["hcn_formation"]["charged_branch_aqueous_time_exposure"] = (
            evaluate_aqueous_hcn_exposure(
                charged_hcn_integration_steps, result["operating_conditions"]
            )
        )
        result["hcn_formation"]["charged_branch_aqueous_hydrolysis_screen"] = (
            evaluate_aqueous_hcn_hydrolysis_screen(
                result["hcn_formation"]["charged_branch_aqueous_time_exposure"]
            )
        )
    result["aqueous_cyanide_speciation"] = evaluate_aqueous_cyanide_speciation(result)
    return result
