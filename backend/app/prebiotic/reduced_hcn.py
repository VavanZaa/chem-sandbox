"""Conditional, source-driven 0-D HCN gas/water model.

The effective fragment-production hazards are inputs, not inferred from an HCN
assay or a plasma calculation. This is deliberately separate from the Origins
EEDF model and does not claim quantitative validation of GRI at low T.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any, Literal

import cantera as ct
import numpy as np
from pydantic import Field, model_validator
from scipy.constants import R
from scipy.integrate import solve_ivp

from app.operations.schema import StrictModel
from app.thermochemistry.schema import EvidenceQuantity
from .aqueous_speciation import _hcn_pka, load_acid_base_subset
from .hcn_ceiling import hcn_atom_ceiling_from_atoms


class ScalarEvidence(StrictModel):
    value: float = Field(ge=0, le=14)
    source: str = Field(min_length=3, max_length=500)
    method: str = Field(min_length=3, max_length=200)
    basis: Literal["experimental", "literature", "hypothetical"]


class EffectiveFragmentSource(StrictModel):
    id: str = Field(min_length=1, max_length=80)
    parent: str = Field(min_length=1, max_length=64)
    products: dict[str, int] = Field(min_length=1, max_length=8)
    parent_conversion_hazard: EvidenceQuantity
    minimum_event_energy: EvidenceQuantity | None = None
    product_state_source: str = Field(min_length=3, max_length=500)
    product_state_basis: Literal["experimental", "literature", "hypothetical"]

    @model_validator(mode="after")
    def check_channel(self) -> "EffectiveFragmentSource":
        if any(count <= 0 or count > 8 for count in self.products.values()):
            raise ValueError("Effective-source product counts must be positive integers")
        if sum(self.products.values()) < 2:
            raise ValueError("Effective fragment sources require at least two product particles")
        if self.parent in self.products:
            raise ValueError("Effective-source parent cannot also be a product")
        if self.parent_conversion_hazard.as_unit("1/s") < 0:
            raise ValueError("Effective-source hazard cannot be negative")
        if self.minimum_event_energy is not None and self.minimum_event_energy.as_unit("J/mol") <= 0:
            raise ValueError("Minimum fragment event energy must be positive")
        return self


class ReducedHCNRequest(StrictModel):
    """Inputs are independently declared effective sources, not EEDF settings."""

    initial_gas_mol: dict[str, float] = Field(min_length=1, max_length=53)
    gas_volume: EvidenceQuantity
    water_volume: EvidenceQuantity
    temperature: EvidenceQuantity
    duration: EvidenceQuantity
    fixed_pH: ScalarEvidence
    hcn_henry_hcp: EvidenceQuantity
    hcn_liquid_side_kla: EvidenceQuantity
    deposited_energy: EvidenceQuantity | None = None
    sources: list[EffectiveFragmentSource] = Field(default_factory=list, max_length=40)
    sample_count: int = Field(default=24, ge=2, le=200)
    maximum_solver_step_s: float | None = Field(default=None, gt=0)
    hcn_assay_used_to_set_inputs: Literal[False]

    @model_validator(mode="after")
    def check_dimensions(self) -> "ReducedHCNRequest":
        if any(value < 0 for value in self.initial_gas_mol.values()) or sum(self.initial_gas_mol.values()) <= 0:
            raise ValueError("Initial gas amounts must be nonnegative with a positive total")
        if min(
            self.gas_volume.as_unit("m^3"), self.water_volume.as_unit("m^3"),
            self.temperature.as_unit("K"), self.duration.as_unit("s"),
            self.hcn_henry_hcp.as_unit("mol/(m^3*Pa)"),
        ) <= 0:
            raise ValueError("Volume, temperature, duration and Henry coefficient must be positive")
        if self.hcn_liquid_side_kla.as_unit("1/s") < 0:
            raise ValueError("kLa cannot be negative")
        if self.deposited_energy is not None and self.deposited_energy.as_unit("J") < 0:
            raise ValueError("Deposited energy cannot be negative")
        lower, upper = load_acid_base_subset()["applicability_temperature_K"]
        if not lower <= self.temperature.as_unit("K") <= upper:
            raise ValueError("HCN pKa temperature is outside its reviewed domain")
        ids = [item.id for item in self.sources]
        if len(ids) != len(set(ids)):
            raise ValueError("Effective-source ids must be unique")
        return self


def simulate_reduced_hcn(request: ReducedHCNRequest) -> dict[str, Any]:
    """Advance GRI chemistry plus measured/declared source and HCN contact.

    No numerical source strength is inferred from HCN. Source and gas/water
    parameters must be supplied even for a hypothetical sensitivity study.
    """
    gas = ct.Solution("gri30.yaml")
    mechanism_path = next(
        (Path(directory) / "gri30.yaml" for directory in ct.get_data_directories()
         if (Path(directory) / "gri30.yaml").is_file()),
        None,
    )
    mechanism_sha256 = (
        hashlib.sha256(mechanism_path.read_bytes()).hexdigest()
        if mechanism_path is not None else None
    )
    names = gas.species_names
    index = {name: i for i, name in enumerate(names)}
    unknown = set(request.initial_gas_mol) - set(index)
    if unknown:
        raise ValueError(f"Initial gas contains species absent from GRI: {sorted(unknown)}")
    initial = np.zeros(gas.n_species, dtype=float)
    for name, amount in request.initial_gas_mol.items():
        initial[index[name]] = amount
    channels: list[tuple[EffectiveFragmentSource, int, np.ndarray]] = []
    for channel in request.sources:
        unknown = {channel.parent, *channel.products} - set(index)
        if unknown:
            raise ValueError(f"Source {channel.id} contains species absent from GRI: {sorted(unknown)}")
        # A measured effective source may make elemental fragments, but may not
        # inject the target or a carbon-nitrogen-bonded precursor as a recipe.
        if any(
            gas.species(name).composition.get("C", 0) > 0
            and gas.species(name).composition.get("N", 0) > 0
            for name in channel.products
        ):
            raise ValueError(f"Source {channel.id} may not inject a C-N product")
        if any("(" in name or ")" in name for name in channel.products):
            raise ValueError(f"Source {channel.id} has a state-labelled product not resolved by this model")
        vector = np.zeros(gas.n_species, dtype=float)
        vector[index[channel.parent]] -= 1
        for name, count in channel.products.items():
            vector[index[name]] += count
        atom_residual = np.array([
            sum(gas.n_atoms(i, element) * vector[i] for i in range(gas.n_species))
            for element in gas.element_names
        ])
        if np.max(np.abs(atom_residual)) > 1e-12 or abs(float(gas.charges @ vector)) > 1e-12:
            raise ValueError(f"Source {channel.id} does not conserve atoms and charge")
        channels.append((channel, index[channel.parent], vector))

    temperature = request.temperature.as_unit("K")
    gas_volume = request.gas_volume.as_unit("m^3")
    water_volume = request.water_volume.as_unit("m^3")
    duration = request.duration.as_unit("s")
    henry = request.hcn_henry_hcp.as_unit("mol/(m^3*Pa)")
    kla = request.hcn_liquid_side_kla.as_unit("1/s")
    pka = _hcn_pka(temperature, load_acid_base_subset())
    neutral_fraction = 1 / (1 + 10 ** (request.fixed_pH.value - pka))
    hazards = np.array([item.parent_conversion_hazard.as_unit("1/s") for item, _, _ in channels])
    hcn_position = index["HCN"]
    y0 = np.concatenate([initial, np.zeros(1 + len(channels))])

    def rhs(_time: float, state: np.ndarray) -> np.ndarray:
        gas_amount = np.maximum(state[:gas.n_species], 0.0)
        total = float(gas_amount.sum())
        if total <= 0:
            raise ValueError("Gas inventory vanished during integration")
        pressure = total * R * temperature / gas_volume
        gas.TPX = temperature, pressure, gas_amount / total
        gas_change = np.asarray(gas.net_production_rates, dtype=float) * gas_volume * 1000.0
        source_rates = np.array([
            hazards[i] * gas_amount[parent]
            for i, (_, parent, _) in enumerate(channels)
        ])
        for rate, (_, _, vector) in zip(source_rates, channels, strict=True):
            gas_change += rate * vector
        aqueous = max(float(state[gas.n_species]), 0.0)
        hcn_partial_pressure = gas_amount[hcn_position] * R * temperature / gas_volume
        transfer_rate = kla * water_volume * (
            henry * hcn_partial_pressure - neutral_fraction * aqueous / water_volume
        )
        gas_change[hcn_position] -= transfer_rate
        return np.concatenate([gas_change, [transfer_rate], source_rates])

    sampled = np.linspace(0.0, duration, request.sample_count)
    integration = solve_ivp(
        rhs, (0.0, duration), y0, method="BDF", t_eval=sampled,
        rtol=1e-7, atol=1e-18,
        max_step=(request.maximum_solver_step_s or math.inf),
    )
    if not integration.success or not np.isfinite(integration.y).all():
        raise ValueError(f"Reduced HCN integration failed: {integration.message}")
    final_state = integration.y[:, -1]
    tolerance = max(float(initial.sum()) * 1e-10, 1e-15)
    if float(np.min(integration.y)) < -tolerance:
        raise ValueError("Reduced HCN integration produced a material negative state")
    final_gas = np.maximum(final_state[:gas.n_species], 0.0)
    aqueous = max(float(final_state[gas.n_species]), 0.0)
    extents = np.maximum(final_state[gas.n_species + 1:], 0.0)
    initial_atoms = {
        element: sum(initial[i] * gas.n_atoms(i, element) for i in range(gas.n_species))
        for element in gas.element_names
    }
    final_atoms = {
        element: sum(final_gas[i] * gas.n_atoms(i, element) for i in range(gas.n_species))
        + (aqueous if element in {"C", "H", "N"} else 0.0)
        for element in gas.element_names
    }
    residual = {element: final_atoms[element] - initial_atoms[element] for element in gas.element_names}
    if max(abs(value) for value in residual.values()) > tolerance:
        raise ValueError(f"Reduced HCN atom ledger did not close: {residual}")
    minimum_source_energy = None
    energy_complete = all(item.minimum_event_energy is not None for item, _, _ in channels)
    if energy_complete:
        minimum_source_energy = sum(
            float(extent) * item.minimum_event_energy.as_unit("J/mol")
            for extent, (item, _, _) in zip(extents, channels, strict=True)
        )
    deposited = request.deposited_energy.as_unit("J") if request.deposited_energy else None
    if deposited is not None and minimum_source_energy is not None and minimum_source_energy > deposited * (1 + 1e-8):
        raise ValueError("Minimum declared fragment energy exceeds deposited electrical energy")
    initial_hcn = float(initial[hcn_position])
    final_hcn = float(final_gas[hcn_position]) + aqueous
    ceiling = hcn_atom_ceiling_from_atoms(
        initial_atoms, initial_hcn,
        atom_inventory_scope="initial_reactive_gas_only_water_is_inert_solvent_in_this_model",
    )
    input_evidence = {
        "gas_volume": request.gas_volume.basis,
        "water_volume": request.water_volume.basis,
        "temperature": request.temperature.basis,
        "duration": request.duration.basis,
        "fixed_pH": request.fixed_pH.basis,
        "hcn_henry_hcp": request.hcn_henry_hcp.basis,
        "hcn_liquid_side_kla": request.hcn_liquid_side_kla.basis,
        **({"deposited_energy": request.deposited_energy.basis} if request.deposited_energy else {}),
    }
    for item, _, _ in channels:
        input_evidence[f"source:{item.id}:hazard"] = item.parent_conversion_hazard.basis
        input_evidence[f"source:{item.id}:product_state"] = item.product_state_basis
        if item.minimum_event_energy is not None:
            input_evidence[f"source:{item.id}:minimum_energy"] = item.minimum_event_energy.basis
    hypothetical_inputs = sorted(key for key, basis in input_evidence.items() if basis == "hypothetical")
    active_source_count = int(np.count_nonzero(hazards))
    return {
        "status": "conditional_reduced_mechanistic_run",
        "engine": "Cantera GRI-3.0 rates + SciPy BDF + finite-rate HCN Henry contact",
        "cantera_version": ct.__version__,
        "gas_mechanism_sha256": mechanism_sha256,
        "hcn_assay_used_to_set_inputs": False,
        "experimental_yield_validated": False,
        "active_effective_source_count": active_source_count,
        "hypothetical_input_fields": hypothetical_inputs,
        "initial_hcn_mol": initial_hcn,
        "final_gas_hcn_mol": float(final_gas[hcn_position]),
        "final_aqueous_hcn_plus_cn_mol": aqueous,
        "final_total_hcn_plus_cn_mol": final_hcn,
        "net_hcn_plus_cn_formed_mol": final_hcn - initial_hcn,
        "aqueous_neutral_hcn_fraction": neutral_fraction,
        "fixed_pH": request.fixed_pH.value,
        "hcn_pKa": pka,
        "atom_ceiling": ceiling,
        "atom_residual_mol": residual,
        "atom_ledger_closed": True,
        "aqueous_charge_speciation_explicitly_balanced": False,
        "source_extents_mol": {item.id: float(extent) for extent, (item, _, _) in zip(extents, channels, strict=True)},
        "minimum_fragment_energy_J": minimum_source_energy,
        "deposited_energy_J": deposited,
        "minimum_fragment_energy_within_deposit": (
            minimum_source_energy <= deposited if minimum_source_energy is not None and deposited is not None else None
        ),
        "full_thermal_energy_balance_closed": False,
        "numerical_convergence_assessed_in_this_run": False,
        "solver": {"method": "SciPy BDF", "rhs_evaluations": integration.nfev,
                   "maximum_step_s": request.maximum_solver_step_s},
        "conditional_inputs": request.model_dump(mode="json"),
        "samples": [
            {"time_s": float(time), "gas_hcn_mol": max(float(integration.y[hcn_position, i]), 0.0),
             "aqueous_hcn_plus_cn_mol": max(float(integration.y[gas.n_species, i]), 0.0)}
            for i, time in enumerate(integration.t)
        ],
        "limitations": [
            *(["All effective fragment hazards are zero; this is a no-source control, not a sparked reaction."]
              if active_source_count == 0 else []),
            "Effective fragment hazards, not spark geometry or EEDF, are supplied externally; no default source is invented.",
            "GRI-3.0 contains HCN reactions and is not independently validated for this low-temperature wet plasma.",
            "Only HCN transfers to water; pH and temperature are fixed; water-phase reactions and full heat flow are omitted.",
            "Aqueous cyanide is one HCN-equivalent pool: CN- and its balancing proton are not separate dynamic species.",
            "This conditional result is not a measured yield or blind HCN discovery.",
            "Repeat a proposed case at tighter solver maximum steps; the single-run result does not prove numerical convergence.",
        ],
    }


def audit_reduced_hcn(request: ReducedHCNRequest) -> dict[str, Any]:
    """Check conditional numerical stability and conservation controls.

    This audit deliberately has no assay input and cannot establish whether an
    effective fragment source represents a real spark.
    """
    duration = request.duration.as_unit("s")
    coarse_step = min(request.maximum_solver_step_s or duration / 16, duration / 16)

    def with_changes(**changes: Any) -> ReducedHCNRequest:
        return ReducedHCNRequest.model_validate({**request.model_dump(mode="json"), **changes})

    coarse = simulate_reduced_hcn(with_changes(maximum_solver_step_s=coarse_step))
    refined_step = coarse_step / 2
    refined = simulate_reduced_hcn(with_changes(maximum_solver_step_s=refined_step))
    absolute_floor_mol = 1e-12
    relative_tolerance = 0.01
    metrics = {
        "final_gas_hcn_mol": (coarse["final_gas_hcn_mol"], refined["final_gas_hcn_mol"]),
        "final_aqueous_hcn_plus_cn_mol": (
            coarse["final_aqueous_hcn_plus_cn_mol"], refined["final_aqueous_hcn_plus_cn_mol"]
        ),
        "final_total_hcn_plus_cn_mol": (
            coarse["final_total_hcn_plus_cn_mol"], refined["final_total_hcn_plus_cn_mol"]
        ),
        **{
            f"source_extent:{source_id}": (coarse["source_extents_mol"][source_id], extent)
            for source_id, extent in refined["source_extents_mol"].items()
        },
    }
    comparisons = {
        name: {
            "coarse_mol": before,
            "refined_mol": after,
            "absolute_difference_mol": abs(after - before),
            "allowed_difference_mol": absolute_floor_mol + relative_tolerance * max(abs(before), abs(after)),
            "within_declared_numerical_tolerance": abs(after - before)
            <= absolute_floor_mol + relative_tolerance * max(abs(before), abs(after)),
        }
        for name, (before, after) in metrics.items()
    }
    solver_step_consistent = all(
        item["within_declared_numerical_tolerance"] for item in comparisons.values()
    )

    source_off = simulate_reduced_hcn(with_changes(sources=[], maximum_solver_step_s=refined_step))
    no_source_background = source_off["final_total_hcn_plus_cn_mol"]
    gas = ct.Solution("gri30.yaml")
    controls: dict[str, dict[str, Any]] = {}
    for element in ("C", "N"):
        remaining = {
            name: amount for name, amount in request.initial_gas_mol.items()
            if gas.species(name).composition.get(element, 0) == 0
        }
        if not remaining or sum(remaining.values()) <= 0:
            controls[f"no_{element.lower()}_feed"] = {
                "status": "not_runnable_without_a_remaining_gas_inventory",
                "final_total_hcn_plus_cn_mol": None,
            }
            continue
        result = simulate_reduced_hcn(with_changes(
            initial_gas_mol=remaining, maximum_solver_step_s=refined_step
        ))
        amount = result["final_total_hcn_plus_cn_mol"]
        controls[f"no_{element.lower()}_feed"] = {
            "status": "below_reporting_floor" if amount <= absolute_floor_mol else "control_failed",
            "final_total_hcn_plus_cn_mol": amount,
            "atom_ledger_closed": result["atom_ledger_closed"],
        }
    controls_passed = all(
        item["status"] != "control_failed" for item in controls.values()
    )
    controls_complete = all(
        item["status"] != "not_runnable_without_a_remaining_gas_inventory"
        for item in controls.values()
    )
    audit_outcome = (
        "incomplete" if not controls_complete else
        "checks_passed" if solver_step_consistent and controls_passed else "check_failed"
    )
    return {
        "status": "conditional_numerical_audit",
        "audit_outcome": audit_outcome,
        "experimental_yield_validated": False,
        "hcn_assay_used_to_set_inputs": False,
        "numerical_checks_passed": audit_outcome == "checks_passed",
        "solver_step_consistent": solver_step_consistent,
        "negative_controls_passed": controls_passed,
        "negative_controls_complete": controls_complete,
        "coarse_maximum_solver_step_s": coarse_step,
        "refined_maximum_solver_step_s": refined_step,
        "relative_tolerance": relative_tolerance,
        "absolute_reporting_floor_mol": absolute_floor_mol,
        "step_comparisons": comparisons,
        "source_off_control": {
            "final_total_hcn_plus_cn_mol": no_source_background,
            "below_reporting_floor": no_source_background <= absolute_floor_mol,
            "interpretation": "Background chemistry or seeded HCN may remain without effective sources; this is not automatically a failure.",
        },
        "atom_removal_controls": controls,
        "refined_result": refined,
        "interpretation": (
            "These are conditional solver and conservation checks only. Agreement cannot validate "
            "fragment hazards, product states, gas rates, kLa, pH, heat flow, or a real HCN yield."
        ),
    }
