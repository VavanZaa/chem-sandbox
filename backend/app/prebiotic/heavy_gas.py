"""Exploratory heavy-species chemistry for coupled plasma cycles.

GRI-Mech is used as a broad, independently published H/C/N/O gas mechanism.
It is not a Miller-Urey mechanism and is not tuned against the historical
product list. The inventory-driving branch still uses a declared neutral-
fragment closure. A separate non-mutating Pavan shadow retains supported ions,
electrons, and two-temperature recombination rates with explicit charge checks.
Both layers remain non-quantitative.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import math
from pathlib import Path
from threading import local
from typing import Any, Callable

import cantera as ct
import numpy as np
from scipy.constants import Avogadro, physical_constants
from scipy.integrate import solve_ivp

from .formaldehyde_hydration import (
    METHANEDIOL_AQ,
    formaldehyde_contact_hydration_rates_mol_s,
)
from .hnc_thermo import nist_hnc_room_temperature_species


FARADAY_J_PER_EV_MOL = physical_constants["Faraday constant"][0]
MECHANISM = "gri30.yaml"
# Internal mechanism-stress controls only. Indices are Cantera's zero-based
# GRI-3.0 order; equation guards below prevent a silent mechanism-version swap.
GRI_HCN_DIAGNOSTIC_REACTIONS = {
    236: "H + HCN (+M) <=> H2CN (+M)",
    237: "H2CN + N <=> CH2 + N2",
    275: "CH3 + N <=> H2 + HCN",
}
# Pearce, Ayers & Pudritz, J. Phys. Chem. A 123, 1861-1873 (2019),
# DOI 10.1021/acs.jpca.8b11323, Table 4 / theoretical case study 1.
# Singlet-surface, barrierless forward rate; their constant fit is 50-400 K.
LOW_T_H2CN_H_TO_HCN_RATE_CM3_MOLECULE_S = 1.8e-11
LOW_T_H2CN_H_TO_HCN_EQUATION = "H + H2CN <=> H2 + HCN"
# Hébrard et al., A&A 541, A21 (2012), DOI 10.1051/0004-6361/201218837,
# Appendix A.5, room-temperature recommended N(4S)+H2CN exit channels.
# The split is estimated from an experimentally constrained overall rate;
# it is not an experimental measurement of either branch in this flask.
HEBRARD_2012_H2CN_N_TO_CH2_N2_RATE_CM3_MOLECULE_S = 3.0e-11
HEBRARD_2012_H2CN_N_TO_HCN_NH_RATE_CM3_MOLECULE_S = 1.0e-11
HEBRARD_2012_H2CN_N_TO_HCN_NH_EQUATION = "H2CN + N <=> HCN + NH"
# Hébrard et al. (2012), Appendix A.4: room-temperature estimated H2CN+H
# exits. Their sum (7.2e-11) is compatible with the experimental total-loss
# lower bound (>7e-11), but their HCN/HNC split is not measured.
HEBRARD_2012_H2CN_H_TO_HCN_H2_RATE_CM3_MOLECULE_S = 6.0e-11
HEBRARD_2012_H2CN_H_TO_HNC_H2_RATE_CM3_MOLECULE_S = 1.2e-11
HEBRARD_2012_H2CN_H_TO_HNC_H2_EQUATION = "H + H2CN <=> H2 + HNC"
# Hébrard et al. (2012), Appendix A.2, estimated relaxed-HNC isomerization
# from the Sumathi & Nguyen potential-energy/rate study with a barrier average.
# k = 4.03e-11 exp(-1200/T) cm3 molecule^-1 s^-1; their stated uncertainty
# factor at 300 K is 10. This is never part of the default GRI phase.
HEBRARD_2012_HNC_H_TO_HCN_H_PREEXP_CM3_MOLECULE_S = 4.03e-11
HEBRARD_2012_HNC_H_TO_HCN_H_BARRIER_K = 1200.0
HEBRARD_2012_HNC_H_TO_HCN_H_EQUATION = "H + HNC <=> H + HCN"
PAVAN_MECHANISM = "example_data/methane-plasma-pavan-2023.yaml"
PAVAN_PHASE = "gri30-plasma"
PAVAN_HCN_DIAGNOSTIC_REACTIONS = {
    126: "CH + H2O <=> CH2O + H",
    236: "H + HCN (+M) <=> H2CN (+M)",
    240: "CH + N2 (+M) <=> HCNN (+M)",
    257: "HCNN + O <=> HCN + NO",
    260: "H + HCNN <=> CH2 + N2",
}
KAMARINOPOULOU_MECHANISM = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "reference"
    / "kamarinopoulou_2024"
    / "kamarinopoulou_2024_neutral.yaml"
)
KAMARINOPOULOU_MECHANISM_SHA256 = (
    "33c2c3044248972da20ebff6f7db90f97b8bab26c4f4089f275c42add883172c"
)
KAMARINOPOULOU_PUBLISHED_RESIDENCE_TIME_S = 0.018
_PAVAN_EXCITED_N2 = {
    **{f"N2(v{level})": f"N2(v{level})" for level in range(1, 9)},
    "N2(A3,v0-4)": "N2(A)",
    "N2(A3,v5-9)": "N2(A)",
    "N2(A3,v10-)": "N2(A)",
    "N2(B3)": "N2(B)",
    "N2(C3)": "N2(C)",
    "N2(a'1)": "N2(a')",
}

_thread_phases = local()


def _bounded_rate_jacobian(
    rhs: Callable[[float, np.ndarray], np.ndarray],
    time_s: float,
    amounts: np.ndarray,
    dynamic_columns: list[int],
) -> np.ndarray:
    """Difference active rate coordinates; passive/inactive columns are exactly zero.

    BDF's adaptive finite-difference multiplier can overflow when a state
    coordinate cannot affect any rate. A fixed representable mole perturbation
    avoids that failure without relaxing the ODE tolerances or balance checks.
    """
    base = rhs(time_s, amounts)
    jacobian = np.zeros((amounts.size, amounts.size), dtype=float)
    for position in dynamic_columns:
        perturbed = amounts.copy()
        step = max(
            math.sqrt(np.finfo(float).eps) * abs(float(amounts[position])),
            1e-22,
        )
        perturbed[position] += step
        actual_step = float(perturbed[position] - amounts[position])
        jacobian[:, position] = (rhs(time_s, perturbed) - base) / actual_step
    return jacobian


def _cached_phase(key: str, factory: Any) -> Any:
    """Reuse mutable Cantera mechanism objects only within one worker thread."""
    cache = getattr(_thread_phases, "cache", None)
    if cache is None:
        cache = {}
        _thread_phases.cache = cache
    if key not in cache:
        cache[key] = factory()
    return cache[key]


def _gri_phase_with_diagnostic_rates(
    multipliers: tuple[tuple[int, float], ...],
    *,
    include_literature_h2cn_h_to_hcn: bool = False,
    include_hebrard_h2cn_n_branches: bool = False,
    include_hebrard_h2cn_h_isomer_branches: bool = False,
    include_hebrard_hnc_h_isomerization: bool = False,
) -> ct.Solution:
    """Isolate an opt-in rate stress from the cached unmodified GRI phase."""
    if include_literature_h2cn_h_to_hcn and include_hebrard_h2cn_h_isomer_branches:
        raise ValueError("Choose either Pearce or Hebrard H2CN + H diagnostic branches")
    if include_hebrard_hnc_h_isomerization and not include_hebrard_h2cn_h_isomer_branches:
        raise ValueError("HNC + H isomerization requires the opt-in HNC species branch")
    if (
        not multipliers and not include_literature_h2cn_h_to_hcn
        and not include_hebrard_h2cn_n_branches
        and not include_hebrard_h2cn_h_isomer_branches
        and not include_hebrard_hnc_h_isomerization
    ):
        return _cached_phase("gri", lambda: ct.Solution(MECHANISM))
    if len({index for index, _ in multipliers}) != len(multipliers):
        raise ValueError("Diagnostic GRI reaction indices must be unique")
    normalized = tuple(sorted((int(index), float(value)) for index, value in multipliers))
    for index, value in normalized:
        if index not in GRI_HCN_DIAGNOSTIC_REACTIONS:
            raise ValueError(f"Unsupported diagnostic GRI reaction index: {index}")
        if not math.isfinite(value) or value < 0:
            raise ValueError("Diagnostic GRI rate multipliers must be finite and nonnegative")

    def factory() -> ct.Solution:
        phase = ct.Solution(MECHANISM)
        if include_hebrard_h2cn_h_isomer_branches:
            phase.add_species(nist_hnc_room_temperature_species(298.15))
        for index, value in normalized:
            if phase.reaction_equations()[index] != GRI_HCN_DIAGNOSTIC_REACTIONS[index]:
                raise ValueError("Diagnostic GRI reaction index no longer matches its equation")
            phase.set_multiplier(value, index)
        if include_literature_h2cn_h_to_hcn:
            if LOW_T_H2CN_H_TO_HCN_EQUATION in phase.reaction_equations():
                raise ValueError("The literature H2CN + H reaction is already in the mechanism")
            phase.add_reaction(ct.Reaction(
                equation=LOW_T_H2CN_H_TO_HCN_EQUATION,
                rate=ct.ArrheniusRate(
                    LOW_T_H2CN_H_TO_HCN_RATE_CM3_MOLECULE_S * Avogadro / 1000.0,
                    0.0,
                    0.0,
                ),
            ))
        if include_hebrard_h2cn_n_branches:
            if phase.reaction_equations()[237] != GRI_HCN_DIAGNOSTIC_REACTIONS[237]:
                raise ValueError("GRI H2CN + N reaction index no longer matches")
            if HEBRARD_2012_H2CN_N_TO_HCN_NH_EQUATION in phase.reaction_equations():
                raise ValueError("The literature H2CN + N branch is already present")
            competing = phase.reaction(237)
            competing.rate = ct.ArrheniusRate(
                HEBRARD_2012_H2CN_N_TO_CH2_N2_RATE_CM3_MOLECULE_S
                * Avogadro / 1000.0,
                0.0,
                0.0,
            )
            phase.modify_reaction(237, competing)
            phase.add_reaction(ct.Reaction(
                equation=HEBRARD_2012_H2CN_N_TO_HCN_NH_EQUATION,
                rate=ct.ArrheniusRate(
                    HEBRARD_2012_H2CN_N_TO_HCN_NH_RATE_CM3_MOLECULE_S
                    * Avogadro / 1000.0,
                    0.0,
                    0.0,
                ),
            ))
        if include_hebrard_h2cn_h_isomer_branches:
            phase.add_reaction(ct.Reaction(
                equation=LOW_T_H2CN_H_TO_HCN_EQUATION,
                rate=ct.ArrheniusRate(
                    HEBRARD_2012_H2CN_H_TO_HCN_H2_RATE_CM3_MOLECULE_S
                    * Avogadro / 1000.0, 0.0, 0.0,
                ),
            ))
            phase.add_reaction(ct.Reaction(
                equation=HEBRARD_2012_H2CN_H_TO_HNC_H2_EQUATION,
                rate=ct.ArrheniusRate(
                    HEBRARD_2012_H2CN_H_TO_HNC_H2_RATE_CM3_MOLECULE_S
                    * Avogadro / 1000.0, 0.0, 0.0,
                ),
            ))
        if include_hebrard_hnc_h_isomerization:
            phase.add_reaction(ct.Reaction(
                equation=HEBRARD_2012_HNC_H_TO_HCN_H_EQUATION,
                rate=ct.ArrheniusRate(
                    HEBRARD_2012_HNC_H_TO_HCN_H_PREEXP_CM3_MOLECULE_S
                    * Avogadro / 1000.0,
                    0.0,
                    HEBRARD_2012_HNC_H_TO_HCN_H_BARRIER_K * ct.gas_constant,
                ),
            ))
        return phase

    return _cached_phase(
        f"gri-diagnostic:{normalized!r}:h2cn-h={include_literature_h2cn_h_to_hcn}"
        f":h2cn-n={include_hebrard_h2cn_n_branches}"
        f":h2cn-h-isomers={include_hebrard_h2cn_h_isomer_branches}"
        f":hnc-h-isomerization={include_hebrard_hnc_h_isomerization}",
        factory,
    )


def _formula(composition: dict[str, int]) -> str:
    return "".join(
        f"{element}{count if count != 1 else ''}"
        for element, count in composition.items()
    )


def _cantera_error_summary(error: Exception) -> str:
    lines = [
        line.strip()
        for line in str(error).splitlines()
        if line.strip() and set(line.strip()) != {"*"}
    ]
    preferred = [
        line
        for line in lines
        if line.startswith(("CVodes error", "CanteraError thrown", "ydot contains"))
    ]
    return " | ".join((preferred or lines)[:2])[:500]


def gri_ideal_gas_internal_energy_inventory(
    amounts_mol: np.ndarray, prepared: dict[str, Any], temperature_K: float,
    *, diagnostic_phase: ct.Solution | None = None,
) -> dict[str, Any]:
    """Inventory gas internal energy from GRI (and opt-in HNC) thermo.

    Cantera reports partial molar internal energies in J/kmol. This audit
    multiplies each by its actual amount in kmol; it does not silently assign
    gas-phase thermochemistry to ions or dissolved/condensed material.
    Absolute values depend on the GRI formation-energy reference, while the
    difference between two inventories on that same reference is useful.
    """
    if not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("Gas internal-energy temperature must be positive")
    amounts = np.asarray(amounts_mol, dtype=float)
    if amounts.shape != (len(prepared["ids"]),) or not np.all(np.isfinite(amounts)):
        raise ValueError("Gas internal-energy inventory has invalid amounts")
    if np.min(amounts, initial=0.0) < -1e-14:
        raise ValueError("Gas internal-energy inventory contains negative amounts")

    phase = diagnostic_phase or _cached_phase("gri", lambda: ct.Solution(MECHANISM))
    phase.TP = temperature_K, ct.one_atm
    molar_internal_energy_J_per_mol = np.asarray(
        phase.partial_molar_int_energies, dtype=float
    ) / 1000.0
    covered = np.asarray(
        [max(float(amounts[prepared["index"][name]]), 0.0) for name in phase.species_names],
        dtype=float,
    )
    supported = set(phase.species_names)
    unsupported = {
        name: float(amounts[position])
        for position, name in enumerate(prepared["ids"])
        if name not in supported and amounts[position] > 0
    }
    return {
        "basis": (
            "GRI-Mech plus NIST room-temperature HNC ideal-gas internal energies"
            if "HNC" in phase.species_names else
            "GRI-Mech ideal-gas species partial molar internal energies"
        ),
        "temperature_K": temperature_K,
        "covered_internal_energy_J": float(
            np.dot(covered, molar_internal_energy_J_per_mol)
        ),
        "covered_amount_mol": float(np.sum(covered)),
        "unsupported_gas_species_mol": unsupported,
        "all_nonzero_species_covered": not unsupported,
        "excluded_phases": ["liquid", "aqueous_ions", "walls", "radiation", "plasma_excitation"],
    }


def gri_ideal_gas_enthalpy_inventory(
    amounts_mol: np.ndarray, prepared: dict[str, Any], temperature_K: float,
    *, diagnostic_phase: ct.Solution | None = None,
) -> dict[str, Any]:
    """Ideal-gas stream enthalpy on the same GRI formation-energy reference.

    This is a gas-only steady-flow diagnostic, not a plasma or wall heat
    balance. Cantera partial molar enthalpies are J/kmol and the inventory
    amounts are mol, hence the factor of 1/1000.
    """
    if not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("Gas enthalpy temperature must be positive")
    amounts = np.asarray(amounts_mol, dtype=float)
    if amounts.shape != (len(prepared["ids"]),) or not np.all(np.isfinite(amounts)):
        raise ValueError("Gas enthalpy inventory has invalid amounts")
    if np.min(amounts, initial=0.0) < -1e-14:
        raise ValueError("Gas enthalpy inventory contains negative amounts")

    phase = diagnostic_phase or _cached_phase("gri", lambda: ct.Solution(MECHANISM))
    phase.TP = temperature_K, ct.one_atm
    molar_enthalpy_J_per_mol = np.asarray(phase.partial_molar_enthalpies, dtype=float) / 1000.0
    covered = np.asarray(
        [max(float(amounts[prepared["index"][name]]), 0.0) for name in phase.species_names],
        dtype=float,
    )
    supported = set(phase.species_names)
    unsupported = {
        name: float(amounts[position])
        for position, name in enumerate(prepared["ids"])
        if name not in supported and amounts[position] > 0
    }
    return {
        "basis": "GRI-Mech ideal-gas species partial molar enthalpies",
        "temperature_K": temperature_K,
        "covered_enthalpy_J": float(np.dot(covered, molar_enthalpy_J_per_mol)),
        "covered_amount_mol": float(np.sum(covered)),
        "unsupported_gas_species_mol": unsupported,
        "all_nonzero_species_covered": not unsupported,
        "excluded_phases": ["liquid", "aqueous_ions", "walls", "radiation", "plasma_excitation"],
    }


def augment_prepared_with_gri(
    prepared: dict[str, Any], *, diagnostic_hnc_temperature_K: float | None = None,
) -> None:
    """Append zero-inventory GRI species so Cantera products enter the ledger."""

    phase = _cached_phase("gri", lambda: ct.Solution(MECHANISM))
    original_count = len(prepared["ids"])
    species_to_add = phase.species()
    if diagnostic_hnc_temperature_K is not None:
        species_to_add.append(nist_hnc_room_temperature_species(diagnostic_hnc_temperature_K))
    for species in species_to_add:
        composition = {
            element: int(round(count))
            for element, count in species.composition.items()
            if count
        }
        charge = int(round(species.charge))
        if species.name in prepared["index"]:
            identity = prepared["identities"][prepared["index"][species.name]]
            if identity["composition"] != composition or identity["charge"] != charge:
                raise ValueError(
                    f"Declared species {species.name} disagrees with {MECHANISM} composition"
                )
            continue
        prepared["ids"].append(species.name)
        prepared["identities"].append(
            {
                "id": species.name,
                "name": species.name,
                "composition": composition,
                "formula": _formula(composition),
                "charge": charge,
                "isomeric_smiles": None,
                "stereo_status": "mechanism_formula_only",
                "stereochemistry": [],
                "phase_behavior": "permanent_gas",
                "enantiomer_of": None,
                "identity_note": (
                    "Opt-in NIST room-temperature HNC diagnostic; no reviewed water fate"
                    if species.name == "HNC" else
                    f"Automatically imported from Cantera {MECHANISM}"
                ),
            }
        )
        prepared["transfer"].append(
            {
                "rates_s_inv": {"evaporation": 0.0, "condensation": 0.0, "wall_loss": 0.0},
                "evidence": {"evaporation": None, "condensation": None, "wall_loss": None},
            }
        )

    prepared["index"] = {key: position for position, key in enumerate(prepared["ids"])}
    added = len(prepared["ids"]) - original_count
    if added:
        prepared["stoichiometry"] = np.vstack(
            (
                prepared["stoichiometry"],
                np.zeros((added, prepared["stoichiometry"].shape[1]), dtype=float),
            )
        )
    inventories = []
    for identity in prepared["identities"]:
        inventory = Counter(identity["composition"])
        inventory["charge"] = identity["charge"]
        inventories.append(inventory)
    prepared["conserved"] = sorted({key for item in inventories for key in item})
    prepared["atomic"] = np.asarray(
        [
            [inventory.get(key, 0) for inventory in inventories]
            for key in prepared["conserved"]
        ],
        dtype=float,
    )


def _source_energy_annotations_match_threshold(source: dict[str, Any]) -> bool:
    """Reject product-label eV notes that disagree with the collision threshold."""

    notes = [
        float(product["energy_annotation_eV"])
        for product in source.get("products", [])
        if product.get("energy_annotation_eV") is not None
    ]
    if not notes:
        return True
    threshold = source.get("threshold_eV")
    if threshold is None:
        return False
    return all(
        math.isclose(note, float(threshold), rel_tol=1e-4, abs_tol=1e-4)
        for note in notes
    )


def _neutral_fragment_extents(
    gas_amounts_mol: np.ndarray,
    prepared: dict[str, Any],
    plasma_snapshot: dict[str, Any],
    deposited_energy_J: float,
) -> tuple[np.ndarray, list[dict[str, Any]], float, dict[str, float], list[str]]:
    phase = _cached_phase("gri", lambda: ct.Solution(MECHANISM))
    supported = set(phase.species_names)
    candidates = []
    excluded_fragment_channels: list[dict[str, Any]] = []
    for source in plasma_snapshot.get("source_terms", []):
        products = source.get("products", [])
        product_ids = [item["formula"] for item in products]
        if (
            not source.get("reactive")
            or not source.get("atoms_and_product_charge_balanced")
            or source.get("unresolved_dissociation")
            or len(products) < 2
            or source["target"] not in supported
            or any(product_id not in supported for product_id in product_ids)
        ):
            continue
        unresolved_product_state = bool(source.get("product_state_unresolved"))
        excited_product = any(product.get("excited") for product in products)
        charged_product = any(product.get("charge") for product in products)
        annotation_mismatch = not _source_energy_annotations_match_threshold(source)
        if unresolved_product_state or excited_product or charged_product or annotation_mismatch:
            excluded_fragment_channels.append(
                {
                    "channel_id": str(source.get("id")),
                    "product_label": str(source.get("product_label")),
                    "reason": (
                        "unresolved_electronic_product_state" if unresolved_product_state
                        else "unresolved_excited_state" if excited_product
                        else "charged_products_require_explicit_ion_chemistry" if charged_product
                        else "energy_annotation_threshold_mismatch"
                    ),
                    "raw_opportunity_mol": float(
                        source.get("uncoupled_active_window_opportunity_mol", 0.0)
                    ),
                }
            )
            continue
        event_threshold_eV = max(
            float(source.get("threshold_eV") or 0.0),
            float(source.get("first_nonzero_energy_eV") or 0.0),
            0.1,
        )
        energy_eV = 0.0 if source.get("outcome_only") else event_threshold_eV
        raw_extent = float(source["uncoupled_active_window_opportunity_mol"])
        candidates.append(
            {
                "channel_id": source["id"],
                "target": source["target"],
                "products": product_ids,
                "raw_opportunity_mol": raw_extent,
                "energy_cost_eV_per_event": energy_eV,
                "physical_event_threshold_eV": event_threshold_eV,
                "energy_loss_accounted_by_parent_manifold": bool(source.get("outcome_only")),
                "energy_cost_J_mol": energy_eV * FARADAY_J_PER_EV_MOL,
            }
        )

    resolved_fragment_raw_demand_J = sum(
        item["raw_opportunity_mol"] * item["energy_cost_J_mol"] for item in candidates
    )
    all_collision_raw_demand_J = sum(
        max(float(source.get("independent_threshold_energy_lower_bound_J",
                             source.get("threshold_energy_lower_bound_J")) or 0.0), 0.0)
        for source in plasma_snapshot.get("source_terms", [])
    )
    overlapping_outcome_raw_threshold_J = sum(
        max(float(source.get("threshold_energy_lower_bound_J") or 0.0), 0.0)
        for source in plasma_snapshot.get("source_terms", [])
        if source.get("outcome_only")
    )
    energy_scale = (
        min(1.0, deposited_energy_J / all_collision_raw_demand_J)
        if all_collision_raw_demand_J > 0
        else 0.0
    )
    by_target: dict[str, float] = {}
    for item in candidates:
        proposed = item["raw_opportunity_mol"] * energy_scale
        item["energy_scaled_extent_mol"] = proposed
        by_target[item["target"]] = by_target.get(item["target"], 0.0) + proposed
    target_scales = {
        target: min(
            1.0,
            float(gas_amounts_mol[prepared["index"][target]]) / total,
        )
        if total > 0
        else 0.0
        for target, total in by_target.items()
    }

    injected = gas_amounts_mol.copy()
    used_energy_J = 0.0
    accepted = []
    for item in candidates:
        extent = item["energy_scaled_extent_mol"] * target_scales[item["target"]]
        item["accepted_extent_mol"] = extent
        if extent <= 0:
            continue
        injected[prepared["index"][item["target"]]] -= extent
        for product in item["products"]:
            injected[prepared["index"][product]] += extent
        used_energy_J += extent * item["energy_cost_J_mol"]
        accepted.append(item)

    tolerance = max(float(np.max(gas_amounts_mol, initial=0)) * 1e-10, 1e-15)
    if np.min(injected, initial=0) < -tolerance:
        raise ValueError("Plasma neutral-fragment injection produced negative inventory")
    warnings = [
        "Charged collision products are excluded from GRI neutral-fragment injection; the separate charged branch handles supported ion-electron chemistry.",
        "The deposited-energy scale is shared across every imported collision opportunity, including "
        "nonreactive excitation, rather than assigning the full deposit only to resolved fragments.",
    ]
    if excluded_fragment_channels:
        warnings.append(
            f"{len(excluded_fragment_channels)} electron fragment channels with "
            "unresolved excited or unspecified electronic states, charged products, "
            "or inconsistent energy annotations were "
            "excluded from ground-state chemistry."
        )
    if all_collision_raw_demand_J > deposited_energy_J * (1 + 1e-9):
        warnings.append(
            f"All collision opportunities were scaled by {energy_scale:.3e} to respect deposited energy."
        )
    if overlapping_outcome_raw_threshold_J > 0:
        warnings.append(
            "Outcome-only dissociation thresholds are excluded from independent energy "
            "demand because the parent electronic excitation manifold already carries "
            "their energy loss."
        )
    allocated_collision_energy_J = min(deposited_energy_J, all_collision_raw_demand_J)
    allocation = {
        "all_collision_raw_threshold_demand_J": all_collision_raw_demand_J,
        "overlapping_outcome_raw_threshold_J": overlapping_outcome_raw_threshold_J,
        "resolved_fragment_raw_threshold_demand_J": resolved_fragment_raw_demand_J,
        "global_collision_energy_scale": energy_scale,
        "allocated_collision_energy_J": allocated_collision_energy_J,
        "resolved_fragment_energy_used_J": used_energy_J,
        "nonfragment_or_inventory_limited_energy_J": max(
            allocated_collision_energy_J - used_energy_J, 0.0
        ),
        "excluded_fragment_channels": excluded_fragment_channels,
    }
    return np.maximum(injected, 0), accepted, used_energy_J, allocation, warnings


def _create_pavan_neutral_shadow_phase() -> tuple[ct.Solution, int, int, str]:
    """Load Pavan's mechanism as a neutral heavy-chemistry shadow phase.

    Cantera reactors do not currently integrate its ``thermo: plasma`` phase.
    The externally solved electron-collision source terms are therefore injected
    first, while Cantera's electron-collision and two-temperature recombination
    reactions are excluded from this explicitly labelled shadow calculation.
    """

    packaged = next(
        Path(directory) / PAVAN_MECHANISM
        for directory in ct.get_data_directories()
        if (Path(directory) / PAVAN_MECHANISM).is_file()
    )
    plasma = ct.Solution(str(packaged), PAVAN_PHASE)
    excluded_rate_types = {"ElectronCollisionPlasmaRate", "TwoTempPlasmaRate"}
    retained = [
        reaction
        for reaction in plasma.reactions()
        if type(reaction.rate).__name__ not in excluded_rate_types
    ]
    phase = ct.Solution(
        thermo="ideal-gas",
        kinetics="gas",
        species=plasma.species(),
        reactions=retained,
        name="pavan-neutral-shadow",
    )
    return (
        phase,
        len(retained),
        plasma.n_reactions - len(retained),
        hashlib.sha256(packaged.read_bytes()).hexdigest(),
    )


def _pavan_neutral_shadow_phase() -> tuple[ct.Solution, int, int, str]:
    return _cached_phase("pavan-neutral", _create_pavan_neutral_shadow_phase)


def _run_pavan_neutral_shadow(
    injected_neutral_mol: np.ndarray,
    prepared: dict[str, Any],
    plasma_snapshot: dict[str, Any],
    energy_allocation: dict[str, float],
    temperature_K: float,
    gas_volume_L: float,
    duration_s: float,
) -> dict[str, Any]:
    """Run a correlated Pavan-mechanism shadow without mutating the ledger."""

    phase, retained_reactions, excluded_reactions, mechanism_sha256 = (
        _pavan_neutral_shadow_phase()
    )
    supported = set(phase.species_names)
    initial = {
        species: float(injected_neutral_mol[prepared["index"][species]])
        for species in phase.species_names
        if species in prepared["index"]
        and injected_neutral_mol[prepared["index"][species]] > 0
    }
    energy_scale = float(energy_allocation["global_collision_energy_scale"])
    candidates = []
    for source in plasma_snapshot.get("source_terms", []):
        mapped = _PAVAN_EXCITED_N2.get(source.get("product_label", ""))
        if (
            source.get("role") != "excitation"
            or source.get("target") != "N2"
            or not source.get("atoms_and_product_charge_balanced")
            or mapped not in supported
        ):
            continue
        energy_eV = max(
            float(source.get("threshold_eV") or 0.0),
            float(source.get("first_nonzero_energy_eV") or 0.0),
            0.1,
        )
        candidates.append(
            {
                "channel_id": source["id"],
                "target": "N2",
                "product": mapped,
                "source_product_label": source["product_label"],
                "raw_opportunity_mol": float(
                    source["uncoupled_active_window_opportunity_mol"]
                ),
                "energy_cost_eV_per_event": energy_eV,
                "energy_cost_J_mol": energy_eV * FARADAY_J_PER_EV_MOL,
            }
        )

    proposed_by_target: dict[str, float] = {}
    for item in candidates:
        item["energy_scaled_extent_mol"] = item["raw_opportunity_mol"] * energy_scale
        proposed_by_target[item["target"]] = (
            proposed_by_target.get(item["target"], 0.0)
            + item["energy_scaled_extent_mol"]
        )
    target_scales = {
        target: min(1.0, initial.get(target, 0.0) / proposed)
        if proposed > 0
        else 0.0
        for target, proposed in proposed_by_target.items()
    }
    accepted = []
    excitation_energy_used_J = 0.0
    for item in candidates:
        extent = item["energy_scaled_extent_mol"] * target_scales[item["target"]]
        item["accepted_extent_mol"] = extent
        if extent <= 0:
            continue
        initial[item["target"]] = initial.get(item["target"], 0.0) - extent
        initial[item["product"]] = initial.get(item["product"], 0.0) + extent
        excitation_energy_used_J += extent * item["energy_cost_J_mol"]
        accepted.append(item)

    total_mol = sum(initial.values())
    if total_mol <= 0:
        return {
            "status": "no_supported_pavan_species",
            "mechanism": PAVAN_MECHANISM,
            "mechanism_sha256": mechanism_sha256,
            "mapped_excitation_source_extents": accepted,
            "formed_mol": {},
            "warnings": ["No species supported by the Pavan shadow mechanism were present."],
        }

    volume_m3 = gas_volume_L * 1e-3
    mass_kg = sum(
        amount * phase.molecular_weights[phase.species_index(species)] / 1000.0
        for species, amount in initial.items()
    )
    phase.TDX = temperature_K, mass_kg / volume_m3, initial
    initial_pressure_Pa = float(phase.P)
    reactor = ct.IdealGasReactor(phase, energy="off", clone=True)
    reactor.volume = volume_m3
    network = ct.ReactorNet([reactor])
    network.rtol = 1e-9
    network.atol = 1e-18
    network.advance(duration_s)
    final_phase = reactor.phase
    final_total_mol = reactor.mass / final_phase.mean_molecular_weight * 1000.0
    final = {
        species: final_total_mol * float(final_phase[species].X[0])
        for species in phase.species_names
    }
    formed = {
        species: amount - initial.get(species, 0.0)
        for species, amount in final.items()
        if amount - initial.get(species, 0.0) > 1e-14
    }
    consumed = {
        species: initial.get(species, 0.0) - amount
        for species, amount in final.items()
        if initial.get(species, 0.0) - amount > 1e-14
    }
    net_change = {
        species: amount - initial.get(species, 0.0)
        for species, amount in final.items()
        if abs(amount - initial.get(species, 0.0)) > 1e-14
    }
    return {
        "status": "correlated_pavan_neutral_shadow_complete",
        "engine": "Cantera",
        "cantera_version": ct.__version__,
        "mechanism": PAVAN_MECHANISM,
        "phase_source": PAVAN_PHASE,
        "mechanism_sha256": mechanism_sha256,
        "mechanism_species": phase.n_species,
        "state_species_names": list(phase.species_names),
        "retained_heavy_reactions": retained_reactions,
        "excluded_plasma_rate_reactions": excluded_reactions,
        "temperature_K": temperature_K,
        "duration_s": duration_s,
        "initial_pressure_Pa_from_inventory": initial_pressure_Pa,
        "final_pressure_Pa": float(final_phase.P),
        "mapped_excitation_source_extents": accepted,
        "mapped_excitation_energy_used_J": excitation_energy_used_J,
        "input_fragment_energy_used_J": energy_allocation[
            "resolved_fragment_energy_used_J"
        ],
        "total_mapped_source_energy_used_J": (
            excitation_energy_used_J
            + energy_allocation["resolved_fragment_energy_used_J"]
        ),
        "formed_mol": formed,
        "consumed_mol": consumed,
        "net_change_mol": net_change,
        "warnings": [
            "This is a correlated mechanism-extension shadow: Pavan extends GRI-Mech, so agreement is not independent evidence.",
            "The externally solved EEDF supplies energy-scaled neutral fragments and mapped N2 excited states.",
            "This neutral-only comparison excludes electron-collision and two-temperature ion/electron rates; the separate charge-balanced Pavan shadow retains the two-temperature reactions.",
            "Phelps spectroscopic bands are lumped into the closest Pavan N2(A), N2(B), N2(C), or vibrational state; unmapped bands remain energy demand only.",
            "The shadow result cannot mutate the primary apparatus inventory or water-phase yield.",
        ],
        "scientific_firewall": {
            "known_products_accessible_to_solver": True,
            "historical_yields_accessible_to_solver": False,
            "mechanism_tuned_for_miller_urey": False,
            "primary_inventory_mutated_by_shadow": False,
            "quantitative_yield_claim_allowed": False,
        },
    }


def _pavan_source_product_name(
    source: dict[str, Any], product: dict[str, Any], supported: set[str]
) -> str | None:
    """Map a parsed external-EEDF product to the Pavan plasma species set."""

    mapped_excited = _PAVAN_EXCITED_N2.get(source.get("product_label", ""))
    if mapped_excited in supported and len(source.get("products", [])) == 1:
        return mapped_excited
    if product.get("excited"):
        return None
    charge = int(product.get("charge", 0))
    if charge > 0:
        label = str(product.get("label", "")).replace("^", "")
        # LXCat labels such as CO2(+,X) and O2(+,X) denote cations. Taking
        # only the text before '(' would silently map them to neutral CO2/O2
        # while source injection still creates an electron.
        parent = label.split("(", 1)[0]
        candidate = parent if parent.endswith("+") else f"{parent}+"
        return candidate if candidate in supported else None
    if charge < 0:
        return None
    formula = str(product.get("formula", ""))
    return formula if formula in supported else None


def _conservative_clip_pavan_ions(
    raw: np.ndarray, phase: ct.Solution
) -> tuple[np.ndarray, dict[str, float], dict[str, float]]:
    """Project negative molecular cations, exposing sub-resolution roundoff.

    Replacing a negative X+ amount by zero alone creates the atoms of X and
    positive charge. Remove the same tiny amount of neutral X and add an
    electron so the projection preserves Cantera's elemental/E ledger.
    If the same neutral is unavailable, a very small ion (at most 1e-15 mol)
    can instead be compensated across already-positive species. The correction
    must preserve every Cantera element and formal charge. If no such bounded
    correction exists, at most 1e-17 mol remains explicitly reported
    unprojected roundoff; anything larger is rejected.
    This is a bounded numerical correction, not a chemical reaction rate.
    """
    if not np.all(np.isfinite(raw)):
        raise ValueError("Pavan integration returned non-finite species amounts")
    final = np.maximum(raw, 0.0)
    clipped = {
        phase.species_names[position]: float(-amount)
        for position, amount in enumerate(raw)
        if amount < 0
    }
    projected: dict[str, float] = {}
    fallback_ions: dict[str, float] = {}
    electron_position = phase.species_index(phase.electron_species_name)
    for ion, amount in clipped.items():
        if not ion.endswith("+") or phase.species(ion).charge != 1:
            continue
        neutral = ion[:-1]
        if neutral not in phase.species_names:
            raise ValueError(f"Cannot conservatively project unsupported ion {ion}")
        neutral_position = phase.species_index(neutral)
        ion_atoms = {key: value for key, value in phase.species(ion).composition.items() if key != "E"}
        neutral_atoms = {key: value for key, value in phase.species(neutral).composition.items() if key != "E"}
        if ion_atoms != neutral_atoms:
            raise ValueError(f"Ion {ion} does not match neutral parent {neutral}")
        if amount > 2e-11:
            raise ValueError(
                f"Negative ion {ion} amount {amount:.3e} mol exceeds the "
                f"conservative projection bound or neutral {neutral} inventory "
                f"{final[neutral_position]:.3e} mol"
            )
        if final[neutral_position] < amount:
            if amount <= 1e-17:
                continue
            if amount <= 1e-15:
                fallback_ions[ion] = amount
                continue
            raise ValueError(
                f"Negative ion {ion} amount {amount:.3e} mol exceeds the "
                f"conservative projection bound or neutral {neutral} inventory "
                f"{final[neutral_position]:.3e} mol"
            )
        final[neutral_position] -= amount
        final[electron_position] += amount
        projected[ion] = amount
    if fallback_ions:
        element_matrix = np.asarray([
            [phase.n_atoms(position, element) for position in range(phase.n_species)]
            for element in phase.element_names
        ], dtype=float)
        charges = np.asarray([species.charge for species in phase.species()], dtype=float)
        conserved = np.vstack((element_matrix, charges))
        target = -sum(
            (amount * conserved[:, phase.species_index(ion)]
             for ion, amount in fallback_ions.items()),
            np.zeros(conserved.shape[0]),
        )
        available = np.flatnonzero(final > 1e-18)
        while len(available):
            scales = final[available]
            correction = scales * np.linalg.lstsq(
                conserved[:, available] * scales, target, rcond=1e-14,
            )[0]
            feasible = correction >= -scales
            if np.all(feasible):
                residual = conserved[:, available] @ correction - target
                if np.max(np.abs(residual), initial=0.0) <= max(
                    5e-18, float(np.max(np.abs(target), initial=0.0)) * 1e-4,
                ):
                    final[available] += correction
                    projected.update(fallback_ions)
                    break
                raise ValueError(
                    "Pavan ion projection could not close atoms and charge"
                )
            available = available[feasible]
        else:
            raise ValueError(
                "No nonnegative atom/charge-conserving Pavan ion projection exists"
            )
    return final, clipped, projected


def _run_pavan_charged_shadow(
    gas_amounts_mol: np.ndarray,
    prepared: dict[str, Any],
    plasma_snapshot: dict[str, Any],
    energy_allocation: dict[str, float],
    temperature_K: float,
    gas_volume_L: float,
    duration_s: float,
    initial_state_mol: dict[str, float] | None = None,
    coupled_henry_transport: dict[str, dict[str, float]] | None = None,
    boiler_species_mol: dict[str, float] | None = None,
    formaldehyde_hydration: dict[str, Any] | None = None,
    declared_absent_elements: set[str] | None = None,
    integrated_pathway_species: tuple[str, ...] = (),
    diagnostic_reaction_multipliers: tuple[tuple[int, float], ...] = (),
) -> dict[str, Any]:
    """Integrate supported ions/electrons without mutating the apparatus ledger.

    Electron-collision rates embedded in the Pavan mechanism are disabled so
    that the external mixture-specific EEDF is not counted twice. Supported
    source products are applied as a time-averaged, charge-balanced source with
    the same collision-wide energy scale as the primary branch. The 19 Pavan
    two-temperature recombination reactions and all neutral/heavy reactions are
    then integrated directly at fixed gas temperature and volume. Distributing
    each cycle's source over its residence step avoids the unphysical ion spike
    caused by batch injection. This also avoids Cantera's unimplemented plasma-
    phase internal-energy reactor equation while retaining its published rates.
    """

    packaged = next(
        Path(directory) / PAVAN_MECHANISM
        for directory in ct.get_data_directories()
        if (Path(directory) / PAVAN_MECHANISM).is_file()
    )
    phase = _cached_phase(
        "pavan-charged", lambda: ct.Solution(str(packaged), PAVAN_PHASE)
    )
    supported = set(phase.species_names)
    electron_collision_indices = []
    two_temperature_indices = []
    for position, reaction in enumerate(phase.reactions()):
        rate_type = type(reaction.rate).__name__
        if rate_type == "ElectronCollisionPlasmaRate":
            phase.set_multiplier(0.0, position)
            electron_collision_indices.append(position)
        elif rate_type == "TwoTempPlasmaRate":
            two_temperature_indices.append(position)
    multiplier_indices = [index for index, _ in diagnostic_reaction_multipliers]
    if len(multiplier_indices) != len(set(multiplier_indices)):
        raise ValueError("Charged diagnostic reaction multipliers must be unique")
    for index, equation in PAVAN_HCN_DIAGNOSTIC_REACTIONS.items():
        if phase.reaction(index).equation != equation:
            raise ValueError(f"Unexpected Pavan diagnostic reaction at index {index}")
        # The Pavan Solution is thread-cached. Restore the baseline before
        # every cycle so an opt-in stress cannot leak into subsequent runs.
        phase.set_multiplier(1.0, index)
    for index, factor in diagnostic_reaction_multipliers:
        if index not in PAVAN_HCN_DIAGNOSTIC_REACTIONS:
            raise ValueError(f"Unsupported Pavan diagnostic reaction index {index}")
        if not math.isfinite(factor) or factor < 0:
            raise ValueError("Charged diagnostic reaction multipliers must be finite and nonnegative")
        phase.set_multiplier(factor, index)

    baseline = np.zeros(phase.n_species, dtype=float)
    if initial_state_mol is not None:
        for species, amount in initial_state_mol.items():
            if species in supported:
                baseline[phase.species_index(species)] = max(float(amount), 0.0)
    else:
        for species in phase.species_names:
            if species in prepared["index"]:
                baseline[phase.species_index(species)] = max(
                    float(gas_amounts_mol[prepared["index"][species]]), 0.0
                )
    initial = baseline.copy()
    energy_scale = float(energy_allocation["global_collision_energy_scale"])
    candidates: list[dict[str, Any]] = []
    unsupported: list[dict[str, str]] = []
    for source in plasma_snapshot.get("source_terms", []):
        if (
            not source.get("reactive")
            or not source.get("atoms_and_product_charge_balanced")
            or source.get("unresolved_dissociation")
            or source.get("target") not in supported
            or not _source_energy_annotations_match_threshold(source)
        ):
            continue
        if source.get("product_state_unresolved"):
            unsupported.append(
                {
                    "channel_id": str(source.get("id")),
                    "product_label": str(source.get("product_label")),
                    "reason": "electronic product-state branching unresolved",
                }
            )
            continue
        mapped_products: list[str] = []
        product_charge = 0
        mapping_problem = "charged or excited products are absent from the Pavan species set"
        for product in source.get("products", []):
            mapped = _pavan_source_product_name(source, product, supported)
            if mapped is None:
                mapped_products = []
                break
            mapped_species = phase.species(mapped)
            mapped_atoms = {
                element: count for element, count in mapped_species.composition.items()
                if element != "E"
            }
            if (
                int(mapped_species.charge) != int(product.get("charge", 0))
                or mapped_atoms != product.get("composition", {})
            ):
                mapped_products = []
                mapping_problem = "mapped Pavan product would change source charge or atoms"
                break
            mapped_products.append(mapped)
            product_charge += int(product.get("charge", 0))
        if not mapped_products:
            unsupported.append(
                {
                    "channel_id": str(source.get("id")),
                    "product_label": str(source.get("product_label")),
                    "reason": mapping_problem,
                }
            )
            continue
        if (
            len(mapped_products) == 1
            and mapped_products[0] == source["target"]
            and product_charge == 0
        ):
            unsupported.append(
                {
                    "channel_id": str(source.get("id")),
                    "product_label": str(source.get("product_label")),
                    "reason": "the contributor label does not resolve a distinct product species",
                }
            )
            continue
        candidates.append(
            {
                "channel_id": source["id"],
                "target": source["target"],
                "source_product_label": source["product_label"],
                "mapped_products": mapped_products,
                "product_charge": product_charge,
                "electron_delta_per_event": product_charge,
                "raw_opportunity_mol": float(
                    source["uncoupled_active_window_opportunity_mol"]
                ),
                "energy_scaled_extent_mol": float(
                    source["uncoupled_active_window_opportunity_mol"]
                )
                * energy_scale,
                "independent_threshold_energy_cost_J_mol": (
                    max(float(source.get(
                        "independent_threshold_energy_lower_bound_J",
                        source.get("threshold_energy_lower_bound_J", 0.0),
                    )), 0.0)
                    / float(source["uncoupled_active_window_opportunity_mol"])
                    if float(source["uncoupled_active_window_opportunity_mol"]) > 0
                    else 0.0
                ),
            }
        )

    proposed_by_target: dict[str, float] = {}
    for item in candidates:
        proposed_by_target[item["target"]] = (
            proposed_by_target.get(item["target"], 0.0)
            + item["energy_scaled_extent_mol"]
        )
    target_scales = {
        target: min(
            1.0,
            initial[phase.species_index(target)] / proposed,
        )
        if proposed > 0
        else 0.0
        for target, proposed in proposed_by_target.items()
    }
    accepted = []
    electron_index = phase.species_index(phase.electron_species_name)
    source_vector_mol_s = np.zeros(phase.n_species, dtype=float)
    source_channels: list[tuple[int, float, np.ndarray]] = []
    for item in candidates:
        extent = item["energy_scaled_extent_mol"] * target_scales[item["target"]]
        item["nominal_inventory_capped_extent_mol"] = extent
        item["nominal_time_averaged_source_rate_mol_s"] = (
            extent / duration_s if duration_s > 0 else 0.0
        )
        if extent <= 0:
            continue
        target_position = phase.species_index(item["target"])
        rate = item["nominal_time_averaged_source_rate_mol_s"]
        channel_stoichiometry = np.zeros(phase.n_species, dtype=float)
        channel_stoichiometry[target_position] -= 1.0
        for product in item["mapped_products"]:
            channel_stoichiometry[phase.species_index(product)] += 1.0
        channel_stoichiometry[electron_index] += item["electron_delta_per_event"]
        source_vector_mol_s += channel_stoichiometry * rate
        source_channels.append((
            target_position, rate / initial[target_position], channel_stoichiometry,
        ))
        accepted.append(item)

    charges = np.asarray([float(species.charge) for species in phase.species()])
    source_charge_residual_rate = float(np.dot(charges, source_vector_mol_s))
    source_charge_residual = source_charge_residual_rate * max(duration_s, 0.0)
    if abs(source_charge_residual) > max(float(initial.max()) * 1e-9, 1e-14):
        raise ValueError("Charged Pavan source injection failed charge closure")
    if initial.sum() <= 0 or duration_s <= 0:
        return {
            "status": "no_supported_charged_pavan_species",
            "mechanism": PAVAN_MECHANISM,
            "formed_mol": {},
            "source_extents": accepted,
            "unsupported_source_channels": unsupported,
            "warnings": ["No supported charged Pavan state was available to integrate."],
        }

    # An absent element cannot be created by either the external collisions or
    # heavy chemistry. Keep its species exactly zero so BDF roundoff cannot
    # seed a trace ion in feed-removal controls (as in the primary GRI branch).
    element_matrix = np.asarray(
        [
            [phase.n_atoms(species_index, element) for species_index in range(phase.n_species)]
            for element in phase.element_names if element != "E"
        ],
        dtype=float,
    )
    supplied_atoms = element_matrix @ (
        initial + np.maximum(source_vector_mol_s, 0.0) * duration_s
    )
    for element in declared_absent_elements or ():
        if element not in phase.element_names or element == "E":
            continue
        position = [name for name in phase.element_names if name != "E"].index(element)
        if supplied_atoms[position] != 0.0:
            raise ValueError(
                f"Declared absent element {element} entered the charged Pavan input/source "
                f"at {supplied_atoms[position]:.3e} mol atoms"
            )
    inactive_species = np.any(element_matrix[supplied_atoms == 0.0] > 0, axis=0)
    if np.any(initial[inactive_species] != 0.0) or np.any(source_vector_mol_s[inactive_species] != 0.0):
        raise ValueError("Inactive-element mask conflicts with supplied charged gas or fragments")

    volume_m3 = gas_volume_L * 1e-3
    molecular_weights = phase.molecular_weights

    def set_state(amounts: np.ndarray) -> None:
        bounded = np.maximum(amounts, 0.0)
        bounded[inactive_species] = 0.0
        mass_kg = float(np.dot(bounded, molecular_weights) / 1000.0)
        phase.TDX = temperature_K, mass_kg / volume_m3, bounded

    set_state(initial)
    phase.electron_energy_distribution_type = "isotropic"
    mean_electron_energy_eV = max(
        float(plasma_snapshot.get("mean_electron_energy_eV") or 0.0), 1e-6
    )
    phase.mean_electron_energy = mean_electron_energy_eV

    tracked = [
        (species, phase.species_index(species), data)
        for species, data in (coupled_henry_transport or {}).items()
        if species in supported
    ]
    tracked_count = len(tracked)
    initial_boiler = np.asarray(
        [float((boiler_species_mol or {}).get(species, 0.0)) for species, _, _ in tracked],
        dtype=float,
    )
    hydration_state_start = phase.n_species + 2 * tracked_count
    initial_methanediol_mol = 0.0
    if formaldehyde_hydration is not None:
        if not any(species == "CH2O" for species, _, _ in tracked) or "H2O" not in supported:
            raise ValueError("Charged hydration needs tracked CH2O and gas H2O identities")
        if boiler_species_mol is None or float(boiler_species_mol.get("H2O", 0.0)) <= 0:
            raise ValueError("Charged hydration needs positive liquid-water inventory")
        initial_methanediol_mol = float(boiler_species_mol.get(METHANEDIOL_AQ, 0.0))
    initial_integrated_state = (
        np.concatenate([initial, initial_boiler, np.zeros(tracked_count)])
        if tracked_count else initial
    )
    if formaldehyde_hydration is not None:
        initial_integrated_state = np.concatenate([
            initial_integrated_state,
            [initial_methanediol_mol, float(boiler_species_mol["H2O"])],
        ])
    source_extent_start = initial_integrated_state.size
    if source_channels:
        initial_integrated_state = np.concatenate([
            initial_integrated_state, np.zeros(len(source_channels), dtype=float)
        ])

    def rhs(_time_s: float, amounts: np.ndarray) -> np.ndarray:
        bounded = np.maximum(amounts[:phase.n_species], 0.0)
        bounded[inactive_species] = 0.0
        set_state(bounded)
        chemistry = (
            np.asarray(phase.net_production_rates, dtype=float) * 1000.0 * volume_m3
        )
        source_rates = np.asarray(
            [coefficient * bounded[target] for target, coefficient, _ in source_channels],
            dtype=float,
        )
        gas_change = chemistry.copy()
        for (_, _, stoichiometry), rate in zip(source_channels, source_rates, strict=True):
            gas_change += stoichiometry * rate
        gas_change[inactive_species] = 0.0
        if not tracked_count:
            return np.concatenate([gas_change, source_rates]) if source_channels else gas_change
        boiler_change = np.zeros(tracked_count, dtype=float)
        captured_change = np.zeros(tracked_count, dtype=float)
        hydration_change = np.zeros(2, dtype=float)
        for offset, (species, position, data) in enumerate(tracked):
            gas_value = float(bounded[position])
            liquid_value = max(float(amounts[phase.n_species + offset]), 0.0)
            if formaldehyde_hydration is not None and species == "CH2O":
                if data["condenser_capture_rate_per_chemistry_s"] != 0:
                    raise ValueError("Closed-flask formaldehyde hydrate cannot use condenser capture")
                changes = formaldehyde_contact_hydration_rates_mol_s(
                    gas_value, liquid_value,
                    max(float(amounts[hydration_state_start]), 0.0),
                    max(float(amounts[hydration_state_start + 1]), 0.0),
                    free_liquid_to_gas_ratio=float(formaldehyde_hydration["free_liquid_to_gas_ratio"]),
                    liquid_side_kla_per_s=float(formaldehyde_hydration["liquid_side_kla_per_s"]),
                    hydration_rate_per_s=float(formaldehyde_hydration["hydration_rate_per_s"]),
                    dehydration_rate_per_s=float(formaldehyde_hydration["dehydration_rate_per_s"]),
                )
                gas_change[position] += changes[0]
                boiler_change[offset] = changes[1]
                hydration_change[:] = changes[2:]
                continue
            ratio = data["boiler_liquid_to_gas_equilibrium_ratio"]
            target_gas = (gas_value + liquid_value) / (1.0 + ratio)
            exchange = data["boiler_relaxation_rate_per_chemistry_s"] * (
                gas_value - target_gas
            )
            capture = data["condenser_capture_rate_per_chemistry_s"] * gas_value
            gas_change[position] -= exchange + capture
            boiler_change[offset] = exchange
            captured_change[offset] = capture
        return np.concatenate([
            gas_change, boiler_change, captured_change,
            hydration_change if formaldehyde_hydration is not None else np.zeros(0),
            source_rates,
        ])

    # Captured-liquid and integrated source extents are write-only ledger
    # coordinates. Inactive-element species have an identically zero RHS.
    # SciPy's adaptive finite-difference Jacobian repeatedly enlarges the
    # perturbation of these zero columns and can overflow its step multiplier
    # on long wet-flask runs. Supply the exact zero columns and difference only
    # coordinates that can influence a rate; do not change the ODE or tolerances.
    dynamic_columns = [
        position for position in range(phase.n_species)
        if not inactive_species[position]
    ]
    dynamic_columns.extend(range(phase.n_species, phase.n_species + tracked_count))
    if formaldehyde_hydration is not None:
        dynamic_columns.extend(range(hydration_state_start, hydration_state_start + 2))

    def bounded_jacobian(time_s: float, amounts: np.ndarray) -> np.ndarray:
        return _bounded_rate_jacobian(rhs, time_s, amounts, dynamic_columns)

    # A stiff solve can overshoot an exhausted neutral species. Clipping that
    # appreciable negative inventory would manufacture atoms, so retry with a
    # smaller step before the existing conservative projection/closure checks.
    neutral_positions = [
        position for position in range(phase.n_species)
        if phase.species(position).charge == 0 and not inactive_species[position]
    ]
    solved_state: np.ndarray | None = None
    negative_neutral: tuple[str, float] | None = None
    for rtol, atol, step_divisor in (
        (1e-10, 1e-22, 40.0),
        (1e-11, 1e-24, 100.0),
        (1e-12, 1e-26, 200.0),
    ):
        solution = solve_ivp(
            rhs,
            (0.0, duration_s),
            initial_integrated_state,
            method="BDF",
            jac=bounded_jacobian,
            rtol=rtol,
            atol=atol,
            max_step=max(duration_s / step_divisor, 1e-9),
            dense_output=bool(integrated_pathway_species),
        )
        if not solution.success:
            raise ValueError(f"Charged Pavan integration failed: {solution.message}")
        candidate_state = np.asarray(solution.y[:, -1], dtype=float)
        worst_position = min(
            neutral_positions, key=lambda position: candidate_state[position]
        )
        worst_amount = float(candidate_state[worst_position])
        if worst_amount >= -1e-13:
            solved_state = candidate_state
            break
        negative_neutral = (phase.species_names[worst_position], worst_amount)
    if solved_state is None:
        raise ValueError(
            "Charged Pavan integration exhausted neutral inventory after "
            f"bounded solver retries: {negative_neutral}"
        )
    for offset, item in enumerate(accepted):
        realized_extent = float(solved_state[source_extent_start + offset])
        if realized_extent < -1e-18:
            raise ValueError("Charged Pavan integrated a negative electron-collision extent")
        item["accepted_extent_mol"] = max(realized_extent, 0.0)
        item["time_averaged_source_rate_mol_s"] = item["accepted_extent_mol"] / duration_s
    realized_independent_threshold_energy_J = sum(
        item["accepted_extent_mol"] * item["independent_threshold_energy_cost_J_mol"]
        for item in accepted
    )
    allocated_collision_energy_J = float(energy_allocation["allocated_collision_energy_J"])
    if realized_independent_threshold_energy_J > allocated_collision_energy_J * (1 + 1e-7):
        raise ValueError("Charged Pavan source exceeded collision threshold-energy allocation")
    if formaldehyde_hydration is not None and np.min(solved_state[hydration_state_start:hydration_state_start + 2]) < -1e-18:
        raise ValueError("Charged formaldehyde hydration produced negative hydrate or water inventory")
    raw_solved_gas = solved_state[:phase.n_species].copy()
    inactive_roundoff = {
        phase.species_names[position]: float(raw_solved_gas[position])
        for position in np.flatnonzero(inactive_species)
        if raw_solved_gas[position] != 0.0
    }
    if max((abs(value) for value in inactive_roundoff.values()), default=0.0) > 1e-18:
        raise ValueError("Charged Pavan absent-element roundoff exceeded 1e-18 mol")
    # The exact solution has zero inventory for these species. BDF's linear
    # solve can leave sub-attomole residue even with a zero derivative; do not
    # carry that residue into the next apparatus cycle as a new element source.
    solved_state[:phase.n_species][inactive_species] = 0.0
    for offset, (_species, position, _) in enumerate(tracked):
        if inactive_species[position] and initial_boiler[offset] == 0.0:
            if max(
                abs(solved_state[phase.n_species + offset]),
                abs(solved_state[phase.n_species + tracked_count + offset]),
            ) > 1e-18:
                raise ValueError("Charged Pavan absent-element transport roundoff exceeded 1e-18 mol")
            solved_state[phase.n_species + offset] = 0.0
            solved_state[phase.n_species + tracked_count + offset] = 0.0
    final, clipped_negative_species, ion_clipping_projection = (
        _conservative_clip_pavan_ions(solved_state[:phase.n_species], phase)
    )
    integrated_pathway_diagnostics: dict[str, Any] | None = None
    if integrated_pathway_species:
        unknown = set(integrated_pathway_species) - supported
        if unknown:
            raise ValueError(f"Unsupported charged pathway species: {sorted(unknown)}")
        reaction_extents = np.zeros(phase.n_reactions, dtype=float)
        abscissae, weights = np.polynomial.legendre.leggauss(5)
        assert solution.sol is not None
        for left, right in zip(solution.t[:-1], solution.t[1:], strict=True):
            midpoint = 0.5 * (left + right)
            half_width = 0.5 * (right - left)
            for abscissa, weight in zip(abscissae, weights, strict=True):
                amounts_at_time = np.asarray(
                    solution.sol(midpoint + half_width * abscissa), dtype=float
                )[:phase.n_species]
                set_state(amounts_at_time)
                reaction_extents += (
                    half_width * weight
                    * np.asarray(phase.net_rates_of_progress, dtype=float)
                    * 1000.0 * volume_m3
                )
        net_stoichiometry = phase.product_stoich_coeffs - phase.reactant_stoich_coeffs
        equations = phase.reaction_equations()
        pathway_species: dict[str, Any] = {}
        for species in integrated_pathway_species:
            position = phase.species_index(species)
            contributions = np.asarray(
                net_stoichiometry[position] * reaction_extents, dtype=float
            )
            ranked = [
                {
                    "reaction_index": int(reaction_position),
                    "equation": equations[reaction_position],
                    "species_net_mol": float(value),
                }
                for reaction_position, value in enumerate(contributions)
                if abs(value) > 1e-30
            ]
            direct_source_mol = sum(
                float(item["accepted_extent_mol"])
                * (item["mapped_products"].count(species) - int(item["target"] == species))
                for item in accepted
            )
            boiler_delta_mol = captured_mol = 0.0
            # In the closed-flask hydration branch, boiler CH2O converts to
            # methanediol. Include that tracked chemical pool when isolating
            # *gas reaction* progress from transport and liquid hydration.
            methanediol_delta_mol = (
                float(solved_state[hydration_state_start]) - initial_methanediol_mol
                if species == "CH2O" and formaldehyde_hydration is not None else 0.0
            )
            for offset, (tracked_species, _, _) in enumerate(tracked):
                if tracked_species == species:
                    boiler_delta_mol = float(
                        solved_state[phase.n_species + offset] - initial_boiler[offset]
                    )
                    captured_mol = float(
                        solved_state[phase.n_species + tracked_count + offset]
                    )
                    break
            observed_chemical_net_mol = float(
                final[position] - initial[position] - direct_source_mol
                + boiler_delta_mol + captured_mol + methanediol_delta_mol
            )
            integrated_chemical_net_mol = float(np.sum(contributions))
            pathway_species[species] = {
                "integrated_chemical_net_mol": integrated_chemical_net_mol,
                "observed_chemical_net_mol": observed_chemical_net_mol,
                "quadrature_closure_residual_mol": (
                    integrated_chemical_net_mol - observed_chemical_net_mol
                ),
                "direct_electron_source_mol": direct_source_mol,
                "coupled_boiler_delta_mol": boiler_delta_mol,
                "coupled_capture_mol": captured_mol,
                "coupled_methanediol_delta_mol": methanediol_delta_mol,
                "integrated_formation_mol": float(
                    sum(item["species_net_mol"] for item in ranked
                        if item["species_net_mol"] > 0)
                ),
                "integrated_consumption_mol": float(
                    -sum(item["species_net_mol"] for item in ranked
                         if item["species_net_mol"] < 0)
                ),
                "formation": sorted(
                    (item for item in ranked if item["species_net_mol"] > 0),
                    key=lambda item: -item["species_net_mol"],
                )[:7],
                "consumption": sorted(
                    (item for item in ranked if item["species_net_mol"] < 0),
                    key=lambda item: item["species_net_mol"],
                )[:7],
                "reaction_contributions": ranked,
            }
        integrated_pathway_diagnostics = {
            "method": "5-point Gauss quadrature on adaptive charged BDF intervals",
            "adaptive_interval_count": int(len(solution.t) - 1),
            "species": pathway_species,
            "not_a_validated_yield_or_mechanism": True,
        }
    set_state(final)
    net_rates_of_progress = np.asarray(phase.net_rates_of_progress, dtype=float)
    reaction_equations = phase.reaction_equations()
    pathway_diagnostics: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for tracked_species in ("HCN", "H2CN", "CH2NH", "CH2O", "CH3OH"):
        if tracked_species not in supported:
            continue
        species_position = phase.species_index(tracked_species)
        contributions = []
        for reaction_position, net_rate in enumerate(net_rates_of_progress):
            stoichiometric_coefficient = (
                phase.product_stoich_coeff(species_position, reaction_position)
                - phase.reactant_stoich_coeff(species_position, reaction_position)
            )
            contribution_mol_s = (
                stoichiometric_coefficient
                * float(net_rate)
                * 1000.0
                * volume_m3
            )
            if abs(contribution_mol_s) <= 1e-30:
                continue
            contributions.append(
                {
                    "reaction_index": reaction_position,
                    "equation": reaction_equations[reaction_position],
                    "rate_type": type(phase.reaction(reaction_position).rate).__name__,
                    "species_net_rate_mol_s": contribution_mol_s,
                }
            )
        pathway_diagnostics[tracked_species] = {
            "formation": sorted(
                (item for item in contributions if item["species_net_rate_mol_s"] > 0),
                key=lambda item: -item["species_net_rate_mol_s"],
            )[:5],
            "consumption": sorted(
                (item for item in contributions if item["species_net_rate_mol_s"] < 0),
                key=lambda item: item["species_net_rate_mol_s"],
            )[:5],
        }
    elements = phase.element_names
    atom_matrix = np.asarray(
        [
            [phase.n_atoms(species_index, element) for species_index in range(phase.n_species)]
            for element in elements
        ],
        dtype=float,
    )
    transport_vector = np.zeros(phase.n_species, dtype=float)
    coupled_transport_report: dict[str, Any] | None = None
    if tracked_count:
        boiler_delta = {}
        captured = {}
        for offset, (species, position, _) in enumerate(tracked):
            delta_mol = float(solved_state[phase.n_species + offset] - initial_boiler[offset])
            captured_mol = float(max(solved_state[phase.n_species + tracked_count + offset], 0.0))
            boiler_delta[species] = delta_mol
            captured[species] = captured_mol
            transport_vector[position] += delta_mol + captured_mol
        if formaldehyde_hydration is not None:
            hydrate_delta = float(
                solved_state[hydration_state_start] - boiler_species_mol[METHANEDIOL_AQ]
            )
            water_delta = float(
                solved_state[hydration_state_start + 1] - boiler_species_mol["H2O"]
            )
            boiler_delta[METHANEDIOL_AQ] = hydrate_delta
            boiler_delta["H2O"] = water_delta
            transport_vector[phase.species_index("CH2O")] += hydrate_delta
            transport_vector[phase.species_index("H2O")] += hydrate_delta + water_delta
        coupled_transport_report = {
            "status": "mean_field_henry_exchange_coupled_to_pavan_bdf",
            "boiler_delta_mol": boiler_delta,
            "captured_mol": captured,
        }
    raw_atom_residual = atom_matrix @ (
        raw_solved_gas - initial + transport_vector
    )
    atom_residual = atom_matrix @ (final - initial + transport_vector)
    final_charge_residual = float(np.dot(charges, final))
    charge_tolerance = max(float(final.max()) * 1e-8, 1e-13)
    atom_tolerance = max(float((atom_matrix @ initial).max(initial=0.0)) * 1e-8, 1e-13)
    if abs(final_charge_residual) > charge_tolerance:
        raise ValueError("Charged Pavan integration failed charge closure")
    if float(np.max(np.abs(atom_residual), initial=0.0)) > atom_tolerance:
        raise ValueError(
            "Charged Pavan integration failed atom closure: "
            f"max_residual_mol={float(np.max(np.abs(atom_residual), initial=0.0)):.6g}, "
            f"tolerance_mol={atom_tolerance:.6g}, "
            f"residual_by_element_mol="
            f"{dict(zip(elements, map(float, atom_residual), strict=True))}, "
            f"raw_residual_by_element_mol="
            f"{dict(zip(elements, map(float, raw_atom_residual), strict=True))}, "
            f"largest_clipped_species_mol="
            f"{sorted(clipped_negative_species.items(), key=lambda item: -item[1])[:5]}"
        )

    delta = final - baseline
    formed = {
        species: float(delta[phase.species_index(species)])
        for species in phase.species_names
        if delta[phase.species_index(species)] > 1e-14
    }
    consumed = {
        species: float(-delta[phase.species_index(species)])
        for species in phase.species_names
        if delta[phase.species_index(species)] < -1e-14
    }
    net_change = {
        species: float(delta[phase.species_index(species)])
        for species in phase.species_names
        if abs(delta[phase.species_index(species)]) > 1e-14
    }
    return {
        "status": "correlated_pavan_charge_balanced_shadow_complete",
        "engine": "Cantera rates + SciPy BDF fixed-T/V plasma state",
        "cantera_version": ct.__version__,
        "mechanism": PAVAN_MECHANISM,
        "phase_source": PAVAN_PHASE,
        "mechanism_sha256": hashlib.sha256(packaged.read_bytes()).hexdigest(),
        "diagnostic_reaction_multipliers": {
            PAVAN_HCN_DIAGNOSTIC_REACTIONS[index]: factor
            for index, factor in diagnostic_reaction_multipliers
        },
        "diagnostic_rate_stress_is_validated_kinetics": False,
        "mechanism_species": phase.n_species,
        "state_species_names": list(phase.species_names),
        "retained_reactions": phase.n_reactions - len(electron_collision_indices),
        "retained_two_temperature_reactions": len(two_temperature_indices),
        "disabled_embedded_electron_collision_reactions": len(
            electron_collision_indices
        ),
        "temperature_K": temperature_K,
        "duration_s": duration_s,
        "coupled_henry_transport": coupled_transport_report,
        "mean_electron_energy_eV": mean_electron_energy_eV,
        "electron_energy_closure": "isotropic distribution with externally solved mean energy; embedded electron collisions disabled",
        "source_extents": accepted,
        "realized_independent_threshold_energy_J": realized_independent_threshold_energy_J,
        "allocated_collision_energy_J": allocated_collision_energy_J,
        "threshold_energy_is_full_thermal_balance": False,
        "unsupported_source_channels": unsupported,
        "source_charge_residual_mol_e": source_charge_residual,
        "source_charge_residual_mol_e_per_s": source_charge_residual_rate,
        "final_charge_residual_mol_e": final_charge_residual,
        "max_atom_residual_mol": float(np.max(np.abs(atom_residual), initial=0.0)),
        "max_raw_atom_residual_mol": float(np.max(np.abs(raw_atom_residual), initial=0.0)),
        "atom_residual_by_element_mol": {
            element: float(value)
            for element, value in zip(elements, atom_residual, strict=True)
        },
        "clipped_negative_species_mol": clipped_negative_species,
        "inactive_element_roundoff_projection_mol": inactive_roundoff,
        "conservative_ion_clipping_projection_mol": ion_clipping_projection,
        "unprojected_ion_roundoff_clip_mol": {
            species: amount for species, amount in clipped_negative_species.items()
            if species.endswith("+") and species not in ion_clipping_projection
        },
        "solver": {
            "method": "SciPy BDF",
            "success": bool(solution.success),
            "function_evaluations": int(solution.nfev),
            "relative_tolerance": rtol,
            "absolute_tolerance_mol": atol,
        },
        "formed_mol": formed,
        "consumed_mol": consumed,
        "net_change_mol": net_change,
        "final_instantaneous_pathway_diagnostics": pathway_diagnostics,
        "integrated_pathway_diagnostics": integrated_pathway_diagnostics,
        "final_state_mol": {
            species: float(final[phase.species_index(species)])
            for species in phase.species_names
            if final[phase.species_index(species)] > 0
        },
        "warnings": [
            "This is a correlated Pavan shadow and cannot mutate the primary apparatus inventory or water-phase yield.",
            "External mixture-EEDF source terms replace the mechanism's embedded electron-collision rates to avoid double counting.",
            "Nominal collision extents set first-order, target-inventory-limited source rates during the fixed-T/V residence step; reported accepted extents are integrated, and individual discharge pulses are not resolved.",
            "Only ions and excited states explicitly present in the Pavan species set are injected; unsupported anions and fragment ions remain reported but unresolved.",
            "Excited products without an explicit matching Pavan state are excluded; no unmodelled relaxation to a ground-state formula is performed.",
            "The fixed-temperature/volume plasma ODE retains published two-temperature recombination rates but does not solve electron energy or circuit dynamics during the step.",
        ],
        "scientific_firewall": {
            "known_products_accessible_to_solver": True,
            "historical_yields_accessible_to_solver": False,
            "embedded_electron_collision_rates_used": False,
            "external_eedf_sources_double_counted": False,
            "primary_inventory_mutated_by_shadow": False,
            "quantitative_yield_claim_allowed": False,
        },
    }


def run_pavan_charged_species_cycle(
    gas_amounts_mol: np.ndarray,
    prepared: dict[str, Any],
    plasma_snapshot: dict[str, Any],
    temperature_K: float,
    gas_volume_L: float,
    duration_s: float,
    deposited_energy_J: float,
    initial_state_mol: dict[str, float] | None = None,
    coupled_henry_transport: dict[str, dict[str, float]] | None = None,
    boiler_species_mol: dict[str, float] | None = None,
    formaldehyde_hydration: dict[str, Any] | None = None,
    declared_absent_elements: set[str] | None = None,
    integrated_pathway_species: tuple[str, ...] = (),
    diagnostic_reaction_multipliers: tuple[tuple[int, float], ...] = (),
) -> dict[str, Any]:
    """Advance one independent charged-Pavan branch cycle.

    The collision-wide energy allocation is recomputed from this branch's own
    self-consistent EEDF snapshot. The returned final state can be passed back
    on the next apparatus cycle; it never mutates the primary GRI inventory.
    """

    _, _, _, energy_allocation, _ = _neutral_fragment_extents(
        gas_amounts_mol, prepared, plasma_snapshot, deposited_energy_J
    )
    result = _run_pavan_charged_shadow(
        gas_amounts_mol,
        prepared,
        plasma_snapshot,
        energy_allocation,
        temperature_K,
        gas_volume_L,
        duration_s,
        initial_state_mol=initial_state_mol,
        coupled_henry_transport=coupled_henry_transport,
        boiler_species_mol=boiler_species_mol,
        formaldehyde_hydration=formaldehyde_hydration,
        declared_absent_elements=declared_absent_elements,
        integrated_pathway_species=integrated_pathway_species,
        diagnostic_reaction_multipliers=diagnostic_reaction_multipliers,
    )
    result["collision_energy_allocation"] = energy_allocation
    result["propagated_across_apparatus_cycles"] = initial_state_mol is not None
    return result


def pavan_internal_conserved_mol(state_mol: dict[str, float]) -> dict[str, float]:
    """Return element and charge inventories for Pavan-only carried species."""

    packaged = next(
        Path(directory) / PAVAN_MECHANISM
        for directory in ct.get_data_directories()
        if (Path(directory) / PAVAN_MECHANISM).is_file()
    )
    phase = _cached_phase(
        "pavan-charged", lambda: ct.Solution(str(packaged), PAVAN_PHASE)
    )
    conserved = {element: 0.0 for element in phase.element_names}
    conserved["charge"] = 0.0
    for species, amount in state_mol.items():
        if species not in phase.species_names or amount == 0:
            continue
        position = phase.species_index(species)
        for element in phase.element_names:
            conserved[element] += float(amount) * phase.n_atoms(position, element)
        conserved["charge"] += float(amount) * float(phase.species(position).charge)
    return conserved


def _create_kamarinopoulou_neutral_shadow_phase() -> tuple[ct.Solution, int, int, str]:
    """Load the published CH4/N2 network with GRI species added as inert bath gas."""

    mechanism_sha256 = hashlib.sha256(KAMARINOPOULOU_MECHANISM.read_bytes()).hexdigest()
    if mechanism_sha256 != KAMARINOPOULOU_MECHANISM_SHA256:
        raise ValueError("Kamarinopoulou 2024 mechanism SHA-256 verification failed")
    published = ct.Solution(str(KAMARINOPOULOU_MECHANISM))
    species = list(published.species())
    published_names = set(published.species_names)
    # The paper's model is CHN-only. Retaining unsupported GRI inventory as
    # chemically inert bath gas preserves H2O dilution and third-body density.
    # GRI CH2 is mapped to the paper's triplet CH2(T) below instead of being
    # added as a separate inert species.
    bath_species = [
        item
        for item in ct.Solution(MECHANISM).species()
        if item.name not in published_names and item.name != "CH2"
    ]
    phase = ct.Solution(
        thermo="ideal-gas",
        kinetics="gas",
        species=species + bath_species,
        reactions=published.reactions(),
        name="kamarinopoulou-2024-neutral-shadow",
    )
    return phase, published.n_species, len(bath_species), mechanism_sha256


def _kamarinopoulou_neutral_shadow_phase() -> tuple[ct.Solution, int, int, str]:
    mechanism_sha256 = hashlib.sha256(KAMARINOPOULOU_MECHANISM.read_bytes()).hexdigest()
    if mechanism_sha256 != KAMARINOPOULOU_MECHANISM_SHA256:
        raise ValueError("Kamarinopoulou 2024 mechanism SHA-256 verification failed")
    return _cached_phase("kamarinopoulou-neutral", _create_kamarinopoulou_neutral_shadow_phase)


def _run_kamarinopoulou_neutral_shadow(
    injected_neutral_mol: np.ndarray,
    prepared: dict[str, Any],
    temperature_K: float,
    gas_volume_L: float,
    duration_s: float,
) -> dict[str, Any]:
    """Run the full published neutral CH4/N2 network as a non-mutating shadow."""

    phase, published_species, bath_species, mechanism_sha256 = (
        _kamarinopoulou_neutral_shadow_phase()
    )
    initial: dict[str, float] = {}
    source_mapping: dict[str, str] = {}
    for species in phase.species_names:
        source = "CH2" if species == "CH2(T)" else species
        if source not in prepared["index"]:
            continue
        amount = float(injected_neutral_mol[prepared["index"][source]])
        if amount <= 0:
            continue
        initial[species] = amount
        if source != species:
            source_mapping[source] = species

    if not initial:
        return {
            "status": "no_supported_kamarinopoulou_species",
            "mechanism": KAMARINOPOULOU_MECHANISM.name,
            "mechanism_sha256": mechanism_sha256,
            "formed_mol": {},
            "net_change_mol": {},
            "warnings": ["No species supported by the literature HCN shadow were present."],
        }

    volume_m3 = gas_volume_L * 1e-3
    mass_kg = sum(
        amount * phase.molecular_weights[phase.species_index(species)] / 1000.0
        for species, amount in initial.items()
    )
    phase.TDX = temperature_K, mass_kg / volume_m3, initial
    initial_pressure_Pa = float(phase.P)
    reactor = ct.IdealGasReactor(phase, energy="off", clone=True)
    reactor.volume = volume_m3
    network = ct.ReactorNet([reactor])
    network.rtol = 1e-9
    network.atol = 1e-18
    integrated_duration_s = min(
        duration_s, KAMARINOPOULOU_PUBLISHED_RESIDENCE_TIME_S
    )
    network.max_time_step = max(integrated_duration_s / 20.0, 1e-9)
    numerical_retry_used = False
    try:
        network.advance(integrated_duration_s)
    except ct.CanteraError as first_error:
        retry_phase, _, _, _ = _kamarinopoulou_neutral_shadow_phase()
        retry_phase.TDX = temperature_K, mass_kg / volume_m3, initial
        retry_reactor = ct.IdealGasReactor(retry_phase, energy="off", clone=True)
        retry_reactor.volume = volume_m3
        retry_network = ct.ReactorNet([retry_reactor])
        retry_network.rtol = 1e-7
        retry_network.atol = 1e-16
        retry_network.max_time_step = max(integrated_duration_s / 100.0, 1e-10)
        try:
            retry_network.advance(integrated_duration_s)
            reactor = retry_reactor
            network = retry_network
            numerical_retry_used = True
        except ct.CanteraError as retry_error:
            return {
                "status": "kamarinopoulou_2024_neutral_shadow_failed",
                "engine": "Cantera",
                "cantera_version": ct.__version__,
                "mechanism": KAMARINOPOULOU_MECHANISM.name,
                "mechanism_sha256": mechanism_sha256,
                "published_species": published_species,
                "retained_heavy_reactions": phase.n_reactions,
                "inert_bath_species_added": bath_species,
                "fixed_condition_electron_reactions_excluded": 62,
                "temperature_K": temperature_K,
                "requested_duration_s": duration_s,
                "integrated_duration_s": integrated_duration_s,
                "published_residence_time_s": KAMARINOPOULOU_PUBLISHED_RESIDENCE_TIME_S,
                "initial_pressure_Pa_from_inventory": initial_pressure_Pa,
                "initial_mole_fractions": {
                    species: float(phase[species].X[0])
                    for species in initial
                    if phase[species].X[0] > 1e-12
                },
                "formed_mol": {},
                "consumed_mol": {},
                "net_change_mol": {},
                "error": _cantera_error_summary(retry_error),
                "first_attempt_error": _cantera_error_summary(first_error),
                "warnings": [
                    "The literature shadow integrator failed both bounded numerical attempts for this cycle; the primary GRI inventory was not affected.",
                    "No zero or fabricated literature-shadow result is used as a yield or probability.",
                    "The published HNC thermo record conflicts with evaluated ground-state HNC enthalpy; HNC-sensitive shadow rates are not independently validated.",
                ],
                "scientific_firewall": {
                    "primary_inventory_mutated_by_shadow": False,
                    "quantitative_yield_claim_allowed": False,
                    "hnc_ground_state_thermo_independently_validated": False,
                },
            }
    final_phase = reactor.phase
    final_total_mol = reactor.mass / final_phase.mean_molecular_weight * 1000.0
    final = {
        species: final_total_mol * float(final_phase[species].X[0])
        for species in phase.species_names
    }
    formed = {
        species: amount - initial.get(species, 0.0)
        for species, amount in final.items()
        if amount - initial.get(species, 0.0) > 1e-14
    }
    consumed = {
        species: initial.get(species, 0.0) - amount
        for species, amount in final.items()
        if initial.get(species, 0.0) - amount > 1e-14
    }
    net_change = {
        species: amount - initial.get(species, 0.0)
        for species, amount in final.items()
        if abs(amount - initial.get(species, 0.0)) > 1e-14
    }
    return {
        "status": "kamarinopoulou_2024_neutral_shadow_complete",
        "engine": "Cantera",
        "cantera_version": ct.__version__,
        "mechanism": KAMARINOPOULOU_MECHANISM.name,
        "mechanism_sha256": mechanism_sha256,
        "source": "Kamarinopoulou et al., Science Advances 10 (2024), adl4246",
        "doi": "10.1126/sciadv.adl4246",
        "published_species": published_species,
        "retained_heavy_reactions": phase.n_reactions,
        "inert_bath_species_added": bath_species,
        "input_species_mapping": source_mapping,
        "fixed_condition_electron_reactions_excluded": 62,
        "temperature_K": temperature_K,
        "requested_duration_s": duration_s,
        "integrated_duration_s": integrated_duration_s,
        "published_residence_time_s": KAMARINOPOULOU_PUBLISHED_RESIDENCE_TIME_S,
        "residence_time_domain_limited": duration_s > integrated_duration_s,
        "numerical_retry_used": numerical_retry_used,
        "initial_pressure_Pa_from_inventory": initial_pressure_Pa,
        "final_pressure_Pa": float(final_phase.P),
        "formed_mol": formed,
        "consumed_mol": consumed,
        "net_change_mol": net_change,
        "warnings": [
            "This literature HCN shadow shares GRI and Glarborg ancestry with the primary and Pavan branches; it is a mechanism-sensitivity result, not independent confidence.",
            "All 189 published neutral/heavy reactions are retained; the paper's 62 fixed 300 Td, ne=1e11 cm^-3 electron reactions are excluded in favor of the sandbox's self-consistent EEDF source terms.",
            "Integration is limited to the paper's 18 ms microkinetic residence-time domain; this cycle-local shadow is not propagated as an apparatus yield.",
            "GRI species absent from the paper's CHN network are present only as inert bath gas so water dilution and third-body density are retained.",
            "The published HNC thermo record conflicts with evaluated ground-state HNC enthalpy; HNC-sensitive shadow rates are not independently validated.",
            "This shadow cannot mutate the primary apparatus inventory, aqueous discovery, or benchmark score.",
        ],
        "scientific_firewall": {
            "known_miller_urey_products_accessible_to_solver": False,
            "historical_miller_urey_yields_accessible_to_solver": False,
            "paper_experimental_yields_accessible_to_solver": False,
            "mechanism_tuned_for_miller_urey": False,
            "primary_inventory_mutated_by_shadow": False,
            "quantitative_yield_claim_allowed": False,
            "hnc_ground_state_thermo_independently_validated": False,
        },
    }


def _conservative_project_nonnegative_neutrals(
    raw: np.ndarray, element_matrix: np.ndarray, species_names: list[str],
    *, maximum_negative_mol: float = 1e-12,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Remove small BDF negatives without creating their atoms.

    The signed BDF state is the conserved target. A minimum-relative-change
    correction to *already positive* species replaces simple zero clipping.
    This is a bounded numerical projection, not a chemical reaction; if the
    existing inventory cannot span the lost atoms, fail closed.
    """
    if not np.all(np.isfinite(raw)):
        raise ValueError("GRI integration returned non-finite species amounts")
    negative = sorted(
        ((species_names[i], float(value)) for i, value in enumerate(raw) if value < 0),
        key=lambda item: item[1],
    )
    if negative and -negative[0][1] > maximum_negative_mol:
        raise ValueError(
            "Time-averaged GRI integration produced negative inventory: "
            f"{negative[0][0]}={negative[0][1]:.6g} mol"
        )
    projected = np.maximum(raw, 0.0)
    correction_target = element_matrix @ (raw - projected)
    if negative:
        available = np.flatnonzero(projected > 1e-18)
        while len(available):
            scales = projected[available]
            scaled_matrix = element_matrix[:, available] * scales
            relative_change = np.linalg.lstsq(
                scaled_matrix, correction_target, rcond=1e-14,
            )[0]
            correction = scales * relative_change
            infeasible = correction < -scales
            if np.any(infeasible):
                available = available[~infeasible]
                continue
            projected[available] += correction
            break
        else:
            raise ValueError("No nonnegative conservative projection exists for GRI clipping")
    residual = element_matrix @ (projected - raw)
    projection_tolerance = max(
        5e-16, float(np.max(np.abs(correction_target), initial=0.0)) * 1e-3,
    )
    if np.any(projected < 0) or np.max(np.abs(residual), initial=0.0) > projection_tolerance:
        raise ValueError("GRI negative-species projection failed atom conservation")
    return projected, {
        "negative_species_clipped_mol": dict(negative[:10]),
        "conservative_projection_total_absolute_mol": float(
            np.sum(np.abs(projected - np.maximum(raw, 0.0)))
        ),
        "conservative_projection_adjustment_mol_by_species": {
            species_names[position]: float(value)
            for position, value in enumerate(projected - np.maximum(raw, 0.0))
            if value != 0.0
        },
        "conservative_projection_maximum_element_residual_mol": float(
            np.max(np.abs(residual), initial=0.0)
        ),
    }


def run_gri_heavy_species_cycle(
    gas_amounts_mol: np.ndarray,
    prepared: dict[str, Any],
    plasma_snapshot: dict[str, Any],
    temperature_K: float,
    gas_volume_L: float,
    duration_s: float,
    deposited_energy_J: float,
    propagated_charged_shadow: dict[str, Any] | None = None,
    source_forcing: str = "impulse",
    include_diagnostic_shadows: bool = True,
    coupled_henry_transport: dict[str, dict[str, float]] | None = None,
    boiler_species_mol: dict[str, float] | None = None,
    formaldehyde_hydration: dict[str, Any] | None = None,
    integrated_pathway_species: tuple[str, ...] = (),
    integration_max_step_fraction: float = 0.05,
    diagnostic_reaction_multipliers: tuple[tuple[int, float], ...] = (),
    diagnostic_literature_h2cn_h_to_hcn: bool = False,
    diagnostic_hebrard_h2cn_n_branches: bool = False,
    diagnostic_hebrard_h2cn_h_isomer_branches: bool = False,
    diagnostic_hebrard_hnc_h_isomerization: bool = False,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Integrate GRI chemistry with energy-capped neutral fragment sources."""
    if source_forcing not in {"impulse", "time_averaged", "inventory_limited"}:
        raise ValueError(f"Unsupported GRI source forcing: {source_forcing}")
    if not 0 < integration_max_step_fraction <= 1:
        raise ValueError("Integration maximum-step fraction must lie in (0, 1]")
    if diagnostic_literature_h2cn_h_to_hcn and not 50 <= temperature_K <= 400:
        raise ValueError("Literature H2CN + H rate is only supported from 50 to 400 K")
    if diagnostic_hebrard_h2cn_n_branches and not 290 <= temperature_K <= 310:
        raise ValueError("Room-temperature H2CN + N branch audit requires 290 to 310 K")
    if diagnostic_hebrard_h2cn_h_isomer_branches and not 290 <= temperature_K <= 310:
        raise ValueError("Room-temperature H2CN + H isomer audit requires 290 to 310 K")
    if diagnostic_hebrard_hnc_h_isomerization and not 290 <= temperature_K <= 310:
        raise ValueError("Room-temperature HNC + H isomer audit requires 290 to 310 K")
    phase = _gri_phase_with_diagnostic_rates(
        diagnostic_reaction_multipliers,
        include_literature_h2cn_h_to_hcn=diagnostic_literature_h2cn_h_to_hcn,
        include_hebrard_h2cn_n_branches=diagnostic_hebrard_h2cn_n_branches,
        include_hebrard_h2cn_h_isomer_branches=diagnostic_hebrard_h2cn_h_isomer_branches,
        include_hebrard_hnc_h_isomerization=diagnostic_hebrard_hnc_h_isomerization,
    )

    initial_gas_energy = gri_ideal_gas_internal_energy_inventory(
        gas_amounts_mol, prepared, temperature_K, diagnostic_phase=phase,
    )

    injected, source_extents, used_energy_J, energy_allocation, warnings = _neutral_fragment_extents(
        gas_amounts_mol, prepared, plasma_snapshot, deposited_energy_J
    )
    skipped_shadow = {
        "status": "deferred_until_declared_cycle_boundary",
        "formed_mol": {},
        "net_change_mol": {},
        "warnings": [
            "This non-mutating mechanism diagnostic is evaluated at the declared cycle boundary."
        ],
    }
    pavan_shadow = (
        _run_pavan_neutral_shadow(
            injected,
            prepared,
            plasma_snapshot,
            energy_allocation,
            temperature_K,
            gas_volume_L,
            duration_s,
        )
        if include_diagnostic_shadows
        else skipped_shadow
    )
    pavan_charged_shadow = propagated_charged_shadow
    if pavan_charged_shadow is None and not include_diagnostic_shadows:
        pavan_charged_shadow = skipped_shadow
    if pavan_charged_shadow is None:
        try:
            pavan_charged_shadow = _run_pavan_charged_shadow(
                gas_amounts_mol,
                prepared,
                plasma_snapshot,
                energy_allocation,
                temperature_K,
                gas_volume_L,
                duration_s,
            )
        except Exception as error:  # pragma: no cover - mechanism failures are data-dependent
            pavan_charged_shadow = {
                "status": "pavan_charged_shadow_failed",
                "mechanism": PAVAN_MECHANISM,
                "error": _cantera_error_summary(error),
                "formed_mol": {},
                "warnings": [
                    "The charge-balanced Pavan shadow failed numerically; the primary inventory was not changed."
                ],
                "scientific_firewall": {
                    "primary_inventory_mutated_by_shadow": False,
                    "quantitative_yield_claim_allowed": False,
                },
            }
    literature_hcn_shadow = (
        _run_kamarinopoulou_neutral_shadow(
            injected,
            prepared,
            temperature_K,
            gas_volume_L,
            duration_s,
        )
        if include_diagnostic_shadows
        else skipped_shadow
    )
    composition = {
        species: float(injected[prepared["index"][species]])
        for species in phase.species_names
        if injected[prepared["index"][species]] > 0
    }
    total_mol = sum(composition.values())
    if total_mol <= 0:
        return injected, {
            "status": "no_supported_heavy_species",
            "mechanism": MECHANISM,
            "source_extents": source_extents,
            "collision_energy_allocation": energy_allocation,
            "mechanism_comparison": {
                "status": "primary_unavailable_shadow_reported",
                "primary_mechanism": MECHANISM,
                "shadow": pavan_shadow,
                "charged_shadow": pavan_charged_shadow,
                "literature_hcn_shadow": literature_hcn_shadow,
                "scientific_firewall": {
                    "primary_inventory_mutated_by_shadow": False,
                    "comparison_is_independent_confidence": False,
                },
            },
            "formed_mol": {},
            "warnings": warnings,
        }

    volume_m3 = gas_volume_L * 1e-3
    mass_kg = sum(
        amount * phase.molecular_weights[phase.species_index(species)] / 1000.0
        for species, amount in composition.items()
    )
    phase.TDX = temperature_K, mass_kg / volume_m3, composition
    initial_pressure_Pa = float(phase.P)
    final = injected.copy()
    before = (
        gas_amounts_mol.copy() if source_forcing == "inventory_limited" else injected.copy()
    )
    coupled_transport_report: dict[str, Any] | None = None
    integrated_pathway_diagnostics: dict[str, Any] | None = None
    if source_forcing in {"time_averaged", "inventory_limited"} and duration_s > 0:
        starting_amounts = np.asarray(
            [gas_amounts_mol[prepared["index"][species]] for species in phase.species_names],
            dtype=float,
        )
        source_rate_mol_s = np.asarray(
            [
                (injected[prepared["index"][species]] - gas_amounts_mol[prepared["index"][species]])
                / duration_s
                for species in phase.species_names
            ],
            dtype=float,
        )
        inventory_limited = source_forcing == "inventory_limited"
        dynamic_channels = []
        if inventory_limited:
            for item in source_extents:
                target = phase.species_index(item["target"])
                initial_target = float(starting_amounts[target])
                if initial_target <= 0:
                    continue
                dynamic_channels.append({
                    "item": item,
                    "target": target,
                    "products": [phase.species_index(species) for species in item["products"]],
                    "hazard_s_inv": (
                        float(item["accepted_extent_mol"])
                        / (initial_target * duration_s)
                    ),
                })
        molecular_weights = phase.molecular_weights
        # An element absent from both the initial gas and external fragments
        # cannot be created by chemical kinetics. Freeze its species exactly
        # rather than letting BDF roundoff seed a false carbon/nitrogen product
        # in feed-removal control experiments.
        element_matrix = np.asarray(
            [
                [phase.n_atoms(species_index, element) for species_index in range(phase.n_species)]
                for element in phase.element_names
            ],
            dtype=float,
        )
        supplied_atoms = element_matrix @ (
            np.maximum(starting_amounts, 0.0)
            + np.maximum(source_rate_mol_s, 0.0) * duration_s
        )
        inactive_species = np.any(element_matrix[supplied_atoms == 0.0] > 0, axis=0)
        if np.any(starting_amounts[inactive_species] != 0.0) or np.any(source_rate_mol_s[inactive_species] != 0.0):
            raise ValueError("Inactive-element mask conflicts with supplied gas or fragments")

        def set_state(amounts: np.ndarray) -> None:
            bounded = np.maximum(amounts, 0.0)
            bounded[inactive_species] = 0.0
            state_mass_kg = float(np.dot(bounded, molecular_weights) / 1000.0)
            phase.TDX = temperature_K, state_mass_kg / volume_m3, bounded

        set_state(starting_amounts)
        initial_pressure_Pa = float(phase.P)

        tracked = [
            (species, phase.species_index(species), data)
            for species, data in (coupled_henry_transport or {}).items()
            if species in phase.species_names
        ]
        tracked_count = len(tracked)
        initial_boiler = np.asarray(
            [float((boiler_species_mol or {}).get(species, 0.0)) for species, _, _ in tracked],
            dtype=float,
        )
        hydration_state_start = phase.n_species + 2 * tracked_count
        if formaldehyde_hydration is not None:
            if not any(species == "CH2O" for species, _, _ in tracked):
                raise ValueError("GRI hydration needs tracked CH2O")
            if boiler_species_mol is None or float(boiler_species_mol.get("H2O", 0.0)) <= 0:
                raise ValueError("GRI hydration needs positive liquid-water inventory")
        initial_state = (
            np.concatenate([starting_amounts, initial_boiler, np.zeros(tracked_count)])
            if tracked_count else starting_amounts
        )
        if formaldehyde_hydration is not None:
            initial_state = np.concatenate([
                initial_state,
                [float(boiler_species_mol.get(METHANEDIOL_AQ, 0.0)), float(boiler_species_mol["H2O"])],
            ])
        if inventory_limited:
            initial_state = np.concatenate([initial_state, np.zeros(len(dynamic_channels))])

        def rhs(_time_s: float, amounts: np.ndarray) -> np.ndarray:
            gas_amounts = amounts[:phase.n_species]
            set_state(gas_amounts)
            chemistry_mol_s = (
                np.asarray(phase.net_production_rates, dtype=float) * 1000.0 * volume_m3
            )
            if inventory_limited:
                source_change = np.zeros(phase.n_species, dtype=float)
                source_channel_rates = np.zeros(len(dynamic_channels), dtype=float)
                for channel_index, channel in enumerate(dynamic_channels):
                    rate = channel["hazard_s_inv"] * max(
                        float(gas_amounts[channel["target"]]), 0.0
                    )
                    source_channel_rates[channel_index] = rate
                    source_change[channel["target"]] -= rate
                    for product in channel["products"]:
                        source_change[product] += rate
            else:
                source_change = source_rate_mol_s
            gas_change = chemistry_mol_s + source_change
            gas_change[inactive_species] = 0.0
            if not tracked_count:
                return (
                    np.concatenate([gas_change, source_channel_rates])
                    if inventory_limited else gas_change
                )
            boiler_change = np.zeros(tracked_count, dtype=float)
            captured_change = np.zeros(tracked_count, dtype=float)
            hydration_change = np.zeros(2, dtype=float)
            for offset, (species, position, data) in enumerate(tracked):
                gas_value = max(float(gas_amounts[position]), 0.0)
                liquid_value = max(float(amounts[phase.n_species + offset]), 0.0)
                if formaldehyde_hydration is not None and species == "CH2O":
                    if data["condenser_capture_rate_per_chemistry_s"] != 0:
                        raise ValueError("Closed-flask formaldehyde hydrate cannot use condenser capture")
                    changes = formaldehyde_contact_hydration_rates_mol_s(
                        gas_value, liquid_value,
                        max(float(amounts[hydration_state_start]), 0.0),
                        max(float(amounts[hydration_state_start + 1]), 0.0),
                        free_liquid_to_gas_ratio=float(formaldehyde_hydration["free_liquid_to_gas_ratio"]),
                        liquid_side_kla_per_s=float(formaldehyde_hydration["liquid_side_kla_per_s"]),
                        hydration_rate_per_s=float(formaldehyde_hydration["hydration_rate_per_s"]),
                        dehydration_rate_per_s=float(formaldehyde_hydration["dehydration_rate_per_s"]),
                    )
                    gas_change[position] += changes[0]
                    boiler_change[offset] = changes[1]
                    hydration_change[:] = changes[2:]
                    continue
                ratio = data["boiler_liquid_to_gas_equilibrium_ratio"]
                target_gas = (gas_value + liquid_value) / (1.0 + ratio)
                exchange = data["boiler_relaxation_rate_per_chemistry_s"] * (
                    gas_value - target_gas
                )
                capture = data["condenser_capture_rate_per_chemistry_s"] * gas_value
                gas_change[position] -= exchange + capture
                boiler_change[offset] = exchange
                captured_change[offset] = capture
            if inventory_limited:
                return np.concatenate([
                    gas_change, boiler_change, captured_change,
                    hydration_change if formaldehyde_hydration is not None else np.zeros(0),
                    source_channel_rates,
                ])
            return np.concatenate([
                gas_change, boiler_change, captured_change,
                hydration_change if formaldehyde_hydration is not None else np.zeros(0),
            ])

        integration_attempts = [(1e-8, 1e-18, integration_max_step_fraction)]
        if tracked_count:
            # The coupled gas/liquid state is stiffer than gas chemistry alone.
            # Retry an endpoint that violates the existing positivity bound
            # with stricter integration, never with a looser clip threshold.
            integration_attempts.append((1e-10, 1e-21, min(integration_max_step_fraction, 0.01)))
        # Captured amounts and source extents integrate rate histories but do
        # not feed back into the gas or liquid chemistry. Inactive-element gas
        # species are identically zero. Leave those Jacobian columns exactly
        # zero instead of asking SciPy to enlarge a futile perturbation.
        dynamic_columns = [
            position for position in range(phase.n_species)
            if not inactive_species[position]
        ]
        dynamic_columns.extend(range(phase.n_species, phase.n_species + tracked_count))
        if formaldehyde_hydration is not None:
            dynamic_columns.extend(range(hydration_state_start, hydration_state_start + 2))

        def bounded_jacobian(time_s: float, amounts: np.ndarray) -> np.ndarray:
            return _bounded_rate_jacobian(rhs, time_s, amounts, dynamic_columns)

        retry_reason = None
        for attempt_index, (relative_tolerance, absolute_tolerance, step_fraction) in enumerate(integration_attempts):
            integration = solve_ivp(
                rhs,
                (0.0, duration_s),
                initial_state,
                method="BDF",
                jac=bounded_jacobian,
                rtol=relative_tolerance,
                atol=absolute_tolerance,
                max_step=max(duration_s * step_fraction, 1e-9),
                dense_output=bool(integrated_pathway_species),
            )
            if not integration.success:
                retry_reason = f"BDF failed: {integration.message}"
                if attempt_index + 1 < len(integration_attempts):
                    continue
                raise ValueError(f"Time-averaged GRI integration failed: {integration.message}")
            solved_state = np.asarray(integration.y[:, -1], dtype=float)
            if (
                formaldehyde_hydration is not None
                and np.min(solved_state[hydration_state_start:hydration_state_start + 2]) < -1e-18
            ):
                retry_reason = "negative formaldehyde hydrate or water inventory"
                if attempt_index + 1 < len(integration_attempts):
                    continue
                raise ValueError(retry_reason)
            raw_amounts = solved_state[:phase.n_species].copy()
            raw_amounts[inactive_species] = 0.0
            raw_gas_conserved_residual = (
                element_matrix @ (raw_amounts - starting_amounts)
                if not tracked_count else None
            )
            try:
                solved_amounts, projection_report = _conservative_project_nonnegative_neutrals(
                    raw_amounts, element_matrix, phase.species_names,
                )
            except ValueError as error:
                retry_reason = str(error)
                if (
                    attempt_index + 1 < len(integration_attempts)
                    and "negative inventory" in retry_reason
                ):
                    continue
                raise
            break
        integrated_source_change = source_rate_mol_s * duration_s
        if inventory_limited:
            integrated_source_change = np.zeros(phase.n_species, dtype=float)
            channel_start = phase.n_species + 2 * tracked_count + (
                2 if formaldehyde_hydration is not None else 0
            )
            actual_extents = solved_state[channel_start:]
            used_energy_J = 0.0
            for channel, actual in zip(dynamic_channels, actual_extents, strict=True):
                actual = max(float(actual), 0.0)
                channel["item"]["actual_integrated_extent_mol"] = actual
                integrated_source_change[channel["target"]] -= actual
                for product in channel["products"]:
                    integrated_source_change[product] += actual
                used_energy_J += actual * float(channel["item"]["energy_cost_J_mol"])
            if used_energy_J > energy_allocation["allocated_collision_energy_J"] * (1 + 1e-7):
                raise ValueError("Inventory-limited source exceeded collision energy allocation")
            energy_allocation["resolved_fragment_energy_used_J"] = used_energy_J
            energy_allocation["nonfragment_or_inventory_limited_energy_J"] = max(
                energy_allocation["allocated_collision_energy_J"] - used_energy_J, 0.0
            )
        if integrated_pathway_species:
            # Integrate every reaction's signed progress over the solver's adaptive
            # time intervals. H2CN turnover is extremely cancellation-prone, so
            # use higher-order read-only quadrature when it is traced.
            reaction_extents = np.zeros(phase.n_reactions, dtype=float)
            quadrature_order = 5 if "H2CN" in integrated_pathway_species else 2
            abscissae, weights = np.polynomial.legendre.leggauss(quadrature_order)
            assert integration.sol is not None
            for left, right in zip(integration.t[:-1], integration.t[1:], strict=True):
                midpoint = 0.5 * (left + right)
                half_width = 0.5 * (right - left)
                for abscissa, weight in zip(abscissae, weights, strict=True):
                    amounts_at_time = np.asarray(
                        integration.sol(midpoint + half_width * abscissa),
                        dtype=float,
                    )[:phase.n_species]
                    set_state(amounts_at_time)
                    reaction_extents += (
                        half_width * weight
                        * np.asarray(phase.net_rates_of_progress, dtype=float)
                        * 1000.0
                        * volume_m3
                    )
            net_stoichiometry = phase.product_stoich_coeffs - phase.reactant_stoich_coeffs
            equations = phase.reaction_equations()
            tracked_pathways: dict[str, Any] = {}
            for species in integrated_pathway_species:
                if species not in phase.species_names:
                    continue
                species_position = phase.species_index(species)
                contributions = net_stoichiometry[species_position] * reaction_extents
                ranked = [
                    {
                        "reaction_index": int(position),
                        "equation": equations[position],
                        "species_net_mol": float(value),
                    }
                    for position, value in enumerate(contributions)
                    if abs(value) > 1e-20
                ]
                boiler_delta = 0.0
                captured = 0.0
                for offset, (tracked_species, _, _) in enumerate(tracked):
                    if tracked_species == species:
                        boiler_delta = float(
                            solved_state[phase.n_species + offset] - initial_boiler[offset]
                        )
                        captured = float(
                            solved_state[phase.n_species + tracked_count + offset]
                        )
                        if formaldehyde_hydration is not None and species == "CH2O":
                            boiler_delta += float(
                                solved_state[hydration_state_start]
                                - boiler_species_mol[METHANEDIOL_AQ]
                            )
                        break
                observed_chemistry_change = float(
                    solved_amounts[species_position]
                    - starting_amounts[species_position]
                    - integrated_source_change[species_position]
                    + boiler_delta
                    + captured
                )
                tracked_pathways[species] = {
                    "integrated_chemical_net_mol": float(np.sum(contributions)),
                    "observed_chemical_net_mol": observed_chemistry_change,
                    "quadrature_closure_residual_mol": float(
                        np.sum(contributions) - observed_chemistry_change
                    ),
                    "integrated_formation_mol": float(
                        sum(item["species_net_mol"] for item in ranked
                            if item["species_net_mol"] > 0)
                    ),
                    "integrated_consumption_mol": float(
                        -sum(item["species_net_mol"] for item in ranked
                             if item["species_net_mol"] < 0)
                    ),
                    "reaction_contributions": ranked,
                    "formation": sorted(
                        (item for item in ranked if item["species_net_mol"] > 0),
                        key=lambda item: -item["species_net_mol"],
                    )[:7],
                    "consumption": sorted(
                        (item for item in ranked if item["species_net_mol"] < 0),
                        key=lambda item: item["species_net_mol"],
                    )[:7],
                }
            integrated_pathway_diagnostics = {
                "method": f"{quadrature_order}-point Gauss quadrature on adaptive BDF intervals",
                "adaptive_interval_count": int(len(integration.t) - 1),
                "species": tracked_pathways,
                "not_a_validated_yield_or_mechanism": True,
            }
        set_state(solved_amounts)
        final_phase = phase
        for species, amount in zip(phase.species_names, solved_amounts, strict=True):
            final[prepared["index"][species]] = float(amount)
        if tracked_count:
            boiler_delta = {
                species: float(solved_state[phase.n_species + offset] - initial_boiler[offset])
                for offset, (species, _, _) in enumerate(tracked)
            }
            if formaldehyde_hydration is not None:
                boiler_delta[METHANEDIOL_AQ] = float(
                    solved_state[hydration_state_start] - boiler_species_mol[METHANEDIOL_AQ]
                )
                boiler_delta["H2O"] = float(
                    solved_state[hydration_state_start + 1] - boiler_species_mol["H2O"]
                )
            coupled_transport_report = {
                "status": "mean_field_henry_exchange_coupled_to_gri_bdf",
                "boiler_delta_mol": boiler_delta,
                "captured_mol": {
                    species: float(max(solved_state[phase.n_species + tracked_count + offset], 0.0))
                    for offset, (species, _, _) in enumerate(tracked)
                },
            }
        solver_details = {
            "method": (
                "SciPy BDF with Cantera GRI rates and inventory-limited plasma channels"
                if inventory_limited else "SciPy BDF with Cantera GRI rates"
            ),
            "relative_tolerance": relative_tolerance,
            "absolute_tolerance_mol": absolute_tolerance,
            "positivity_retry_used": attempt_index > 0,
            "initial_attempt_failure": retry_reason,
            "function_evaluations": int(integration.nfev),
            **projection_report,
            "raw_gas_maximum_element_residual_mol": (
                float(np.max(np.abs(raw_gas_conserved_residual), initial=0.0))
                if raw_gas_conserved_residual is not None else None
            ),
        }
    else:
        reactor = ct.IdealGasReactor(phase, energy="off", clone=True)
        reactor.volume = volume_m3
        network = ct.ReactorNet([reactor])
        network.rtol = 1e-9
        network.atol = 1e-18
        network.advance(duration_s)
        final_phase = reactor.phase
        final_total_mol = reactor.mass / final_phase.mean_molecular_weight * 1000.0
        for species in phase.species_names:
            final[prepared["index"][species]] = (
                final_total_mol * final_phase[species].X[0]
            )
        solver_details = {
            "method": "Cantera ReactorNet after instantaneous source",
            "relative_tolerance": 1e-9,
            "absolute_tolerance_mol": 1e-18,
        }
    delta = final - before
    formed = {
        prepared["ids"][position]: float(value)
        for position, value in enumerate(delta)
        if value > 1e-14
    }
    consumed = {
        prepared["ids"][position]: float(-value)
        for position, value in enumerate(delta)
        if value < -1e-14
    }
    primary_net_change = {
        prepared["ids"][position]: float(value)
        for position, value in enumerate(delta)
        if abs(value) > 1e-14
    }
    shadow_net_change = pavan_shadow.get("net_change_mol", {})
    charged_shadow_net_change = pavan_charged_shadow.get("net_change_mol", {})
    literature_net_change = literature_hcn_shadow.get("net_change_mol", {})
    comparison_species = sorted(
        set(primary_net_change)
        | set(shadow_net_change)
        | set(charged_shadow_net_change)
        | set(literature_net_change)
    )
    species_comparison = {}
    for species in comparison_species:
        primary_amount = float(primary_net_change.get(species, 0.0))
        shadow_amount = float(shadow_net_change.get(species, 0.0))
        charged_shadow_amount = float(charged_shadow_net_change.get(species, 0.0))
        literature_amount = float(literature_net_change.get(species, 0.0))
        species_comparison[species] = {
            "gri_net_change_mol": primary_amount,
            "pavan_shadow_net_change_mol": shadow_amount,
            "pavan_to_gri_ratio": (
                shadow_amount / primary_amount if abs(primary_amount) > 1e-30 else None
            ),
            "pavan_charged_shadow_net_change_mol": charged_shadow_amount,
            "pavan_charged_to_gri_ratio": (
                charged_shadow_amount / primary_amount
                if abs(primary_amount) > 1e-30
                else None
            ),
            "literature_hcn_shadow_net_change_mol": literature_amount,
            "literature_hcn_to_gri_ratio": (
                literature_amount / primary_amount
                if abs(primary_amount) > 1e-30
                else None
            ),
        }
    final_pathway_diagnostics: dict[str, dict[str, list[dict[str, Any]]]] = {}
    final_reaction_rates = np.asarray(final_phase.net_rates_of_progress, dtype=float)
    final_reaction_equations = final_phase.reaction_equations()
    for species in ("HCN", "CH2O"):
        if species not in final_phase.species_names:
            continue
        species_position = final_phase.species_index(species)
        contributions = []
        for reaction_position, net_rate in enumerate(final_reaction_rates):
            coefficient = (
                final_phase.product_stoich_coeff(species_position, reaction_position)
                - final_phase.reactant_stoich_coeff(species_position, reaction_position)
            )
            contribution = float(coefficient * net_rate * 1000.0 * volume_m3)
            if abs(contribution) > 1e-30:
                contributions.append({
                    "reaction_index": reaction_position,
                    "equation": final_reaction_equations[reaction_position],
                    "species_net_rate_mol_s": contribution,
                })
        final_pathway_diagnostics[species] = {
            "formation": sorted(
                (item for item in contributions if item["species_net_rate_mol_s"] > 0),
                key=lambda item: -item["species_net_rate_mol_s"],
            )[:5],
            "consumption": sorted(
                (item for item in contributions if item["species_net_rate_mol_s"] < 0),
                key=lambda item: item["species_net_rate_mol_s"],
            )[:5],
        }
    final_nonnegative = np.maximum(final, 0)
    final_gas_energy = gri_ideal_gas_internal_energy_inventory(
        final_nonnegative, prepared, temperature_K, diagnostic_phase=phase,
    )
    gas_energy_audit = {
        "basis": initial_gas_energy["basis"],
        "initial_J": initial_gas_energy["covered_internal_energy_J"],
        "final_J": final_gas_energy["covered_internal_energy_J"],
        "change_J": (
            final_gas_energy["covered_internal_energy_J"]
            - initial_gas_energy["covered_internal_energy_J"]
        ),
        "all_nonzero_gas_species_covered": (
            initial_gas_energy["all_nonzero_species_covered"]
            and final_gas_energy["all_nonzero_species_covered"]
        ),
        "includes_coupled_gas_water_transfer": coupled_henry_transport is not None,
        "complete_apparatus_energy_balance": False,
    }
    return final_nonnegative, {
        "status": "exploratory_neutral_gas_network_integrated",
        "engine": "Cantera",
        "cantera_version": ct.__version__,
        "mechanism": MECHANISM,
        "diagnostic_reaction_multipliers": {
            str(index): value for index, value in diagnostic_reaction_multipliers
        },
        "diagnostic_literature_h2cn_h_to_hcn": diagnostic_literature_h2cn_h_to_hcn,
        "diagnostic_hebrard_h2cn_n_branches": diagnostic_hebrard_h2cn_n_branches,
        "diagnostic_hebrard_h2cn_h_isomer_branches": diagnostic_hebrard_h2cn_h_isomer_branches,
        "diagnostic_hebrard_hnc_h_isomerization": diagnostic_hebrard_hnc_h_isomerization,
        "hnc_gas_water_fate_reviewed": False if diagnostic_hebrard_h2cn_h_isomer_branches else None,
        "mechanism_species": phase.n_species,
        "mechanism_reactions": phase.n_reactions,
        "temperature_K": temperature_K,
        "duration_s": duration_s,
        "source_forcing": source_forcing,
        "coupled_henry_transport": coupled_transport_report,
        "integrated_pathway_diagnostics": integrated_pathway_diagnostics,
        "final_instantaneous_pathway_diagnostics": final_pathway_diagnostics,
        "diagnostic_shadows_evaluated": include_diagnostic_shadows,
        "solver": solver_details,
        "initial_pressure_Pa_from_inventory": initial_pressure_Pa,
        "final_pressure_Pa": float(final_phase.P),
        "source_extents": source_extents,
        "source_energy_used_J": used_energy_J,
        "ideal_gas_internal_energy_audit": gas_energy_audit,
        "collision_energy_allocation": energy_allocation,
        "mechanism_comparison": {
            "status": "multi_mechanism_shadow_complete",
            "primary_mechanism": MECHANISM,
            "shadow": pavan_shadow,
            "charged_shadow": pavan_charged_shadow,
            "literature_hcn_shadow": literature_hcn_shadow,
            "species_comparison": species_comparison,
            "scientific_firewall": {
                "primary_inventory_mutated_by_shadow": False,
                "comparison_is_independent_confidence": False,
                "literature_hcn_network_is_independent_confidence": False,
                "quantitative_probability_claim_allowed": False,
            },
        },
        "formed_mol": formed,
        "consumed_mol": consumed,
        "net_change_mol": primary_net_change,
        "warnings": warnings
        + [
            "GRI-Mech is a broad combustion mechanism, not a plasma or Miller-Urey validation mechanism.",
            (
                "Energy-limited radical sources are time-averaged across the isothermal discharge residence; individual pulses are not resolved."
                if source_forcing == "time_averaged"
                else "The pulse train is lumped into one radical injection followed by an isothermal batch stage."
            ),
            "Only species in the reviewed Henry-law subset can enter the water inventory; other gas "
            "products remain gas-phase pending property evidence.",
        ],
        "scientific_firewall": {
            "known_products_accessible_to_solver": True,
            "historical_yields_accessible_to_solver": False,
            "mechanism_tuned_for_miller_urey": False,
            "quantitative_yield_claim_allowed": False,
        },
    }
