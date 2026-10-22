"""Couple apparatus gas inventories to the installed electron-collision solver.

The EEDF supplies conditional source opportunities. Downstream neutral and
charged branches accept only compatible, balanced channels under an explicit
energy bound; this layer does not by itself predict experimental product yield.
"""
from __future__ import annotations

from collections import Counter
import re
from typing import Any

import numpy as np
from scipy.constants import Avogadro, elementary_charge

from .boltzmann import (
    IST_CO2_O2_DATASET_ID,
    collision_set_targets,
    solve_self_consistent_eedf,
)
from .collision_intake import split_lxcat_products
from .n2d_fate import n2d_fate_diagnostic
from .schema import BoltzmannDischargeProgram, BoltzmannEEDFRequest, PrebioticApparatusRequest
from .state_fate import o1d_fate_diagnostic


_ATOM = re.compile(r"([A-Z][a-z]?)(\d*)")


def _split_products(label: str) -> list[str]:
    return split_lxcat_products(label)


def _parse_product_token(token: str) -> dict[str, Any] | None:
    raw = token.strip()
    charge = (
        1 if re.search(r"\^\+|\(\+,|(?<!\^)\+$", raw)
        else -1 if re.search(r"\^-|\(-,|(?<!\^)\-$", raw)
        else 0
    )
    # Atomic H(1s) is the electronic ground state, not an excited fragment.
    # Other parenthesized states must retain their excited-state identity until
    # an explicit supported state and relaxation/quenching kinetics exist.
    atomic_h_ground = bool(re.fullmatch(r"H\(1[sS]\)", raw))
    # X labels the lowest electronic state in these reviewed contributor
    # records; O(3P) is ground-state atomic oxygen. Do not turn O(1D) or
    # CO(a3P) into a ground-state reagent by stripping their parentheses.
    ground_state_annotation = atomic_h_ground or bool(
        re.fullmatch(r"(?:CO|O2)\(X\)|CO2\(X,v=000\)|O\(3P\)|(?:O|O2|CO2)\([+-],X\)", raw)
    )
    # Morgan H2O uses labels such as OH(7eV) for a channel's threshold, not
    # for an OH electronic state. Retain that annotation for a later threshold
    # consistency check instead of suppressing a genuine ground-state channel.
    energy_note = re.search(r"\((\d+(?:\.\d+)?)\s*eV\)$", raw, flags=re.IGNORECASE)
    excited = ("*" in raw or "(" in raw) and not ground_state_annotation and not energy_note
    formula = re.sub(r"\([^)]*\)", "", raw)
    formula = formula.replace("^+", "").replace("^-", "").replace("*", "").rstrip("+-").strip()
    matches = list(_ATOM.finditer(formula))
    if not matches or "".join(match.group(0) for match in matches) != formula:
        return None
    composition: Counter[str] = Counter()
    for match in matches:
        composition[match.group(1)] += int(match.group(2) or "1")
    return {
        "label": raw,
        "formula": formula,
        "composition": dict(composition),
        "charge": charge,
        "excited": excited,
        "atomic_h_ground_state": atomic_h_ground,
        "energy_annotation_eV": float(energy_note.group(1)) if energy_note else None,
    }


def _channel_identity(
    target_composition: dict[str, int], kind: str, product_label: str
) -> dict[str, Any]:
    tokens = _split_products(product_label)
    products = [_parse_product_token(token) for token in tokens]
    parsed = [item for item in products if item is not None]
    product_atoms: Counter[str] = Counter()
    for product in parsed:
        product_atoms.update(product["composition"])
    expected_charge = 1 if kind == "ionization" else -1 if kind == "attachment" else 0
    explicit_fragments = len(tokens) > 1
    unresolved_dissociation = "dissoc" in product_label.lower() and not explicit_fragments
    balanced = (
        len(parsed) == len(tokens)
        and dict(product_atoms) == target_composition
        and sum(item["charge"] for item in parsed) == expected_charge
    )
    if kind == "ionization":
        role = "ionization"
    elif kind == "attachment":
        role = "attachment"
    elif explicit_fragments or "dissoc" in product_label.lower():
        role = "dissociation"
    else:
        role = "excitation"
    return {
        "role": role,
        "reactive": role != "excitation",
        "products": parsed,
        "label_parsed": len(parsed) == len(tokens),
        "atoms_and_product_charge_balanced": balanced,
        "unresolved_dissociation": unresolved_dissociation,
        "electron_bookkeeping": (
            "one incident electron, two outgoing electrons"
            if kind == "ionization"
            else "incident electron captured"
            if kind == "attachment"
            else "electron is catalytic and omitted from net stoichiometry"
        ),
    }


def _uniform_drift_power_diagnostic(
    *,
    reduced_mobility: float,
    reduced_field_Td: float,
    neutral_amount_mol: float,
    volume_m3: float,
    electron_density_m3: float,
    active_duration_s: float,
    deposited_energy_J: float,
) -> dict[str, Any]:
    """Compare declared electron density with a uniform-drift power estimate.

    This is a consistency audit, not a spark circuit model: sheaths, ions,
    streamer geometry and non-electron energy deposition are not resolved.
    """

    neutral_number_density_m3 = neutral_amount_mol * Avogadro / volume_m3
    electric_field_V_m = reduced_field_Td * 1e-21 * neutral_number_density_m3
    mobility_m2_V_s = reduced_mobility / neutral_number_density_m3
    drift_power_W = (
        electron_density_m3
        * volume_m3
        * elementary_charge
        * mobility_m2_V_s
        * electric_field_V_m**2
    )
    drift_energy_J = drift_power_W * active_duration_s
    implied_density_m3 = (
        electron_density_m3 * deposited_energy_J / drift_energy_J
        if drift_energy_J > 0
        else None
    )
    return {
        "model": "uniform_two_term_electron_drift_power_audit",
        "neutral_number_density_m3": neutral_number_density_m3,
        "electric_field_V_m": electric_field_V_m,
        "electron_mobility_m2_V_s": mobility_m2_V_s,
        "declared_electron_density_m3": electron_density_m3,
        "estimated_electron_drift_power_W": drift_power_W,
        "estimated_electron_drift_energy_J": drift_energy_J,
        "deposited_energy_J": deposited_energy_J,
        "drift_to_deposited_energy_ratio": (
            drift_energy_J / deposited_energy_J if deposited_energy_J > 0 else None
        ),
        "electron_density_if_all_deposit_went_to_drift_m3": implied_density_m3,
        "circuit_and_electrode_geometry_resolved": False,
        "electron_energy_fraction_known": False,
        "used_to_scale_source_terms": False,
    }


def solve_apparatus_plasma_snapshot(
    request: PrebioticApparatusRequest,
    identities: list[dict[str, Any]],
    gas_amounts_mol: np.ndarray,
    gas_volume_L: float,
    discharge_duration_s: float,
    duty_fraction: float,
    deposited_energy_J: float,
) -> dict[str, Any]:
    """Solve one apparatus gas-loop EEDF and conditional primary source terms."""
    if request.spark.boltzmann is None:
        raise ValueError("Boltzmann apparatus coupling requires solver controls")
    return solve_mixture_plasma_snapshot(
        controls=request.spark.boltzmann,
        identities=identities,
        gas_amounts_mol=gas_amounts_mol,
        gas_volume_L=gas_volume_L,
        gas_temperature_K=request.apparatus.discharge_temperature.as_unit("K"),
        discharge_duration_s=discharge_duration_s,
        duty_fraction=duty_fraction,
        deposited_energy_J=deposited_energy_J,
    )


def solve_mixture_plasma_snapshot(
    *,
    controls: BoltzmannDischargeProgram,
    identities: list[dict[str, Any]],
    gas_amounts_mol: np.ndarray,
    gas_volume_L: float,
    gas_temperature_K: float,
    discharge_duration_s: float,
    duty_fraction: float,
    deposited_energy_J: float,
) -> dict[str, Any]:
    """Solve a fixed-volume gas EEDF without imposing apparatus geometry."""

    if min(gas_volume_L, gas_temperature_K, discharge_duration_s, deposited_energy_J) <= 0:
        raise ValueError("Plasma snapshot needs positive volume, temperature, duration and energy")
    if not 0 < duty_fraction <= 1:
        raise ValueError("Plasma duty fraction must be in (0, 1]")
    if gas_amounts_mol.shape != (len(identities),):
        raise ValueError("Gas inventory length must match species identities")
    if not np.all(np.isfinite(gas_amounts_mol)) or np.any(gas_amounts_mol < 0):
        raise ValueError("Gas inventory must be finite and nonnegative")
    supported_targets = collision_set_targets(controls.collision_set)
    by_id = {item["id"]: (position, item) for position, item in enumerate(identities)}
    supported_amounts = {
        target: float(max(gas_amounts_mol[by_id[target][0]], 0.0))
        for target in supported_targets
        if target in by_id
    }
    supported_total = sum(supported_amounts.values())
    total_material = sum(
        float(max(gas_amounts_mol[position], 0.0))
        for position, identity in enumerate(identities)
        if identity["phase_behavior"] != "electron"
    )
    neutral_material = sum(
        float(max(gas_amounts_mol[position], 0.0))
        for position, identity in enumerate(identities)
        if identity["charge"] == 0 and identity["phase_behavior"] != "electron"
    )
    if supported_total <= 0:
        return {
            "status": "no_supported_collision_targets",
            "covered_gas_fraction": 0.0,
            "mole_fractions": {},
            "source_terms": [],
            "warnings": [
                "The gas loop contains no supported collision target from "
                f"{', '.join(supported_targets)}."
            ],
        }

    fractions = {
        target: amount / supported_total
        for target, amount in supported_amounts.items()
        if amount > 0
    }
    eedf_request = BoltzmannEEDFRequest(
        solver=controls.solver,
        collision_set=controls.collision_set,
        mole_fractions=fractions,
        gas_temperature_K=gas_temperature_K,
        reduced_electric_field_Td=controls.reduced_electric_field_Td,
        energy_max_eV=controls.energy_max_eV,
        energy_cells=controls.energy_cells,
        initial_electron_temperature_eV=controls.initial_electron_temperature_eV,
        maximum_iterations=controls.maximum_iterations,
        relative_tolerance=controls.relative_tolerance,
        acknowledge_two_term_limit=True,
    )
    result = solve_self_consistent_eedf(
        eedf_request, include_same_set_cantera_reference=False
    )
    declared_electron_density_m3 = controls.electron_number_density.as_unit("1/m^3")
    volume_m3 = gas_volume_L * 1e-3
    active_duration_s = discharge_duration_s * duty_fraction
    drift_power = _uniform_drift_power_diagnostic(
        reduced_mobility=result["transport"]["reduced_electron_mobility_m_inv_V_inv_s_inv"],
        reduced_field_Td=controls.reduced_electric_field_Td,
        neutral_amount_mol=neutral_material,
        volume_m3=volume_m3,
        electron_density_m3=declared_electron_density_m3,
        active_duration_s=active_duration_s,
        deposited_energy_J=deposited_energy_J,
    )
    power_limited_density_m3 = drift_power[
        "electron_density_if_all_deposit_went_to_drift_m3"
    ]
    electron_density_m3 = min(declared_electron_density_m3, power_limited_density_m3)
    density_was_capped = electron_density_m3 < declared_electron_density_m3 * (1 - 1e-9)
    drift_power["effective_electron_density_upper_bound_m3"] = electron_density_m3
    drift_power["used_to_scale_source_terms"] = density_was_capped
    source_terms: list[dict[str, Any]] = []
    for rate in result["rate_coefficients"]:
        target = rate["target"]
        target_mol = supported_amounts.get(target, 0.0)
        target_concentration_mol_m3 = target_mol / volume_m3
        source_mol_s = (
            rate["rate_coefficient_m3_s"]
            * electron_density_m3
            * target_concentration_mol_m3
            * volume_m3
        )
        opportunity_mol = source_mol_s * active_duration_s
        minimum_energy_J = (
            opportunity_mol
            * Avogadro
            * rate["threshold_eV"]
            * elementary_charge
        )
        identity = _channel_identity(
            by_id[target][1]["composition"], rate["kind"], rate["product_label"]
        )
        source_terms.append(
            {
                **rate,
                **identity,
                "target_amount_mol": target_mol,
                "conditional_source_rate_mol_s": source_mol_s,
                "uncoupled_active_window_opportunity_mol": opportunity_mol,
                "threshold_energy_lower_bound_J": minimum_energy_J,
                "independent_threshold_energy_lower_bound_J": (
                    0.0 if rate.get("outcome_only") else minimum_energy_J
                ),
            }
        )

    source_terms.sort(
        key=lambda item: item["conditional_source_rate_mol_s"], reverse=True
    )
    reactive = [item for item in source_terms if item["reactive"]]
    lower_bound_J = sum(
        item["independent_threshold_energy_lower_bound_J"] for item in source_terms
    )
    overlapping_outcome_threshold_J = sum(
        item["threshold_energy_lower_bound_J"] for item in source_terms
        if item.get("outcome_only")
    )
    warnings = list(result["numerical_warnings"])
    covered_fraction = supported_total / total_material if total_material > 0 else 0.0
    if covered_fraction < 1 - 1e-9:
        warnings.append(
            "The imported collision set covers only "
            f"{covered_fraction:.3%} of the current gas inventory; source terms exclude other gases."
        )
    unresolved = sum(item["unresolved_dissociation"] for item in reactive)
    if unresolved:
        warnings.append(
            f"{unresolved} dissociative collision labels do not identify fragments and remain unresolved."
        )
    if lower_bound_J > deposited_energy_J * (1 + 1e-9):
        warnings.append(
            "Uncoupled collision opportunities exceed the deposited-energy ledger; they are diagnostic "
            "source rates, not simultaneous reaction extents."
        )
    if overlapping_outcome_threshold_J > 0:
        warnings.append(
            "Outcome-only N2 dissociation overlaps the Phelps electronic energy-loss manifold. "
            "Its event threshold is reported separately, not charged again to the "
            "independent collision-energy budget. Its electronic product-state branching "
            "is unresolved and must not be assumed to form two ground-state N atoms."
        )
    if density_was_capped:
        warnings.append(
            "The declared electron density and reduced field imply more uniform electron-drift "
            "energy than the deposited-energy budget. Collision source terms were capped at "
            "the upper-bound density if all deposited energy powered electron drift. Actual "
            "electron energy fraction, circuit and electrode geometry remain unresolved."
        )
    neutral_bath = {
        identity["id"]: float(gas_amounts_mol[position])
        for position, identity in enumerate(identities)
        if identity["charge"] == 0
        and identity["phase_behavior"] != "electron"
        and gas_amounts_mol[position] > 0
    }
    n2d_fate = None
    if any(
        item["target"] == "N2"
        and item.get("product_state_unresolved")
        and item["uncoupled_active_window_opportunity_mol"] > 0
        for item in source_terms
    ):
        n2d_fate = n2d_fate_diagnostic(
            neutral_bath,
            gas_volume_L=gas_volume_L,
            gas_temperature_K=gas_temperature_K,
            active_discharge_duration_s=active_duration_s,
        )
        if n2d_fate["status"] != "conditional_first_collision_fate_resolved":
            warnings.append(
                "Conditional N(2D) first-collision fate is unavailable: "
                f"{n2d_fate['reason']}. The total N2 dissociation rate remains "
                "excluded from state-resolved chemistry."
            )
    o1d_fate = None
    if controls.collision_set == IST_CO2_O2_DATASET_ID and any(
        "O(1D)" in item["product_label"]
        and item["uncoupled_active_window_opportunity_mol"] > 0
        for item in source_terms
    ):
        o1d_fate = o1d_fate_diagnostic(
            neutral_bath,
            gas_volume_L=gas_volume_L,
            gas_temperature_K=gas_temperature_K,
            active_discharge_duration_s=active_duration_s,
        )
        if o1d_fate["status"] != "conditional_fate_resolved":
            warnings.append(
                "The O(1D) first-collision fate audit is incomplete: "
                f"{o1d_fate['reason']}. Excited products remain excluded from "
                "ground-state chemistry."
            )
    return {
        "status": "apparatus_mixture_eedf_solved",
        "mole_fractions": result["conditions"]["mole_fractions"],
        "covered_gas_fraction": covered_fraction,
        "electron_number_density_m3": electron_density_m3,
        "declared_electron_number_density_m3": declared_electron_density_m3,
        "active_discharge_duration_s": active_duration_s,
        "mean_electron_energy_eV": result["eedf"]["mean_energy_eV"],
        "electron_temperature_K": result["eedf"]["electron_temperature_K"],
        "solver": result["solver"],
        "collision_set": {
            "dataset_id": result["collision_set"]["dataset_id"],
            "sha256": result["collision_set"]["sha256"],
            "source": result["collision_set"]["source"],
            "selected_process_count": result["collision_set"]["selected_process_count"],
        },
        "source_terms": source_terms,
        "n2d_fate_diagnostic": n2d_fate,
        "o1d_fate_diagnostic": o1d_fate,
        "reactive_channel_count": len(reactive),
        "balanced_reactive_channel_count": sum(
            item["atoms_and_product_charge_balanced"]
            and not item["unresolved_dissociation"]
            for item in reactive
        ),
        "threshold_energy_lower_bound_J": lower_bound_J,
        "overlapping_outcome_threshold_energy_J": overlapping_outcome_threshold_J,
        "deposited_energy_J": deposited_energy_J,
        "drift_power_diagnostic": drift_power,
        "warnings": warnings,
        "scientific_firewall": {
            "known_products_accessible_to_solver": False,
            "reaction_recipe_used": False,
            "heavy_species_inventory_mutated": False,
            "quantitative_yield_claim_allowed": False,
        },
    }
