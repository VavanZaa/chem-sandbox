"""Evidence-labelled gas/water partitioning for generated prebiotic feedstocks."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.constants import R

import cantera as ct


DATA_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "thermo"
    / "prebiotic_henry_water.json"
)
WATER_KG_PER_MOL = 0.01801528


def load_henry_subset() -> dict[str, Any]:
    return json.loads(DATA_PATH.read_text(encoding="utf-8"))


def enforce_water_vapor_capacity(
    gas_mol: np.ndarray,
    water_phase_mol: np.ndarray,
    prepared: dict[str, Any],
    temperature_K: float,
    gas_volume_L: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Cap water vapor at pure-water saturation in a fixed gas volume.

    The former first-order evaporation step could place many liquid-water moles
    into a small gas loop. This equilibrium cap returns supersaturated water to
    the contacting liquid before plasma and heavy-gas kinetics are evaluated.
    Non-condensable inventory is never removed, so this is not a pressure
    controller and any remaining total-pressure mismatch stays explicit.
    """

    index = prepared["index"]
    if "H2O" not in index:
        return gas_mol, water_phase_mol, {
            "status": "no_water_species",
            "recondensed_water_mol": 0.0,
        }
    if temperature_K >= ct.Water().critical_temperature:
        return gas_mol, water_phase_mol, {
            "status": "above_water_critical_temperature",
            "temperature_K": temperature_K,
            "recondensed_water_mol": 0.0,
        }
    water = ct.Water()
    water.TQ = temperature_K, 0.0
    saturation_pressure_Pa = float(water.P)
    gas_volume_m3 = gas_volume_L * 1e-3
    maximum_water_vapor_mol = (
        saturation_pressure_Pa * gas_volume_m3 / (R * temperature_K)
    )
    position = index["H2O"]
    gas = gas_mol.copy()
    liquid = water_phase_mol.copy()
    recondensed = max(float(gas[position]) - maximum_water_vapor_mol, 0.0)
    recondensed = min(recondensed, float(gas[position]))
    gas[position] -= recondensed
    liquid[position] += recondensed
    water_partial_pressure_Pa = float(gas[position]) * R * temperature_K / gas_volume_m3
    ideal_total_pressure_Pa = float(np.maximum(gas, 0).sum()) * R * temperature_K / gas_volume_m3
    return gas, liquid, {
        "status": "water_vapor_capacity_applied",
        "temperature_K": temperature_K,
        "gas_volume_L": gas_volume_L,
        "saturation_pressure_Pa": saturation_pressure_Pa,
        "maximum_water_vapor_mol": maximum_water_vapor_mol,
        "water_partial_pressure_Pa": water_partial_pressure_Pa,
        "ideal_total_pressure_Pa": ideal_total_pressure_Pa,
        "recondensed_water_mol": recondensed,
        "method": "Cantera Water saturation state with ideal-gas vapor capacity",
        "limitations": [
            "Pure-water saturation is used; dissolved gases, droplets, pressure drop and nonideal vapor fugacity are not solved.",
            "The cap removes only supersaturated water and does not force total pressure to the declared apparatus pressure.",
        ],
    }


def equilibrate_gas_water(
    gas_mol: np.ndarray,
    water_phase_mol: np.ndarray,
    prepared: dict[str, Any],
    temperature_K: float,
    pressure_Pa: float,
    gas_volume_L: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Apply dilute Henry partitioning in the fixed-volume gas loop.

    The gas partial pressure is nRT/V, not a stale total-pressure estimate
    from before condensation changed the headspace inventory.
    """

    if min(temperature_K, pressure_Pa, gas_volume_L) <= 0:
        raise ValueError("Henry partition requires positive temperature, pressure and gas volume")

    dataset = load_henry_subset()
    index = prepared["index"]
    if "H2O" not in index:
        return gas_mol, water_phase_mol, {
            "status": "no_water_species",
            "transferred_to_water_mol": {},
            "entries_used": [],
        }
    water_mass_kg = float(max(water_phase_mol[index["H2O"]], 0.0)) * WATER_KG_PER_MOL
    gas_total_mol = float(np.maximum(gas_mol, 0).sum())
    if water_mass_kg <= 0 or gas_total_mol <= 0:
        return gas_mol, water_phase_mol, {
            "status": "no_contact_inventory",
            "transferred_to_water_mol": {},
            "entries_used": [],
        }

    gas = gas_mol.copy()
    liquid = water_phase_mol.copy()
    transferred: dict[str, float] = {}
    used = []
    gas_molar_pressure_bar = R * temperature_K / (gas_volume_L * 1e-3) / 1e5
    for species_id, entry in dataset["entries"].items():
        if species_id not in index:
            continue
        position = index[species_id]
        combined = float(max(gas[position], 0.0) + max(liquid[position], 0.0))
        if combined <= 0:
            continue
        kH = float(entry["kH0_mol_kg_bar"]) * math.exp(
            float(entry["temperature_dependence_K"])
            * (1.0 / temperature_K - 1.0 / dataset["reference_temperature_K"])
        )
        liquid_to_gas_ratio = kH * gas_molar_pressure_bar * water_mass_kg
        target_liquid = combined * liquid_to_gas_ratio / (1.0 + liquid_to_gas_ratio)
        movement = target_liquid - liquid[position]
        movement = min(movement, gas[position]) if movement >= 0 else -min(-movement, liquid[position])
        gas[position] -= movement
        liquid[position] += movement
        transferred[species_id] = float(movement)
        used.append(
            {
                "species_id": species_id,
                "kH_mol_kg_bar": kH,
                "water_mass_kg": water_mass_kg,
                "source": entry["source"],
                "selection": entry["selection"],
            }
        )
    return gas, liquid, {
        "status": "henry_partition_applied",
        "temperature_K": temperature_K,
        "pressure_Pa": pressure_Pa,
        "actual_ideal_gas_pressure_Pa": gas_total_mol * gas_molar_pressure_bar * 1e5,
        "gas_volume_L": gas_volume_L,
        "partial_pressure_basis": "fixed_volume_ideal_gas_nRT_over_V",
        "dataset_id": dataset["dataset_id"],
        "transferred_to_water_mol": transferred,
        "entries_used": used,
        "limitations": [
            "Dilute ideal Henry partitioning is used; aqueous activity, pH, hydration and salting effects are not solved.",
            "The supplied pressure is diagnostic only; species partial pressures use the current fixed-volume gas inventory.",
            "Only species with reviewed entries in the local property subset are repartitioned.",
        ],
    }


def transfer_gas_water_finite_rate(
    gas_mol: np.ndarray,
    water_phase_mol: np.ndarray,
    prepared: dict[str, Any],
    temperature_K: float,
    pressure_Pa: float,
    gas_volume_L: float,
    duration_s: float,
    liquid_side_kla_per_s: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Exact two-compartment Henry relaxation for a declared liquid-side kLa.

    For each neutral species, dn_liquid/dt = kLa * (K*n_gas - n_liquid),
    where K is the dimensionless equilibrium liquid/gas amount ratio for
    this fixed gas volume and water inventory. The analytic update conserves
    each species exactly. kLa is an apparatus input, not inferred from Henry
    data or an HCN measurement.
    """

    if min(temperature_K, pressure_Pa, gas_volume_L) <= 0:
        raise ValueError("Finite-rate Henry transfer requires positive temperature, pressure and gas volume")
    if duration_s < 0 or liquid_side_kla_per_s < 0:
        raise ValueError("Finite-rate Henry transfer requires nonnegative duration and kLa")
    index = prepared["index"]
    if "H2O" not in index:
        return gas_mol, water_phase_mol, {
            "status": "no_water_species", "transferred_to_water_mol": {}, "entries_used": [],
        }
    water_mass_kg = float(max(water_phase_mol[index["H2O"]], 0.0)) * WATER_KG_PER_MOL
    if water_mass_kg <= 0:
        return gas_mol, water_phase_mol, {
            "status": "no_contact_water", "transferred_to_water_mol": {}, "entries_used": [],
        }

    dataset = load_henry_subset()
    gas = gas_mol.copy()
    liquid = water_phase_mol.copy()
    gas_molar_pressure_bar = R * temperature_K / (gas_volume_L * 1e-3) / 1e5
    transferred: dict[str, float] = {}
    used: list[dict[str, Any]] = []
    for species_id, entry in dataset["entries"].items():
        if species_id not in index:
            continue
        position = index[species_id]
        total = float(gas[position] + liquid[position])
        if total <= 0:
            continue
        kH = float(entry["kH0_mol_kg_bar"]) * math.exp(
            float(entry["temperature_dependence_K"])
            * (1.0 / temperature_K - 1.0 / dataset["reference_temperature_K"])
        )
        ratio = kH * gas_molar_pressure_bar * water_mass_kg
        target_liquid = total * ratio / (1.0 + ratio)
        relaxation_fraction = -math.expm1(-liquid_side_kla_per_s * (1.0 + ratio) * duration_s)
        movement = (target_liquid - float(liquid[position])) * relaxation_fraction
        movement = min(movement, float(gas[position])) if movement >= 0 else -min(-movement, float(liquid[position]))
        gas[position] -= movement
        liquid[position] += movement
        transferred[species_id] = float(movement)
        used.append({
            "species_id": species_id,
            "kH_mol_kg_bar": kH,
            "liquid_to_gas_equilibrium_ratio": ratio,
            "relaxation_fraction": relaxation_fraction,
            "source": entry["source"],
            "selection": entry["selection"],
        })
    return gas, liquid, {
        "status": "finite_rate_henry_contact_applied",
        "method": "analytic two-compartment liquid-side kLa relaxation",
        "temperature_K": temperature_K,
        "pressure_Pa": pressure_Pa,
        "gas_volume_L": gas_volume_L,
        "duration_s": duration_s,
        "liquid_side_kla_per_s": liquid_side_kla_per_s,
        "water_mass_kg": water_mass_kg,
        "dataset_id": dataset["dataset_id"],
        "transferred_to_water_mol": transferred,
        "entries_used": used,
        "limitations": [
            "The user-declared kLa is not inferred from this apparatus or from an HCN yield.",
            "Each substep applies gas chemistry before transfer; gas/transfer operator-splitting convergence is required.",
            "Only neutral species in the reviewed Henry subset transfer; aqueous ionization is updated separately.",
        ],
    }


def closed_flask_coupled_henry_coefficients(
    prepared: dict[str, Any],
    gas_volume_L: float,
    temperature_K: float,
    water_phase_mol: float,
    liquid_side_kla_per_s: float,
) -> dict[str, dict[str, float]]:
    """Map declared closed-flask kLa onto the existing coupled BDF balance.

    The BDF transport RHS uses ``rate * (n_g - (n_g+n_l)/(1+K))``.
    Setting ``rate = kLa*(1+K)`` makes it exactly
    ``d n_l / dt = kLa*(K*n_g - n_l)`` without gas/transfer splitting.
    The coefficient is supplied, not estimated from a product amount.
    """

    if gas_volume_L <= 0 or temperature_K <= 0:
        raise ValueError("Coupled closed-flask contact requires positive volume and temperature")
    if water_phase_mol < 0 or liquid_side_kla_per_s < 0:
        raise ValueError("Coupled closed-flask contact requires nonnegative water and kLa")
    if water_phase_mol == 0 or liquid_side_kla_per_s == 0:
        return {}
    dataset = load_henry_subset()
    water_mass_kg = water_phase_mol * WATER_KG_PER_MOL
    gas_molar_pressure_bar = R * temperature_K / (gas_volume_L * 1e-3) / 1e5
    output: dict[str, dict[str, float]] = {}
    for species_id, entry in dataset["entries"].items():
        if species_id not in prepared["index"]:
            continue
        kH = float(entry["kH0_mol_kg_bar"]) * math.exp(
            float(entry["temperature_dependence_K"])
            * (1.0 / temperature_K - 1.0 / dataset["reference_temperature_K"])
        )
        ratio = kH * gas_molar_pressure_bar * water_mass_kg
        output[species_id] = {
            "boiler_liquid_to_gas_equilibrium_ratio": ratio,
            "boiler_relaxation_rate_per_chemistry_s": liquid_side_kla_per_s * (1.0 + ratio),
            "condenser_capture_rate_per_chemistry_s": 0.0,
            "condenser_contact_water_kg": 0.0,
        }
    return output


def capture_gas_water_by_flow(
    gas_mol: np.ndarray,
    contact_liquid_mol: np.ndarray,
    prepared: dict[str, Any],
    temperature_K: float,
    pressure_Pa: float,
    gas_volume_L: float,
    circulation_rate_L_s: float,
    duration_s: float,
    contact_water_holdup_L: float,
    available_water_mol: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Upper-bound continuous capture for a well-mixed recirculating gas loop.

    Each loop turnover is assumed to approach Henry equilibrium with a fixed
    water contact hold-up. The gas therefore loses a species at Q/V times the
    single-pass equilibrium capture fraction, integrated analytically over
    the declared duration. The hold-up is an explicit apparatus scenario, not
    created water inventory, and is limited by water actually present.
    """

    if min(gas_volume_L, contact_water_holdup_L, temperature_K) <= 0:
        raise ValueError("Flow-limited capture requires positive geometry and temperature")
    if circulation_rate_L_s < 0 or duration_s < 0 or available_water_mol < 0:
        raise ValueError("Flow, duration and available water cannot be negative")
    dataset = load_henry_subset()
    gas = gas_mol.copy()
    liquid = contact_liquid_mol.copy()
    if circulation_rate_L_s == 0 or duration_s == 0 or available_water_mol == 0:
        return gas, liquid, {
            "status": "no_flow_or_contact_water",
            "method": "well-mixed flow-limited Henry capture upper bound",
            "transferred_to_water_mol": {},
            "contact_water_mass_kg": 0.0,
            "gas_turnovers": circulation_rate_L_s * duration_s / gas_volume_L,
            "entries_used": [],
        }
    water = ct.Water()
    if temperature_K >= water.critical_temperature:
        raise ValueError("Flow-limited aqueous capture requires subcritical condenser water")
    water.TQ = temperature_K, 0.0
    declared_water_mass_kg = contact_water_holdup_L * 1e-3 * float(water.density)
    contact_water_mass_kg = min(declared_water_mass_kg, available_water_mol * WATER_KG_PER_MOL)
    gas_total_mol = float(np.maximum(gas, 0).sum())
    if gas_total_mol <= 0 or contact_water_mass_kg <= 0:
        return gas, liquid, {
            "status": "no_contact_inventory",
            "method": "well-mixed flow-limited Henry capture upper bound",
            "transferred_to_water_mol": {},
            "contact_water_mass_kg": contact_water_mass_kg,
            "gas_turnovers": circulation_rate_L_s * duration_s / gas_volume_L,
            "entries_used": [],
        }
    turnovers = circulation_rate_L_s * duration_s / gas_volume_L
    gas_molar_pressure_bar = R * temperature_K / (gas_volume_L * 1e-3) / 1e5
    transferred: dict[str, float] = {}
    used: list[dict[str, Any]] = []
    for species_id, entry in dataset["entries"].items():
        if species_id not in prepared["index"]:
            continue
        position = prepared["index"][species_id]
        if gas[position] <= 0:
            continue
        kH = float(entry["kH0_mol_kg_bar"]) * math.exp(
            float(entry["temperature_dependence_K"])
            * (1.0 / temperature_K - 1.0 / dataset["reference_temperature_K"])
        )
        equilibrium_ratio = kH * gas_molar_pressure_bar * contact_water_mass_kg
        single_pass_capture_fraction = equilibrium_ratio / (1.0 + equilibrium_ratio)
        capture_fraction = -math.expm1(-turnovers * single_pass_capture_fraction)
        moved = float(gas[position]) * capture_fraction
        gas[position] -= moved
        liquid[position] += moved
        transferred[species_id] = moved
        used.append({
            "species_id": species_id,
            "kH_mol_kg_bar": kH,
            "single_pass_capture_fraction": single_pass_capture_fraction,
            "source": entry["source"],
        })
    return gas, liquid, {
        "status": "flow_limited_henry_capture_applied",
        "method": "well-mixed gas-loop turnover with per-pass Henry equilibrium upper bound",
        "dataset_id": dataset["dataset_id"],
        "contact_water_holdup_L": contact_water_holdup_L,
        "contact_water_mass_kg": contact_water_mass_kg,
        "actual_ideal_gas_pressure_Pa": gas_total_mol * gas_molar_pressure_bar * 1e5,
        "supplied_pressure_Pa": pressure_Pa,
        "partial_pressure_basis": "fixed_volume_ideal_gas_nRT_over_V",
        "gas_turnovers": turnovers,
        "transferred_to_water_mol": transferred,
        "entries_used": used,
        "limitations": [
            "Perfect per-pass approach to Henry equilibrium is an upper bound; interfacial area, residence distribution and mass-transfer resistance are not measured.",
            "The declared contact hold-up is a scenario value and does not add water to the material ledger.",
            "The supplied pressure is diagnostic only; species partial pressures use the current fixed-volume gas inventory.",
            "pH, activities, hydration and reactions inside the trap are not solved.",
        ],
    }


def coupled_henry_transport_coefficients(
    prepared: dict[str, Any],
    gas_volume_L: float,
    circulation_rate_L_s: float,
    boiler_temperature_K: float,
    condenser_temperature_K: float,
    boiler_water_mol: float,
    active_water_mol: float,
    condenser_contact_water_holdup_L: float,
    boiler_stage_s: float,
    condenser_stage_s: float,
    chemistry_residence_s: float,
) -> dict[str, dict[str, float]]:
    """Mean-field gas/liquid rates on the spark-residence integration clock.

    These rates couple the two water contacts to the gas chemistry within one
    stiff integration. The correspondence between a physical stage second and
    a spark-zone residence second is explicit. It is a well-mixed overlapping
    apparatus approximation, not a spatially resolved condenser model.
    """

    if min(gas_volume_L, boiler_temperature_K, condenser_temperature_K) <= 0:
        raise ValueError("Coupled Henry transport requires positive geometry and temperature")
    if min(circulation_rate_L_s, boiler_water_mol, active_water_mol,
           condenser_contact_water_holdup_L, boiler_stage_s, condenser_stage_s,
           chemistry_residence_s) < 0:
        raise ValueError("Coupled Henry transport inputs cannot be negative")
    if chemistry_residence_s == 0 or circulation_rate_L_s == 0:
        return {}
    water = ct.Water()
    if condenser_temperature_K >= water.critical_temperature:
        raise ValueError("Coupled aqueous capture requires subcritical condenser water")
    water.TQ = condenser_temperature_K, 0.0
    condenser_water_kg = min(
        condenser_contact_water_holdup_L * 1e-3 * float(water.density),
        active_water_mol * WATER_KG_PER_MOL,
    )
    boiler_water_kg = boiler_water_mol * WATER_KG_PER_MOL
    loop_turnover_s_inv = circulation_rate_L_s / gas_volume_L
    gas_volume_m3 = gas_volume_L * 1e-3
    dataset = load_henry_subset()
    output: dict[str, dict[str, float]] = {}
    for species_id, entry in dataset["entries"].items():
        if species_id not in prepared["index"]:
            continue
        def ratio(temperature_K: float, water_kg: float) -> float:
            kH = float(entry["kH0_mol_kg_bar"]) * math.exp(
                float(entry["temperature_dependence_K"])
                * (1.0 / temperature_K - 1.0 / dataset["reference_temperature_K"])
            )
            return kH * (R * temperature_K / gas_volume_m3 / 1e5) * water_kg

        boiler_ratio = ratio(boiler_temperature_K, boiler_water_kg)
        condenser_ratio = ratio(condenser_temperature_K, condenser_water_kg)
        output[species_id] = {
            "boiler_liquid_to_gas_equilibrium_ratio": boiler_ratio,
            "boiler_relaxation_rate_per_chemistry_s": (
                loop_turnover_s_inv * boiler_stage_s / chemistry_residence_s
            ),
            "condenser_capture_rate_per_chemistry_s": (
                loop_turnover_s_inv * condenser_stage_s / chemistry_residence_s
                * condenser_ratio / (1.0 + condenser_ratio)
            ),
            "condenser_contact_water_kg": condenser_water_kg,
        }
    return output
