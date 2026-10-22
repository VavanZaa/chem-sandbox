"""Source-bounded aqueous CH2O hydration and explicit contact-state utilities.

The measured rate law is not an HCN mechanism or a mass-transfer coefficient.
An opt-in closed-flask apparatus mode carries the distinct free/hydrated
inventory and liquid-water stoichiometry through both stiff gas solvers.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
import json
import math
from pathlib import Path
from typing import Any

from scipy.constants import R
import numpy as np


DATA_PATH = (
    Path(__file__).resolve().parents[2]
    / "data" / "kinetics" / "formaldehyde_hydration_winkelman_2002_v1.json"
)
METHANEDIOL_AQ = "CH4O2_AQ"
METHANEDIOL_COMPOSITION = {"C": 1, "H": 4, "O": 2}
METHANEDIOL_IDENTITY_NOTE = "Distinct aqueous formaldehyde hydrate; no direct gas transfer"
DIAGNOSTIC_PH_LOWER_BOUND = 4.0


def augment_prepared_with_methanediol(prepared: dict[str, Any]) -> None:
    """Declare a distinct zero-inventory aqueous hydrate in the material ledger.

    This does not itself activate a reaction; the coupled gas/liquid solver
    must explicitly move CH2O and H2O into this state when enabled.
    """
    if METHANEDIOL_AQ in prepared["index"]:
        position = prepared["index"][METHANEDIOL_AQ]
        identity = prepared["identities"][position]
        transfer = prepared["transfer"][position]
        if (
            identity["composition"] != METHANEDIOL_COMPOSITION
            or identity["charge"] != 0
            or identity.get("identity_note") != METHANEDIOL_IDENTITY_NOTE
            or any(float(rate) != 0.0 for rate in transfer["rates_s_inv"].values())
        ):
            raise ValueError("Methanediol identity is reserved for the zero-transfer aqueous hydrate")
        return
    prepared["ids"].append(METHANEDIOL_AQ)
    prepared["identities"].append({
        "id": METHANEDIOL_AQ,
        "name": "Methanediol (aqueous)",
        "composition": METHANEDIOL_COMPOSITION.copy(),
        "formula": "CH4O2",
        "charge": 0,
        "isomeric_smiles": None,
        "stereo_status": "aqueous_formula_only",
        "stereochemistry": [],
        "phase_behavior": "condensable",
        "enantiomer_of": None,
        "identity_note": METHANEDIOL_IDENTITY_NOTE,
    })
    prepared["transfer"].append({
        "rates_s_inv": {"evaporation": 0.0, "condensation": 0.0, "wall_loss": 0.0},
        "evidence": {"evaporation": None, "condensation": None, "wall_loss": None},
    })
    prepared["index"] = {key: position for position, key in enumerate(prepared["ids"])}
    prepared["stoichiometry"] = np.vstack((
        prepared["stoichiometry"],
        np.zeros((1, prepared["stoichiometry"].shape[1]), dtype=float),
    ))
    inventories = []
    for identity in prepared["identities"]:
        inventory = Counter(identity["composition"])
        inventory["charge"] = identity["charge"]
        inventories.append(inventory)
    prepared["conserved"] = sorted({key for item in inventories for key in item})
    prepared["atomic"] = np.asarray([
        [inventory.get(key, 0) for inventory in inventories]
        for key in prepared["conserved"]
    ], dtype=float)


@lru_cache(maxsize=1)
def load_formaldehyde_hydration_source() -> dict[str, Any]:
    return json.loads(DATA_PATH.read_text(encoding="utf-8"))


def formaldehyde_hydration_enthalpy_J_mol() -> float:
    """Return the primary study's Van't Hoff reaction enthalpy, not heat loss.

    For ln(K_h)=B/T+c, delta_H=-R*B. This is a partial isothermal
    liquid-reaction enthalpy term, not an apparatus energy balance.
    """
    B = float(load_formaldehyde_hydration_source()["hydration_equilibrium_ratio"]["inverse_temperature_K"])
    return -R * B


def validate_diagnostic_formaldehyde_range(temperature_K: float, pH: float) -> None:
    """Bound a labeled sensitivity run; this does not extend source validity."""
    source = load_formaldehyde_hydration_source()
    if (
        not math.isfinite(temperature_K) or not math.isfinite(pH)
        or not source["temperature_K"]["minimum"] <= temperature_K <= source["temperature_K"]["maximum"]
        or not DIAGNOSTIC_PH_LOWER_BOUND <= pH <= source["pH"]["maximum"]
    ):
        raise ValueError("Diagnostic formaldehyde extrapolation is bounded to measured temperature and pH 4-7")


def formaldehyde_hydration_parameters(
    temperature_K: float, pH: float, *, allow_extrapolation: bool = False,
) -> dict[str, float | bool | str]:
    """Return measured pseudo-first-order rates, rejecting unreviewed domains."""
    if not math.isfinite(temperature_K) or not math.isfinite(pH) or temperature_K <= 0:
        raise ValueError("Formaldehyde hydration needs finite temperature and pH")
    source = load_formaldehyde_hydration_source()
    domain = (
        source["temperature_K"]["minimum"] <= temperature_K <= source["temperature_K"]["maximum"]
        and source["pH"]["minimum"] <= pH <= source["pH"]["maximum"]
    )
    if not domain and not allow_extrapolation:
        raise ValueError("Formaldehyde hydration conditions leave the measured temperature/pH domain")
    hydration = source["hydration_rate_per_s"]
    equilibrium = source["hydration_equilibrium_ratio"]
    k_h = hydration["prefactor"] * math.exp(-hydration["inverse_temperature_K"] / temperature_K)
    K_h = math.exp(
        equilibrium["inverse_temperature_K"] / temperature_K
        + equilibrium["log_intercept"]
    )
    if not math.isfinite(k_h) or not math.isfinite(K_h) or k_h <= 0 or K_h <= 0:
        raise ValueError("Formaldehyde hydration source evaluated outside numeric range")
    return {
        "source_id": source["dataset_id"],
        "temperature_K": temperature_K,
        "pH": pH,
        "within_measured_domain": domain,
        "explicit_extrapolation": not domain,
        "hydration_rate_per_s": k_h,
        "equilibrium_hydrate_to_free_ratio": K_h,
        "dehydration_rate_per_s": k_h / K_h,
    }


def free_formaldehyde_henry_from_effective(
    effective_henry_mol_kg_bar: float, hydration_equilibrium_ratio: float,
) -> float:
    """Infer a free-molecule Henry value from a *total* HCHO+hydrate value.

    This is a cross-source model-construction estimate, not a measurement.
    It assumes the reported total contains only free HCHO and methanediol.
    """
    if (
        not math.isfinite(effective_henry_mol_kg_bar)
        or not math.isfinite(hydration_equilibrium_ratio)
        or effective_henry_mol_kg_bar <= 0
        or hydration_equilibrium_ratio < 0
    ):
        raise ValueError("Free formaldehyde Henry derivation needs positive finite inputs")
    return effective_henry_mol_kg_bar / (1.0 + hydration_equilibrium_ratio)


def formaldehyde_hydration_rate_mol_s(
    free_ch2o_mol: float, methanediol_mol: float, water_mol: float,
    hydration_rate_per_s: float, dehydration_rate_per_s: float,
) -> float:
    """Net CH2O+H2O -> methanediol rate for a water-rich dilute liquid."""
    values = (free_ch2o_mol, methanediol_mol, water_mol,
              hydration_rate_per_s, dehydration_rate_per_s)
    if any(not math.isfinite(value) or value < 0 for value in values):
        raise ValueError("Formaldehyde hydration amounts and rates must be finite and nonnegative")
    if water_mol <= 0 and free_ch2o_mol > 0:
        raise ValueError("Measured pseudo-first-order hydration requires liquid water")
    return hydration_rate_per_s * free_ch2o_mol - dehydration_rate_per_s * methanediol_mol


def formaldehyde_contact_hydration_rates_mol_s(
    gas_ch2o_mol: float,
    free_ch2o_mol: float,
    methanediol_mol: float,
    water_mol: float,
    *,
    free_liquid_to_gas_ratio: float,
    liquid_side_kla_per_s: float,
    hydration_rate_per_s: float,
    dehydration_rate_per_s: float,
) -> tuple[float, float, float, float]:
    """Return gas/free/hydrate/water changes for an isolated dilute contact.

    The caller supplies a *free-molecule* Henry amount ratio and an
    independent contact rate. This function contains no gas-production
    chemistry and makes no statement about reactor-scale HCN yields.
    """
    if any(
        not math.isfinite(value) or value < 0
        for value in (gas_ch2o_mol, free_liquid_to_gas_ratio, liquid_side_kla_per_s)
    ):
        raise ValueError("Formaldehyde gas amount, ratio, and kLa must be finite and nonnegative")
    if water_mol == 0 and free_ch2o_mol >= 0 and methanediol_mol >= 0:
        # An implicit Jacobian can trial a water-depleted state. Do not
        # create hydrate without water; reverse dehydration restores it.
        hydrate = formaldehyde_hydration_rate_mol_s(
            0.0, methanediol_mol, 0.0,
            hydration_rate_per_s, dehydration_rate_per_s,
        )
    else:
        hydrate = formaldehyde_hydration_rate_mol_s(
            free_ch2o_mol, methanediol_mol, water_mol,
            hydration_rate_per_s, dehydration_rate_per_s,
        )
    transfer = liquid_side_kla_per_s * (
        free_liquid_to_gas_ratio * gas_ch2o_mol - free_ch2o_mol
    )
    return -transfer, transfer - hydrate, hydrate, -hydrate


def source_bounded_formaldehyde_contact(
    *,
    temperature_K: float,
    pH: float,
    gas_volume_L: float,
    liquid_water_mol: float,
    liquid_side_kla_per_s: float,
    allow_extrapolation: bool = False,
) -> dict[str, Any]:
    """Build *candidate* free-CH2O contact parameters from reviewed cards.

    Sander's effective total-pool Henry entry is divided by the measured
    hydration equilibrium ratio. This cross-source construction is not an
    independently measured free-CH2O constant or a calibrated kLa.
    """
    from .phase_partition import WATER_KG_PER_MOL, load_henry_subset

    if any(
        not math.isfinite(value) or value <= 0
        for value in (gas_volume_L, liquid_water_mol)
    ):
        raise ValueError("Formaldehyde contact needs positive finite gas volume and liquid water")
    if not math.isfinite(liquid_side_kla_per_s) or liquid_side_kla_per_s < 0:
        raise ValueError("Formaldehyde contact kLa must be finite and nonnegative")
    rates = formaldehyde_hydration_parameters(
        temperature_K, pH, allow_extrapolation=allow_extrapolation,
    )
    henry_data = load_henry_subset()
    entry = henry_data["entries"]["CH2O"]
    if (
        entry.get("aqueous_pool_definition")
        != "effective_total_unhydrated_formaldehyde_plus_methanediol"
        or entry.get("free_molecular_transfer_valid") is not False
    ):
        raise ValueError("Formaldehyde Henry entry is not a reviewed effective total-pool value")
    H_effective = float(entry["kH0_mol_kg_bar"]) * math.exp(
        float(entry["temperature_dependence_K"])
        * (1.0 / temperature_K - 1.0 / float(henry_data["reference_temperature_K"]))
    )
    H_free = free_formaldehyde_henry_from_effective(
        H_effective, float(rates["equilibrium_hydrate_to_free_ratio"]),
    )
    gas_molar_pressure_bar = R * temperature_K / (gas_volume_L * 1e-3) / 1e5
    K_free = H_free * gas_molar_pressure_bar * liquid_water_mol * WATER_KG_PER_MOL
    return {
        **rates,
        "effective_total_pool_henry_mol_kg_bar": H_effective,
        "derived_free_henry_mol_kg_bar": H_free,
        "free_liquid_to_gas_ratio": K_free,
        "liquid_side_kla_per_s": liquid_side_kla_per_s,
        "effective_henry_source": entry["source"],
        "cross_source_derived_not_measured": True,
        "kla_is_independently_supplied_not_fitted": True,
        "liquid_water_mass_fixed_within_substep": True,
        "not_a_calibrated_hcn_yield": True,
    }
