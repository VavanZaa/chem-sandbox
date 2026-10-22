"""Charge-balanced aqueous HCN/NH3/carbonate equilibrium and diagnostics.

The optional apparatus coupling mutates only the liquid acid-base species at
substep boundaries; neutral HCN and NH3 alone enter Henry exchange. The
post-run diagnostic reads the frozen inventory without feeding a known product
identity or yield into discovery.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any

import cantera as ct
import numpy as np
from scipy.constants import R
from scipy.optimize import brentq


DATA_PATH = (
    Path(__file__).resolve().parents[2]
    / "data" / "thermo" / "aqueous_cyanide_acid_base_v2.json"
)
WATER_KG_PER_MOL = 0.01801528


@lru_cache(maxsize=1)
def load_acid_base_subset() -> dict[str, Any]:
    return json.loads(DATA_PATH.read_text(encoding="utf-8"))


def _hcn_pka(temperature_K: float, data: dict[str, Any]) -> float:
    fit = data["hcn"]["ln_Ka_temperature_fit"]

    def ln_ka(temp: float) -> float:
        return (
            fit["constant"] + fit["linear_T"] * temp
            + fit["ln_T"] * math.log(temp) + fit["inverse_T"] / temp
        )

    reference = data["reference_temperature_K"]
    return float(data["hcn"]["preferred_pKa_298_15K"]) - (
        ln_ka(temperature_K) - ln_ka(reference)
    ) / math.log(10.0)


def _pka_from_constant_heat_capacity(
    temperature_K: float, reference: float, pka_reference: float,
    enthalpy: float, heat_capacity: float,
) -> float:
    ln_ka = -math.log(10.0) * pka_reference
    ln_ka += enthalpy / R * (1.0 / reference - 1.0 / temperature_K)
    ln_ka += heat_capacity / R * (
        math.log(temperature_K / reference) + reference / temperature_K - 1.0
    )
    return -ln_ka / math.log(10.0)


def _ammonium_pka(temperature_K: float, data: dict[str, Any]) -> float:
    item = data["ammonium"]
    return _pka_from_constant_heat_capacity(
        temperature_K, data["reference_temperature_K"], item["pKa_298_15K"],
        item["dissociation_enthalpy_J_mol"],
        item["dissociation_heat_capacity_J_mol_K"],
    )


def _carbonate_pkas(temperature_K: float, data: dict[str, Any]) -> tuple[float, float]:
    item = data["carbonate"]
    reference = data["reference_temperature_K"]
    return (
        _pka_from_constant_heat_capacity(
            temperature_K, reference, item["first_pKa_298_15K"],
            item["first_dissociation_enthalpy_J_mol"],
            item["first_dissociation_heat_capacity_J_mol_K"],
        ),
        _pka_from_constant_heat_capacity(
            temperature_K, reference, item["second_pKa_298_15K"],
            item["second_dissociation_enthalpy_J_mol"],
            item["second_dissociation_heat_capacity_J_mol_K"],
        ),
    )


def _water_pkw(temperature_K: float, density_g_cm3: float, data: dict[str, Any]) -> float:
    item = data["water"]
    a0, a1, a2 = item["alpha"]
    b0, b1, b2 = item["beta"]
    rho = density_g_cm3
    z = rho * math.exp(a0 + a1 / temperature_K + a2 * rho ** (2.0 / 3.0) / temperature_K**2)
    gas_limit = (
        0.61415 + 48251.33 / temperature_K - 67707.93 / temperature_K**2
        + 10102100.0 / temperature_K**3
    )
    return (
        -2.0 * item["coordination_number"] * (
            math.log10(1.0 + z)
            - z / (1.0 + z) * rho * (b0 + b1 / temperature_K + b2 * rho)
        )
        + gas_limit
        + 2.0 * math.log10(item["molar_mass_g_mol"] / 1000.0)
    )


def speciate_hcn_ammonia_water(
    *, hcn_mol: float, ammonia_mol: float, water_mol: float, temperature_K: float,
    cyanide_ion_mol: float = 0.0, ammonium_ion_mol: float = 0.0,
    hydrogen_ion_mol: float = 0.0, hydroxide_ion_mol: float = 0.0,
    carbon_dioxide_mol: float = 0.0, bicarbonate_ion_mol: float = 0.0,
    carbonate_ion_mol: float = 0.0,
) -> dict[str, Any]:
    """Solve ideal HCN, ammonium, carbonate and water charge neutrality.

    The reported molecular and ionic amounts are an ideal-solution equilibrium
    view of the same aqueous totals, not extra ledger inventory. The caller
    decides whether to apply them to a simulated liquid compartment.
    """
    amounts = (hcn_mol, ammonia_mol, water_mol, temperature_K,
               cyanide_ion_mol, ammonium_ion_mol, hydrogen_ion_mol, hydroxide_ion_mol,
               carbon_dioxide_mol, bicarbonate_ion_mol, carbonate_ion_mol)
    if any(not math.isfinite(value) for value in amounts):
        raise ValueError("Aqueous speciation inputs must be finite")
    solutes_and_water = (
        hcn_mol, ammonia_mol, water_mol, cyanide_ion_mol,
        ammonium_ion_mol, hydrogen_ion_mol, hydroxide_ion_mol,
        carbon_dioxide_mol, bicarbonate_ion_mol, carbonate_ion_mol,
    )
    if any(value < 0 for value in solutes_and_water) or water_mol <= 0:
        raise ValueError("Aqueous speciation requires nonnegative solutes and positive water")
    data = load_acid_base_subset()
    lower, upper = data["applicability_temperature_K"]
    if not lower <= temperature_K <= upper:
        raise ValueError("Aqueous speciation temperature is outside the evidence range")
    cyanide_total_mol = hcn_mol + cyanide_ion_mol
    ammonia_total_mol = ammonia_mol + ammonium_ion_mol
    carbonate_total_mol = carbon_dioxide_mol + bicarbonate_ion_mol + carbonate_ion_mol
    # Relative to neutral dissolved CO2, each bicarbonate or carbonate ion
    # contains one oxygen supplied by water. OH- also contains one water O.
    water_total_mol = (
        water_mol + hydroxide_ion_mol + bicarbonate_ion_mol + carbonate_ion_mol
    )
    conserved_ion_charge_mol = (
        ammonium_ion_mol + hydrogen_ion_mol - cyanide_ion_mol
        - hydroxide_ion_mol - bicarbonate_ion_mol - 2.0 * carbonate_ion_mol
    )
    water = ct.Water()
    water.TQ = temperature_K, 0.0  # saturated-liquid density, including the boiling boiler
    density_kg_m3 = float(water.density)
    liquid_volume_L = water_total_mol * WATER_KG_PER_MOL / density_kg_m3 * 1000.0
    cyanide_total_M = cyanide_total_mol / liquid_volume_L
    ammonia_total_M = ammonia_total_mol / liquid_volume_L
    carbonate_total_M = carbonate_total_mol / liquid_volume_L
    pka_hcn = _hcn_pka(temperature_K, data)
    pka_nh4 = _ammonium_pka(temperature_K, data)
    pka_carbonate_1, pka_carbonate_2 = _carbonate_pkas(temperature_K, data)
    pkw = _water_pkw(temperature_K, density_kg_m3 / 1000.0, data)
    ka_hcn = 10.0 ** -pka_hcn
    ka_nh4 = 10.0 ** -pka_nh4
    ka_carbonate_1 = 10.0 ** -pka_carbonate_1
    ka_carbonate_2 = 10.0 ** -pka_carbonate_2
    kw = 10.0 ** -pkw

    def carbonate_ions(hydrogen_M: float) -> tuple[float, float]:
        denominator = (
            hydrogen_M**2 + ka_carbonate_1 * hydrogen_M
            + ka_carbonate_1 * ka_carbonate_2
        )
        return (
            carbonate_total_M * ka_carbonate_1 * hydrogen_M / denominator,
            carbonate_total_M * ka_carbonate_1 * ka_carbonate_2 / denominator,
        )

    def charge_residual(p_h: float) -> float:
        hydrogen_M = 10.0 ** -p_h
        cyanide_M = cyanide_total_M * ka_hcn / (hydrogen_M + ka_hcn)
        ammonium_M = ammonia_total_M * hydrogen_M / (hydrogen_M + ka_nh4)
        bicarbonate_M, carbonate_M = carbonate_ions(hydrogen_M)
        return (
            hydrogen_M + ammonium_M - cyanide_M - bicarbonate_M
            - 2.0 * carbonate_M - kw / hydrogen_M
            - conserved_ion_charge_mol / liquid_volume_L
        )

    p_h = float(brentq(charge_residual, -1.0, 30.0, xtol=1e-12))
    hydrogen_M = 10.0 ** -p_h
    hydroxide_M = kw / hydrogen_M
    cyanide_M = cyanide_total_M * ka_hcn / (hydrogen_M + ka_hcn)
    ammonium_M = ammonia_total_M * hydrogen_M / (hydrogen_M + ka_nh4)
    bicarbonate_M, carbonate_M = carbonate_ions(hydrogen_M)
    cyanide_mol = cyanide_M * liquid_volume_L
    ammonium_mol = ammonium_M * liquid_volume_L
    bicarbonate_mol = bicarbonate_M * liquid_volume_L
    carbonate_mol = carbonate_M * liquid_volume_L
    hydrogen_mol = hydrogen_M * liquid_volume_L
    hydroxide_mol = hydroxide_M * liquid_volume_L
    water_after_mol = water_total_mol - hydroxide_mol - bicarbonate_mol - carbonate_mol
    if water_after_mol < 0:
        raise ValueError("Aqueous speciation would consume more water than available")
    amounts_after = {
        "HCN": cyanide_total_mol - cyanide_mol,
        "CN-": cyanide_mol,
        "NH3": ammonia_total_mol - ammonium_mol,
        "NH4+": ammonium_mol,
        "H+": hydrogen_mol,
        "OH-": hydroxide_mol,
        "H2O": water_after_mol,
        "CO2": carbonate_total_mol - bicarbonate_mol - carbonate_mol,
        "HCO3-": bicarbonate_mol,
        "CO3--": carbonate_mol,
    }
    charge_residual_mol = (
        hydrogen_mol + ammonium_mol - cyanide_mol - hydroxide_mol
        - bicarbonate_mol - 2.0 * carbonate_mol
        - conserved_ion_charge_mol
    )
    # By stoichiometry H_balance equals this charge residual. Evaluate it
    # directly to avoid cancellation against many moles of bulk water.
    hydrogen_residual_mol = charge_residual_mol
    return {
        "status": "liquid_ideal_acid_base_equilibrium_estimated",
        "dataset_id": data["dataset_id"],
        "temperature_K": temperature_K,
        "liquid_volume_L": liquid_volume_L,
        "pH_ideal_aqueous_model": p_h,
        "pKa_HCN": pka_hcn,
        "pKa_NH4": pka_nh4,
        "pKa_carbonate_1": pka_carbonate_1,
        "pKa_carbonate_2": pka_carbonate_2,
        "pKw": pkw,
        "cyanide_total_mol": cyanide_total_mol,
        "ammonia_total_mol": ammonia_total_mol,
        "carbonate_total_mol": carbonate_total_mol,
        "conserved_ion_charge_mol": conserved_ion_charge_mol,
        "equilibrated_species_mol": amounts_after,
        "cyanide_ion_fraction": cyanide_mol / cyanide_total_mol if cyanide_total_mol > 0 else None,
        "charge_balance_residual_mol": charge_residual_mol,
        "hydrogen_balance_residual_mol": hydrogen_residual_mol,
        "prediction_inventory_mutated": False,
        "limitations": data["limitations"] + [
            "This equilibrium view is not a coupled gas-water-chemistry solve; any apparatus application occurs at substep boundaries.",
            "The apparatus has no measured pH or complete dissolved-buffer/ion inventory, so this is not a calibrated cyanide-ion yield."
        ],
    }


_AQUEOUS_IONS: dict[str, tuple[str, dict[str, int], int]] = {
    "CN-": ("Cyanide ion", {"C": 1, "N": 1}, -1),
    "NH4+": ("Ammonium ion", {"N": 1, "H": 4}, 1),
    "H+": ("Hydrogen ion (formal aqueous species)", {"H": 1}, 1),
    "OH-": ("Hydroxide ion", {"O": 1, "H": 1}, -1),
    "HCO3-": ("Bicarbonate ion", {"H": 1, "C": 1, "O": 3}, -1),
    "CO3--": ("Carbonate ion", {"C": 1, "O": 3}, -2),
}
_NEUTRAL_COMPONENTS: dict[str, tuple[str, dict[str, int], int]] = {
    "HCN": ("Hydrogen cyanide", {"H": 1, "C": 1, "N": 1}, 0),
    "NH3": ("Ammonia", {"N": 1, "H": 3}, 0),
    "H2O": ("Water", {"H": 2, "O": 1}, 0),
    "CO2": ("Carbon dioxide", {"C": 1, "O": 2}, 0),
}


def augment_prepared_with_aqueous_ions(prepared: dict[str, Any]) -> None:
    """Add explicit acid-base ions to the shared atom/charge ledger."""
    original_count = len(prepared["ids"])
    for species_id, (name, composition, charge) in {
        **_NEUTRAL_COMPONENTS, **_AQUEOUS_IONS
    }.items():
        if species_id in prepared["index"]:
            identity = prepared["identities"][prepared["index"][species_id]]
            if identity["composition"] != composition or identity["charge"] != charge:
                raise ValueError(f"Declared aqueous ion {species_id} has the wrong atoms or charge")
            continue
        prepared["ids"].append(species_id)
        prepared["identities"].append({
            "id": species_id,
            "name": name,
            "composition": composition,
            "formula": species_id,
            "charge": charge,
            "isomeric_smiles": None,
            "stereo_status": "ion_formula_only",
            "stereochemistry": [],
            "phase_behavior": "dissolved_ion" if charge else "condensable",
            "enantiomer_of": None,
            "identity_note": (
                "Explicit aqueous acid-base ion; no gas-phase transport"
                if charge else "Zero-inventory neutral component of the aqueous acid-base model"
            ),
        })
        prepared["transfer"].append({
            "rates_s_inv": {"evaporation": 0.0, "condensation": 0.0, "wall_loss": 0.0},
            "evidence": {"evaporation": None, "condensation": None, "wall_loss": None},
        })
    prepared["index"] = {key: position for position, key in enumerate(prepared["ids"])}
    added = len(prepared["ids"]) - original_count
    if added:
        prepared["stoichiometry"] = np.vstack((
            prepared["stoichiometry"],
            np.zeros((added, prepared["stoichiometry"].shape[1]), dtype=float),
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


def equilibrate_aqueous_ion_inventory(
    liquid_mol: np.ndarray, prepared: dict[str, Any], temperature_K: float
) -> tuple[np.ndarray, dict[str, Any]]:
    """Equilibrate the explicit HCN/NH3/carbonate/water pool in one liquid phase."""
    index = prepared["index"]
    if liquid_mol[index["H2O"]] <= 0:
        return liquid_mol, {"status": "no_liquid_water_inventory"}
    old = liquid_mol.copy()
    report = speciate_hcn_ammonia_water(
        hcn_mol=float(old[index["HCN"]]),
        ammonia_mol=float(old[index["NH3"]]),
        water_mol=float(old[index["H2O"]]),
        cyanide_ion_mol=float(old[index["CN-"]]),
        ammonium_ion_mol=float(old[index["NH4+"]]),
        hydrogen_ion_mol=float(old[index["H+"]]),
        hydroxide_ion_mol=float(old[index["OH-"]]),
        carbon_dioxide_mol=float(old[index["CO2"]]),
        bicarbonate_ion_mol=float(old[index["HCO3-"]]),
        carbonate_ion_mol=float(old[index["CO3--"]]),
        temperature_K=temperature_K,
    )
    updated = old.copy()
    for species_id, amount in report["equilibrated_species_mol"].items():
        updated[index[species_id]] = amount
    residual = prepared["atomic"] @ (updated - old)
    expected = prepared["atomic"] @ old
    tolerances = np.maximum(np.abs(expected) * 1e-12, 1e-12)
    if np.any(np.abs(residual) > tolerances):
        raise ValueError("Aqueous acid-base update did not conserve atoms and charge")
    report["ledger_maximum_residual_mol"] = float(np.max(np.abs(residual), initial=0.0))
    report["ledger_tolerance_mol_by_conserved_quantity"] = {
        key: float(value)
        for key, value in zip(prepared["conserved"], tolerances, strict=True)
    }
    report["prediction_inventory_mutated"] = True
    return updated, report


def evaluate_aqueous_cyanide_speciation(result: dict[str, Any]) -> dict[str, Any]:
    """Speciate frozen boiler, trap and collected totals without touching them."""
    data = load_acid_base_subset()
    lower, upper = data["applicability_temperature_K"]
    compartments: dict[str, Any] = {}
    for place, temperature_key in (
        ("boiler_liquid", "boiler_temperature_K"),
        ("condensate_trap", "condenser_temperature_K"),
        ("collected_condensate", "condenser_temperature_K"),
    ):
        inventory = result["final_compartments_mol"].get(place, {})
        water_mol = float(inventory.get("H2O", 0.0))
        temperature_K = float(result["operating_conditions"][temperature_key])
        if water_mol <= 0:
            compartments[place] = {
                "status": "no_liquid_water_inventory",
                "prediction_inventory_mutated": False,
            }
        elif not lower <= temperature_K <= upper:
            compartments[place] = {
                "status": "temperature_outside_evidence_range",
                "temperature_K": temperature_K,
                "evidence_range_K": [lower, upper],
                "prediction_inventory_mutated": False,
            }
        else:
            try:
                compartments[place] = speciate_hcn_ammonia_water(
                    hcn_mol=float(inventory.get("HCN", 0.0)),
                    ammonia_mol=float(inventory.get("NH3", 0.0)),
                    water_mol=water_mol,
                    cyanide_ion_mol=float(inventory.get("CN-", 0.0)),
                    ammonium_ion_mol=float(inventory.get("NH4+", 0.0)),
                    hydrogen_ion_mol=float(inventory.get("H+", 0.0)),
                    hydroxide_ion_mol=float(inventory.get("OH-", 0.0)),
                    carbon_dioxide_mol=float(inventory.get("CO2", 0.0)),
                    bicarbonate_ion_mol=float(inventory.get("HCO3-", 0.0)),
                    carbonate_ion_mol=float(inventory.get("CO3--", 0.0)),
                    temperature_K=temperature_K,
                )
            except (ValueError, ct.CanteraError) as error:
                compartments[place] = {
                    "status": "speciation_unavailable_for_inventory",
                    "reason": str(error),
                    "prediction_inventory_mutated": False,
                }
    return {
        "status": "aqueous_speciation_diagnostic_complete",
        "dataset_id": data["dataset_id"],
        "compartments": compartments,
        "prediction_inventory_mutated": False,
        "boiler_liquid_speciated": compartments["boiler_liquid"]["status"] == "liquid_ideal_acid_base_equilibrium_estimated",
        "known_product_or_yield_used_by_simulator": False,
    }
