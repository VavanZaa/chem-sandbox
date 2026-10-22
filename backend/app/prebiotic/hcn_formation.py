"""Post-simulation HCN formation and fate diagnostic.

This module reads a completed apparatus result. It never supplies species,
reactions, rates, or experimental targets to the discovery calculation.
"""
from __future__ import annotations

import math
from typing import Any


_AQUEOUS = ("boiler_liquid", "condensate_trap", "collected_condensate")


def evaluate_aqueous_hcn_exposure(
    integration_steps: list[dict[str, Any]],
    operating_conditions: dict[str, Any],
) -> dict[str, Any]:
    """Read-only, step-endpoint inventory/time diagnostic, not a kinetic loss.

    An operator-split step can move HCN at any point within the step. Its
    endpoint inventory times duration is therefore a quadrature estimate,
    neither a molecular residence-time distribution nor a strict bound.
    """
    from .aqueous_speciation import load_acid_base_subset, speciate_hcn_ammonia_water

    temperature_by_place = {
        "boiler_liquid": float(operating_conditions["boiler_temperature_K"]),
        "condensate_trap": float(operating_conditions["condenser_temperature_K"]),
        "collected_condensate": float(operating_conditions["condenser_temperature_K"]),
    }
    lower, upper = load_acid_base_subset()["applicability_temperature_K"]
    totals = {
        place: {"molecular_hcn_mol_s": 0.0, "cyanide_family_mol_s": 0.0,
                "water_present_s": 0.0, "pH_min": None, "pH_max": None}
        for place in _AQUEOUS
    }
    snapshots: list[dict[str, Any]] = []
    elapsed_s = 0.0
    for step in integration_steps:
        start = float(step["start_time_s"])
        end = float(step["end_time_s"])
        duration = end - start
        if not all(math.isfinite(x) for x in (start, end)) or duration <= 0:
            raise ValueError("HCN exposure requires positive finite step durations")
        if abs(start - elapsed_s) > max(1e-8, end * 1e-10):
            raise ValueError("HCN exposure steps must be contiguous")
        elapsed_s = end
        places: dict[str, Any] = {}
        for place in _AQUEOUS:
            inventory = step["compartments_mol"].get(place, {})
            hcn = float(inventory.get("HCN", 0.0))
            cyanide = float(inventory.get("CN-", 0.0))
            water = float(inventory.get("H2O", 0.0))
            free_formaldehyde = float(inventory.get("CH2O", 0.0))
            methanediol = float(inventory.get("CH4O2_AQ", 0.0))
            if any(not math.isfinite(x) or x < -1e-14 for x in (
                hcn, cyanide, water, free_formaldehyde, methanediol,
            )):
                raise ValueError("HCN exposure contains invalid liquid inventory")
            hcn, cyanide, water = max(hcn, 0.0), max(cyanide, 0.0), max(water, 0.0)
            free_formaldehyde = max(free_formaldehyde, 0.0)
            methanediol = max(methanediol, 0.0)
            totals[place]["molecular_hcn_mol_s"] += hcn * duration
            totals[place]["cyanide_family_mol_s"] += (hcn + cyanide) * duration
            pH: float | None = None
            pKa_hcn: float | None = None
            pKw: float | None = None
            liquid_volume_L: float | None = None
            pH_status = "no_liquid_water_inventory"
            temperature = temperature_by_place[place]
            if water > 0:
                totals[place]["water_present_s"] += duration
                if lower <= temperature <= upper:
                    try:
                        view = speciate_hcn_ammonia_water(
                            hcn_mol=hcn,
                            cyanide_ion_mol=cyanide,
                            ammonia_mol=float(inventory.get("NH3", 0.0)),
                            ammonium_ion_mol=float(inventory.get("NH4+", 0.0)),
                            hydrogen_ion_mol=float(inventory.get("H+", 0.0)),
                            hydroxide_ion_mol=float(inventory.get("OH-", 0.0)),
                            carbon_dioxide_mol=float(inventory.get("CO2", 0.0)),
                            bicarbonate_ion_mol=float(inventory.get("HCO3-", 0.0)),
                            carbonate_ion_mol=float(inventory.get("CO3--", 0.0)),
                            water_mol=water,
                            temperature_K=temperature,
                        )
                        pH = float(view["pH_ideal_aqueous_model"])
                        pKa_hcn = float(view["pKa_HCN"])
                        pKw = float(view["pKw"])
                        liquid_volume_L = float(view["liquid_volume_L"])
                        pH_status = "ideal_limited_species_pH_estimated"
                        previous_min, previous_max = totals[place]["pH_min"], totals[place]["pH_max"]
                        totals[place]["pH_min"] = pH if previous_min is None else min(previous_min, pH)
                        totals[place]["pH_max"] = pH if previous_max is None else max(previous_max, pH)
                    except ValueError:
                        pH_status = "speciation_unavailable_for_inventory"
                else:
                    pH_status = "temperature_outside_evidence_range"
            places[place] = {
                "molecular_hcn_mol": hcn,
                "cyanide_ion_mol": cyanide,
                "water_mol": water,
                "aqueous_free_formaldehyde_mol": free_formaldehyde,
                "aqueous_methanediol_mol": methanediol,
                "liquid_volume_L": liquid_volume_L,
                "temperature_K": temperature,
                "pH_ideal_aqueous_model": pH,
                "pKa_HCN": pKa_hcn,
                "pKw": pKw,
                "pH_status": pH_status,
            }
        snapshots.append({
            "cycle": step["cycle"],
            "substep": step.get("substep"),
            "end_time_s": end,
            "duration_s": duration,
            "compartments": places,
        })
    return {
        "status": "step_endpoint_occupancy_diagnostic_only",
        "quadrature": "right_endpoint_inventory_times_step_duration",
        "step_count": len(snapshots),
        "total_duration_s": elapsed_s,
        "by_compartment": totals,
        "step_snapshots": snapshots,
        "prediction_inventory_mutated": False,
        "aqueous_hydrolysis_calculated": False,
        "interpretation": (
            "Inventory-time and ideal pH are post-simulation diagnostics. Operator-split "
            "within-step arrival times are unresolved; these values are not molecular "
            "residence times, hydrolysis losses, or experimental uncertainty bounds."
        ),
    }


def evaluate_hcn_formation(result: dict[str, Any]) -> dict[str, Any]:
    """Separate numerical formation from experimental yield validation."""

    ledger = result["material_ledger"]
    compartments = result["final_compartments_mol"]
    amounts = {
        place: float(species.get("HCN", 0.0))
        for place, species in compartments.items()
    }
    ionic_amounts = {
        place: float(species.get("CN-", 0.0))
        for place, species in compartments.items()
    }
    if any(not math.isfinite(value) or value < -1e-14 for value in (*amounts.values(), *ionic_amounts.values())):
        raise ValueError("HCN fate contains a non-finite or negative inventory")

    initial_molecular = float(ledger["initial_mol"].get("HCN", 0.0))
    added_molecular = float(ledger["added_mol"].get("HCN", 0.0))
    sampled_molecular = float(ledger["sampled_mol"].get("HCN", 0.0))
    initial_ionic = float(ledger["initial_mol"].get("CN-", 0.0))
    added_ionic = float(ledger["added_mol"].get("CN-", 0.0))
    sampled_ionic = float(ledger["sampled_mol"].get("CN-", 0.0))
    initial = initial_molecular + initial_ionic
    added = added_molecular + added_ionic
    sampled = sampled_molecular + sampled_ionic
    final = sum(amounts.values())
    final_ion = sum(ionic_amounts.values())
    generated = final + final_ion + sampled - initial - added
    carbon_by_species = {
        item["id"]: int(item["composition"].get("C", 0))
        for item in result["species_identity"]
    }
    feed_carbon = sum(
        carbon_by_species.get(species, 0) * (float(amount) + float(ledger["added_mol"].get(species, 0.0)))
        for species, amount in ledger["initial_mol"].items()
    ) + sum(
        carbon_by_species.get(species, 0) * float(amount)
        for species, amount in ledger["added_mol"].items()
        if species not in ledger["initial_mol"]
    )
    # This is a numerical reporting floor, not an instrumental detection limit.
    floor = max(1e-12, feed_carbon * 1e-10)
    aqueous = sum(amounts.get(place, 0.0) for place in _AQUEOUS)
    aqueous_ion = sum(ionic_amounts.get(place, 0.0) for place in _AQUEOUS)
    deposited_J = float(result["energy_ledger"]["deposited_energy_used_J"])
    ledger_closed = bool(ledger["closed"])
    if not ledger_closed:
        status = "material_ledger_open"
    elif initial + added > floor:
        status = "initial_or_added_hcn_present"
    elif generated > floor:
        status = "unseeded_hcn_formed"
    else:
        status = "no_hcn_above_numerical_floor"
    return {
        "status": status,
        "numerical_reporting_floor_mol": floor,
        "initial_hcn_mol": initial_molecular,
        "added_hcn_mol": added_molecular,
        "sampled_hcn_mol": sampled_molecular,
        "initial_cyanide_ion_mol": initial_ionic,
        "added_cyanide_ion_mol": added_ionic,
        "sampled_cyanide_ion_mol": sampled_ionic,
        "final_hcn_by_compartment_mol": amounts,
        "final_cyanide_ion_by_compartment_mol": ionic_amounts,
        "final_hcn_total_mol": final,
        "final_cyanide_ion_total_mol": final_ion,
        "final_hcn_plus_cyanide_total_mol": final + final_ion,
        "final_hcn_aqueous_mol": aqueous,
        "final_hcn_plus_cyanide_aqueous_mol": aqueous + aqueous_ion,
        "net_hcn_generated_mol": generated,
        "net_hcn_generated_basis": "HCN plus CN- conserved cyanide family",
        "net_hcn_generated_nmol_per_declared_J": (
            generated * 1e9 / deposited_J if deposited_J > 0 else None
        ),
        "feed_carbon_mol_atoms": feed_carbon,
        "net_hcn_per_feed_carbon": generated / feed_carbon if feed_carbon > 0 else None,
        "final_hcn_total_nmol_per_declared_J": final * 1e9 / deposited_J if deposited_J > 0 else None,
        "final_hcn_aqueous_nmol_per_declared_J": aqueous * 1e9 / deposited_J if deposited_J > 0 else None,
        "final_hcn_plus_cyanide_total_nmol_per_declared_J": (
            (final + final_ion) * 1e9 / deposited_J if deposited_J > 0 else None
        ),
        "final_hcn_plus_cyanide_aqueous_nmol_per_declared_J": (
            (aqueous + aqueous_ion) * 1e9 / deposited_J if deposited_J > 0 else None
        ),
        "deposited_energy_basis": result["energy_ledger"]["deposited_energy_basis"],
        "spark_mode": result["energy_ledger"]["spark_mode"],
        "material_ledger_closed": ledger_closed,
        "gas_to_water_accounted_in_final_inventory": True,
        "experimental_detection_or_yield_validated": False,
        "interpretation": (
            "A positive unseeded value is a simulated formation result, not an "
            "experimental detection, historical reproduction, or calibrated yield. "
            "Gas, liquid, collected, and wall amounts are kept distinct."
        ),
        "scientific_firewall": {
            "evaluated_after_discovery": True,
            "reference_hcn_value_used_for_simulation": False,
            "reference_hcn_value_used_for_fitting": False,
            "preexisting_gri_mechanism_contains_hcn_reactions": (
                result["energy_ledger"]["spark_mode"] == "boltzmann_eedf"
            ),
            "de_novo_product_identity_discovery_claimed": False,
            "quantitative_yield_claim_allowed": False,
        },
    }
