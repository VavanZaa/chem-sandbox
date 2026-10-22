"""Conditional, read-only aqueous HCN/formaldehyde endpoint screen.

This is deliberately outside the inventory-driving reactor. The published
aqueous equilibrium alone cannot establish reaction progress throughout a
changing spark/water cycle, especially without measured liquid pH.
"""

from __future__ import annotations

from functools import lru_cache
import json
import math
from pathlib import Path
from typing import Any

from .formaldehyde_hydration import METHANEDIOL_AQ


DATA_PATH = (
    Path(__file__).resolve().parents[2] / "data" / "kinetics"
    / "glycolonitrile_schlesinger_miller_1973_v1.json"
)


@lru_cache(maxsize=1)
def load_glycolonitrile_equilibrium_source() -> dict[str, Any]:
    return json.loads(DATA_PATH.read_text(encoding="utf-8"))


def endpoint_glycolonitrile_equilibrium(
    *, aqueous_cyanide_mol: float, aqueous_formaldehyde_pool_mol: float,
    liquid_volume_L: float, temperature_K: float, pH: float, pKa_HCN: float,
) -> dict[str, Any]:
    """Solve an isolated fixed-pH endpoint, not a time-dependent yield.

    The formaldehyde pool is *aqueous* free CH2O plus methanediol, not gas
    CH2O. The original paper distinguishes this from the kinetic free pool.
    """
    data = load_glycolonitrile_equilibrium_source()
    bounds = data["temperature_K"]
    amounts = (aqueous_cyanide_mol, aqueous_formaldehyde_pool_mol)
    if any(not math.isfinite(value) or value < 0 for value in amounts):
        raise ValueError("Endpoint inventories must be finite and nonnegative")
    if not math.isfinite(liquid_volume_L) or liquid_volume_L <= 0:
        raise ValueError("Endpoint liquid volume must be finite and positive")
    if any(not math.isfinite(value) for value in (temperature_K, pH, pKa_HCN)):
        raise ValueError("Endpoint temperature, pH and pKa must be finite")
    if not bounds["minimum"] <= temperature_K <= bounds["maximum"]:
        raise ValueError("Endpoint temperature is outside the measured range")

    k_data = data["log10_dissociation_K_M"]
    dissociation_K_M = 10.0 ** (
        k_data["intercept"] + k_data["inverse_temperature_K"] / temperature_K
    )
    # Avoid overflow if a caller probes an extreme but finite pH.
    hcn_fraction = 1.0 / (1.0 + 10.0 ** min(300.0, pH - pKa_HCN))
    cyanide, formaldehyde = amounts
    limit = min(cyanide, formaldehyde)
    if limit == 0:
        formed = 0.0
    else:
        # K = [HCN][formaldehyde pool]/[glycolonitrile]. At fixed pH,
        # [HCN] = hcn_fraction * (cyanide - formed) / volume.
        low, high = 0.0, limit
        for _ in range(100):
            formed = (low + high) / 2.0
            residual = (
                hcn_fraction * (cyanide - formed) * (formaldehyde - formed)
                - dissociation_K_M * liquid_volume_L * formed
            )
            if residual > 0:
                low = formed
            else:
                high = formed
        formed = (low + high) / 2.0
    return {
        "status": "conditional_fixed_pH_endpoint_equilibrium_screen",
        "dataset_id": data["dataset_id"],
        "dissociation_K_M": dissociation_K_M,
        "hcn_fraction_at_frozen_pH": hcn_fraction,
        "aqueous_cyanide_input_mol": cyanide,
        "aqueous_formaldehyde_pool_input_mol": formaldehyde,
        "conditional_glycolonitrile_mol": formed,
        "conditional_unbound_cyanide_mol": cyanide - formed,
        "conditional_unbound_formaldehyde_pool_mol": formaldehyde - formed,
        "fraction_of_endpoint_cyanide_potentially_bound": (
            formed / cyanide if cyanide > 0 else None
        ),
        "equilibrium_pool_basis_supported_by_original_paper_ocr": True,
        "original_pdf_visually_verified": False,
        "prediction_inventory_mutated": False,
        "is_reactor_yield_or_kinetic_prediction": False,
        "limitation": data["equilibrium_formaldehyde_basis"] + " " + data["validation_use"],
    }


def endpoint_glycolonitrile_association_timescale(
    *, aqueous_free_formaldehyde_mol: float, liquid_volume_L: float,
    temperature_K: float, pH: float, pKa_HCN: float, duration_s: float,
) -> dict[str, Any]:
    """Initial one-way association stress at a frozen endpoint composition.

    It deliberately omits reverse reaction, depletion, pH history and
    production history; it must never be read as a reactor-wide kinetic bound.
    """
    data = load_glycolonitrile_equilibrium_source()
    rate = data["association_rate_CN_minus_plus_unhydrated_formaldehyde"]
    values = (aqueous_free_formaldehyde_mol, liquid_volume_L, temperature_K,
              pH, pKa_HCN, duration_s)
    if any(not math.isfinite(value) for value in values):
        raise ValueError("Association timescale inputs must be finite")
    if aqueous_free_formaldehyde_mol < 0 or duration_s < 0 or liquid_volume_L <= 0:
        raise ValueError("Association timescale needs nonnegative amounts/time and positive volume")
    if abs(temperature_K - 298.15) > 0.01:
        raise ValueError("Association timescale uses only the measured 298.15 K coefficient")
    exponent = min(300.0, pH - pKa_HCN)
    cyanide_to_hcn = 10.0 ** exponent
    cyanide_fraction = cyanide_to_hcn / (1.0 + cyanide_to_hcn)
    free_formaldehyde_M = aqueous_free_formaldehyde_mol / liquid_volume_L
    pseudo_first_order_rate_per_s = (
        float(rate["k2_298_15K_M_inv_s"])
        * cyanide_fraction * free_formaldehyde_M
    )
    within_measured_pH = rate["pH"]["minimum"] <= pH <= rate["pH"]["maximum"]
    return {
        "status": (
            "source_pH_bounded_constant_endpoint_association_timescale"
            if within_measured_pH else
            "outside_measured_pH_conditional_endpoint_association_timescale"
        ),
        "dataset_id": data["dataset_id"],
        "within_measured_pH": within_measured_pH,
        "measured_pH_range": [rate["pH"]["minimum"], rate["pH"]["maximum"]],
        "k2_298_15K_M_inv_s": rate["k2_298_15K_M_inv_s"],
        "aqueous_free_formaldehyde_mol": aqueous_free_formaldehyde_mol,
        "aqueous_free_formaldehyde_M": free_formaldehyde_M,
        "cyanide_fraction_at_frozen_pH": cyanide_fraction,
        "pseudo_first_order_forward_rate_per_s": pseudo_first_order_rate_per_s,
        "characteristic_forward_time_s": (
            1.0 / pseudo_first_order_rate_per_s
            if pseudo_first_order_rate_per_s > 0 else None
        ),
        "constant_endpoint_one_way_conversion_fraction_over_duration": (
            -math.expm1(-pseudo_first_order_rate_per_s * duration_s)
        ),
        "duration_s": duration_s,
        "prediction_inventory_mutated": False,
        "is_reactor_yield_or_trajectory_bound": False,
        "limitation": (
            "CN- and free CH2O use final modeled amounts and frozen pH for the whole declared duration; "
            "reverse reaction and evolving species are omitted. Outside the measured pH range "
            "the reported rate is an explicit extrapolation, not source-bounded kinetics."
        ),
    }


def evaluate_glycolonitrile_step_exposure(aqueous_exposure: dict[str, Any]) -> dict[str, Any]:
    """Read-only forward-event quadrature from existing step-endpoint states.

    This deliberately does not propagate glycolonitrile or change HCN. A
    right-endpoint opportunity is neither an integrated true event count
    nor a bound because within-step arrival and reverse reaction are unknown.
    """
    totals = {
        place: {
            "source_pH_bounded_gross_forward_opportunity_mol": 0.0,
            "out_of_pH_domain_conditional_gross_forward_opportunity_mol": 0.0,
            "cyanide_bearing_steps": 0,
            "outside_measured_pH_steps": 0,
            "unavailable_steps": 0,
            "maximum_step_forward_hazard": 0.0,
        }
        for place in ("boiler_liquid", "condensate_trap", "collected_condensate")
    }
    for step in aqueous_exposure["step_snapshots"]:
        duration = float(step["duration_s"])
        if not math.isfinite(duration) or duration <= 0:
            raise ValueError("Step association screen needs finite positive duration")
        for place, summary in totals.items():
            liquid = step["compartments"][place]
            cyanide = (
                float(liquid["molecular_hcn_mol"])
                + float(liquid["cyanide_ion_mol"])
            )
            free_formaldehyde = float(liquid["aqueous_free_formaldehyde_mol"])
            if any(not math.isfinite(value) or value < 0 for value in (
                cyanide, free_formaldehyde,
            )):
                raise ValueError("Step association screen needs nonnegative finite inventory")
            if cyanide == 0:
                continue
            summary["cyanide_bearing_steps"] += 1
            if liquid["pH_status"] != "ideal_limited_species_pH_estimated":
                summary["unavailable_steps"] += 1
                continue
            try:
                estimate = endpoint_glycolonitrile_association_timescale(
                    aqueous_free_formaldehyde_mol=free_formaldehyde,
                    liquid_volume_L=float(liquid["liquid_volume_L"]),
                    temperature_K=float(liquid["temperature_K"]),
                    pH=float(liquid["pH_ideal_aqueous_model"]),
                    pKa_HCN=float(liquid["pKa_HCN"]),
                    duration_s=duration,
                )
            except (TypeError, ValueError):
                summary["unavailable_steps"] += 1
                continue
            hazard = estimate["pseudo_first_order_forward_rate_per_s"] * duration
            summary["maximum_step_forward_hazard"] = max(
                summary["maximum_step_forward_hazard"], hazard,
            )
            opportunity = (
                cyanide * estimate["constant_endpoint_one_way_conversion_fraction_over_duration"]
            )
            if estimate["within_measured_pH"]:
                summary["source_pH_bounded_gross_forward_opportunity_mol"] += opportunity
            else:
                summary["outside_measured_pH_steps"] += 1
                summary["out_of_pH_domain_conditional_gross_forward_opportunity_mol"] += opportunity
    bearing_steps = sum(item["cyanide_bearing_steps"] for item in totals.values())
    return {
        "status": (
            "no_aqueous_cyanide_to_screen" if bearing_steps == 0
            else "read_only_step_endpoint_association_opportunity_not_yield"
        ),
        "dataset_id": load_glycolonitrile_equilibrium_source()["dataset_id"],
        "by_compartment": totals,
        "cyanide_bearing_step_compartment_count": bearing_steps,
        "all_cyanide_bearing_steps_within_measured_pH": (
            all(
                item["outside_measured_pH_steps"] == 0 and item["unavailable_steps"] == 0
                for item in totals.values()
            ) if bearing_steps else None
        ),
        "prediction_inventory_mutated": False,
        "product_injected_into_solver": False,
        "is_reactor_yield_or_trajectory_bound": False,
        "interpretation": (
            "Right-endpoint one-way association opportunity with modeled pH and free CH2O; "
            "outside-pH-domain contributions are labeled extrapolations. Within-step arrivals, "
            "reaction reversibility, changes in pH and formaldehyde after association, "
            "and other sinks are omitted. This is neither a measured product nor a strict bound."
        ),
    }


def evaluate_glycolonitrile_endpoint_screen(
    result: dict[str, Any], *, duration_s: float | None = None,
) -> dict[str, Any]:
    """Screen each liquid compartment without modifying the completed run."""
    compartments: dict[str, Any] = {}
    for place in ("boiler_liquid", "condensate_trap", "collected_condensate"):
        inventory = result["final_compartments_mol"].get(place, {})
        speciation = result["aqueous_cyanide_speciation"]["compartments"][place]
        if speciation["status"] != "liquid_ideal_acid_base_equilibrium_estimated":
            compartments[place] = {
                "status": "aqueous_speciation_unavailable",
                "prediction_inventory_mutated": False,
            }
            continue
        try:
            compartments[place] = endpoint_glycolonitrile_equilibrium(
                aqueous_cyanide_mol=float(speciation["cyanide_total_mol"]),
                aqueous_formaldehyde_pool_mol=(
                    float(inventory.get("CH2O", 0.0))
                    + float(inventory.get(METHANEDIOL_AQ, 0.0))
                ),
                liquid_volume_L=float(speciation["liquid_volume_L"]),
                temperature_K=float(speciation["temperature_K"]),
                pH=float(speciation["pH_ideal_aqueous_model"]),
                pKa_HCN=float(speciation["pKa_HCN"]),
            )
            if duration_s is not None:
                try:
                    compartments[place]["association_timescale"] = endpoint_glycolonitrile_association_timescale(
                        aqueous_free_formaldehyde_mol=float(inventory.get("CH2O", 0.0)),
                        liquid_volume_L=float(speciation["liquid_volume_L"]),
                        temperature_K=float(speciation["temperature_K"]),
                        pH=float(speciation["pH_ideal_aqueous_model"]),
                        pKa_HCN=float(speciation["pKa_HCN"]),
                        duration_s=duration_s,
                    )
                except ValueError as error:
                    compartments[place]["association_timescale"] = {
                        "status": "outside_measured_rate_screen_domain",
                        "reason": str(error),
                        "prediction_inventory_mutated": False,
                    }
        except ValueError as error:
            compartments[place] = {
                "status": "outside_equilibrium_screen_domain",
                "reason": str(error),
                "prediction_inventory_mutated": False,
            }
    report = {
        "status": "conditional_endpoint_screens_complete",
        "compartments": compartments,
        "prediction_inventory_mutated": False,
        "known_product_or_yield_used_by_simulator": False,
    }
    if "hcn_formation" in result:
        report["step_exposure"] = evaluate_glycolonitrile_step_exposure(
            result["hcn_formation"]["aqueous_time_exposure"]
        )
    return report
