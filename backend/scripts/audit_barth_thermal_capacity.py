"""Conditional water-heating screen for the assay-free Barth feed card.

This is not an apparatus heat balance: actual discharge deposition, glass
heat capacity, gas/liquid reaction enthalpy, evaporation, and wall losses
are not independently measured for the selected cyanide-assayed run.
"""

from __future__ import annotations

import json
import math
from typing import Any

import cantera as ct

from scripts.probe_barth_feed_pathway import FEED_CARD, make_probe_payload


def audit_barth_water_thermal_capacity(*, duration_min: float | None = None) -> dict[str, Any]:
    feed = json.loads(FEED_CARD.read_text(encoding="utf-8"))
    payload = make_probe_payload(duration_s=60.0)
    spark = payload["spark"]
    volume_L = float(feed["reported_water_volume_L"])
    pressure_Pa = float(feed["diagnostic_only_assumptions"]["initial_pressure_Pa"])
    temperature_K = float(feed["diagnostic_only_assumptions"]["gas_and_water_temperature_K"])
    reported_duration_min = float(feed["reported_spark_duration_min"])
    selected_duration_min = reported_duration_min if duration_min is None else duration_min
    if not math.isfinite(selected_duration_min) or not 0 < selected_duration_min <= reported_duration_min:
        raise ValueError("Barth thermal screen needs a positive duration within the reported run")
    duration_s = selected_duration_min * 60.0
    duty = float(spark["pulse_width"]["value"]) * float(spark["frequency"]["value"])
    if not all(math.isfinite(value) and value > 0 for value in (
        volume_L, pressure_Pa, temperature_K, duration_s,
    )) or not 0 < duty <= 1:
        raise ValueError("Barth thermal screen requires finite physical feed and duty inputs")
    nominal_power_W = (
        float(spark["voltage"]["value"])
        * float(spark["current"]["value"])
        * duty
        * float(spark["coupling_efficiency"])
    )
    nominal_energy_J = nominal_power_W * duration_s

    water = ct.Water()
    water.TP = temperature_K, pressure_Pa
    if water.phase_of_matter != "liquid":
        raise ValueError("Barth thermal screen starts outside the liquid-water phase")
    initial_density_kg_m3 = float(water.density)
    initial_cp_J_kg_K = float(water.cp_mass)
    initial_h_J_kg = float(water.enthalpy_mass)
    mass_kg = volume_L * 1e-3 * initial_density_kg_m3
    water.PQ = pressure_Pa, 0.0
    saturation_temperature_K = float(water.T)
    sensible_J = mass_kg * (float(water.enthalpy_mass) - initial_h_J_kg)
    if saturation_temperature_K <= temperature_K or sensible_J <= 0:
        raise ValueError("Barth thermal screen needs subcooled water and positive sensible heat")

    return {
        "status": "conditional_water_thermal_capacity_screen_not_heat_balance",
        "feed_card": feed["id"],
        "feed_card_contains_hcn_assay": False,
        "water_property_basis": (
            "Cantera 3.2 Water pure-fluid equation of state, IAPWS-95; "
            "https://www.iapws.org/relguide/IAPWS-95.html"
        ),
        "nominal_discharge_energy_basis": "hypothetical declared electrical estimate, not measured deposited energy",
        "initial_water_volume_L": volume_L,
        "initial_water_mass_kg": mass_kg,
        "initial_pressure_Pa": pressure_Pa,
        "initial_temperature_K": temperature_K,
        "initial_liquid_cp_J_kg_K": initial_cp_J_kg_K,
        "initial_water_heat_capacity_J_K": mass_kg * initial_cp_J_kg_K,
        "saturation_temperature_at_initial_pressure_K": saturation_temperature_K,
        "liquid_sensible_enthalpy_to_saturation_J": sensible_J,
        "declared_spark_duty_fraction": duty,
        "declared_spark_coupling_efficiency": float(spark["coupling_efficiency"]),
        "declared_average_deposited_power_W": nominal_power_W,
        "declared_duration_s": duration_s,
        "declared_total_deposited_energy_J": nominal_energy_J,
        "declared_input_to_liquid_sensible_capacity_ratio": nominal_energy_J / sensible_J,
        "actual_heat_rejection_inferred_J": None,
        "complete_thermal_energy_balance_closed": False,
        "interpretation": (
            "If the declared electrical estimate were actual deposited heat and all of it "
            "entered only the initial liquid at fixed pressure without rejection or "
            "chemical storage, it would exceed the liquid sensible energy to saturation. "
            "The real reactor also has gas, glass, vaporization, reaction heat, and "
            "environmental transfer; their unknown contributions prevent a temperature "
            "trajectory or heat-loss inference from this ratio."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_barth_water_thermal_capacity(), indent=2))
