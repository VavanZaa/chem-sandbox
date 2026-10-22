from __future__ import annotations

from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from .schema import GasLiquidTransferRequest


R_J_MOL_K = 8.31446261815324


def _simulate(request: GasLiquidTransferRequest, henry_hcp: float, kla_s_inv: float) -> dict[str, Any]:
    temperature_K = request.temperature.as_unit("K")
    total_pressure_Pa = request.total_pressure.as_unit("Pa")
    gas_volume_m3 = request.gas_volume.as_unit("m^3")
    liquid_volume_m3 = request.liquid_volume.as_unit("m^3")
    duration_s = request.duration.as_unit("s")
    initial_gas_mol = request.initial_gas_amount.as_unit("mol")
    initial_liquid_mol = request.initial_dissolved_amount.as_unit("mol")
    values = {
        "temperature": temperature_K,
        "total pressure": total_pressure_Pa,
        "gas volume": gas_volume_m3,
        "liquid volume": liquid_volume_m3,
        "duration": duration_s,
    }
    invalid = [name for name, value in values.items() if value <= 0]
    if invalid:
        raise ValueError(f"Gas-liquid quantities must be greater than zero: {', '.join(invalid)}")
    if initial_gas_mol < 0 or initial_liquid_mol < 0:
        raise ValueError("Initial gas and dissolved amounts cannot be negative")
    if henry_hcp <= 0 or kla_s_inv < 0:
        raise ValueError("Henry Hcp must be positive and kLa cannot be negative")
    initial_partial_pressure = initial_gas_mol * R_J_MOL_K * temperature_K / gas_volume_m3
    if request.mode == "closed_headspace" and initial_partial_pressure > total_pressure_Pa * (1 + 1e-8):
        raise ValueError("Initial solute partial pressure exceeds total pressure")
    ambient_pressure = 0.0
    if request.ambient_partial_pressure is not None:
        ambient_pressure = request.ambient_partial_pressure.as_unit("Pa")
        if ambient_pressure < 0 or ambient_pressure > total_pressure_Pa:
            raise ValueError("Ambient partial pressure must lie between zero and total pressure")
    total_closed = initial_gas_mol + initial_liquid_mol

    def rhs(_time, state):
        dissolved = float(state[0])
        if request.mode == "closed_headspace":
            gas_mol = max(total_closed - dissolved, 0.0)
            partial_pressure = gas_mol * R_J_MOL_K * temperature_K / gas_volume_m3
        else:
            partial_pressure = ambient_pressure
        equilibrium_concentration = henry_hcp * partial_pressure
        actual_concentration = dissolved / liquid_volume_m3
        return [kla_s_inv * liquid_volume_m3 * (equilibrium_concentration - actual_concentration)]

    times = np.linspace(0.0, duration_s, request.points)
    solution = solve_ivp(
        rhs, (0.0, duration_s), [initial_liquid_mol], method="BDF", t_eval=times,
        rtol=1e-9, atol=max((initial_gas_mol + initial_liquid_mol) * 1e-12, 1e-15),
    )
    if not solution.success or not np.isfinite(solution.y).all():
        raise ValueError(f"Gas-liquid integration failed: {solution.message}")
    samples = []
    for time_s, dissolved in zip(solution.t, solution.y[0], strict=True):
        dissolved = max(float(dissolved), 0.0)
        if request.mode == "closed_headspace":
            gas_mol = max(total_closed - dissolved, 0.0)
            partial_pressure = gas_mol * R_J_MOL_K * temperature_K / gas_volume_m3
            external_exchange = 0.0
        else:
            gas_mol = None
            partial_pressure = ambient_pressure
            external_exchange = dissolved - initial_liquid_mol
        equilibrium_concentration = henry_hcp * partial_pressure
        concentration = dissolved / liquid_volume_m3
        samples.append(
            {
                "time_s": float(time_s),
                "dissolved_mol": dissolved,
                "headspace_mol": gas_mol,
                "partial_pressure_Pa": partial_pressure,
                "dissolved_concentration_mol_m3": concentration,
                "equilibrium_concentration_mol_m3": equilibrium_concentration,
                "external_exchange_mol": external_exchange,
            }
        )
    final = samples[-1]
    if request.mode == "closed_headspace":
        residual = final["dissolved_mol"] + final["headspace_mol"] - total_closed
        mass_balance = {
            "closed": abs(residual) <= max(total_closed, 1.0) * 1e-10,
            "residual_mol": residual,
            "basis": "closed solute inventory",
        }
    else:
        residual = final["dissolved_mol"] - initial_liquid_mol - final["external_exchange_mol"]
        mass_balance = {
            "closed": abs(residual) <= max(abs(final["external_exchange_mol"]), initial_liquid_mol, 1.0) * 1e-10,
            "residual_mol": residual,
            "basis": "liquid inventory plus exchange with fixed-pressure gas reservoir",
        }
    return {
        "samples": samples,
        "mass_balance": mass_balance,
        "final_equilibrium_residual_mol_m3": (
            final["dissolved_concentration_mol_m3"] - final["equilibrium_concentration_mol_m3"]
        ),
    }


def simulate_gas_liquid_transfer(request: GasLiquidTransferRequest) -> dict[str, Any]:
    henry = request.henry_hcp.as_unit("mol/(m^3*Pa)")
    kla = request.kla.as_unit("1/s")
    nominal = _simulate(request, henry, kla)
    problems = []
    if request.henry_hcp.basis == "hypothetical":
        problems.append("Henry-law evidence is hypothetical")
    if request.kla.basis == "hypothetical":
        problems.append("kLa evidence is hypothetical")
    scenarios = []
    if request.henry_hcp.uncertainty is not None:
        uncertainty = request.henry_hcp.uncertainty.as_unit("mol/(m^3*Pa)")
        if uncertainty < 0 or henry - uncertainty <= 0:
            raise ValueError("Henry-law uncertainty must be nonnegative and smaller than the nominal value")
        for label, value in (("henry_low", henry - uncertainty), ("henry_high", henry + uncertainty)):
            scenario = _simulate(request, value, kla)
            scenarios.append(
                {
                    "scenario": label,
                    "henry_hcp_mol_m3_Pa": value,
                    "final": scenario["samples"][-1],
                }
            )
    accepted = not problems and nominal["mass_balance"]["closed"]
    return {
        "status": "completed" if accepted else "conditional_only",
        "accepted": accepted,
        "solute_id": request.solute_id,
        "mode": request.mode,
        "temperature_K": request.temperature.as_unit("K"),
        "henry_hcp_mol_m3_Pa": henry,
        "kla_s_inv": kla,
        **nominal,
        "uncertainty_scenarios": scenarios,
        "problems": problems,
        "evidence": request.model_dump(mode="json"),
        "limitations": [
            "This model transfers one declared neutral solute using Henry's law and a supplied volumetric kLa.",
            "It does not calculate acid/base speciation, ionic activities, salting-out, chemical enhancement, or precipitation.",
            "Use the isolated Reaktoro equilibrium service for electrolyte chemistry; missing service results are never replaced by ideal speciation.",
        ],
    }
