from __future__ import annotations

import math
from typing import Any

import numpy as np

from .schema import PressureDependence, RateArtifactRequest


R_J_MOL_K = 8.31446261815324
K_B_J_K = 1.380649e-23
H_J_S = 6.62607015e-34


def _rate_unit(order: int) -> str:
    if order == 1:
        return "1/s"
    if order == 2:
        return "m^3/(mol*s)"
    return "m^6/(mol^2*s)"


def _pressure_multipliers(model: PressureDependence, pressure_Pa: float | None) -> tuple[float, float, list[str]]:
    if model.model == "none":
        return 1.0, 1.0, []
    problems: list[str] = []
    if model.basis == "hypothetical":
        problems.append("pressure-dependence evidence is hypothetical")
    if pressure_Pa is None:
        problems.append("tabulated pressure dependence requires a simulation pressure")
        return 1.0, 1.0, problems
    rows = sorted(
        (
            point.pressure.as_unit("Pa"),
            point.forward_multiplier,
            point.reverse_multiplier,
        )
        for point in model.points
    )
    if pressure_Pa < rows[0][0] or pressure_Pa > rows[-1][0]:
        problems.append("requested pressure lies outside the tabulated pressure range; extrapolation is forbidden")
        return 1.0, 1.0, problems
    pressures = np.log([row[0] for row in rows])
    forward = math.exp(float(np.interp(math.log(pressure_Pa), pressures, np.log([row[1] for row in rows]))))
    reverse = math.exp(float(np.interp(math.log(pressure_Pa), pressures, np.log([row[2] for row in rows]))))
    return forward, reverse, problems


def _uncertainty_kj_mol(quantity) -> float | None:
    if quantity.uncertainty is None:
        return None
    value = quantity.uncertainty.as_unit("kJ/mol")
    if value < 0:
        raise ValueError("Energy uncertainty cannot be negative")
    return value


def build_rate_artifact(request: RateArtifactRequest) -> dict[str, Any]:
    temperature_K = request.temperature.as_unit("K")
    if temperature_K <= 0:
        raise ValueError("Temperature must be greater than 0 K")
    pressure_Pa = None if request.pressure is None else request.pressure.as_unit("Pa")
    if pressure_Pa is not None and pressure_Pa <= 0:
        raise ValueError("Pressure must be greater than zero")
    standard_concentration = request.standard_concentration.as_unit("mol/m^3")
    if standard_concentration <= 0:
        raise ValueError("Standard concentration must be greater than zero")
    tolerance = request.detailed_balance_tolerance.as_unit("kJ/mol")
    if tolerance < 0:
        raise ValueError("Detailed-balance tolerance cannot be negative")

    forward_barrier = request.forward_activation_gibbs.as_unit("kJ/mol")
    reverse_barrier = request.reverse_activation_gibbs.as_unit("kJ/mol")
    reaction_gibbs = request.reaction_gibbs.as_unit("kJ/mol")
    problems: list[str] = []
    if forward_barrier > 500 or reverse_barrier > 500 or abs(reaction_gibbs) > 500:
        raise ValueError("Rate-artifact Gibbs energies are limited to ±500 kJ/mol for numerical safety")
    for label, evidence in (
        ("forward activation Gibbs energy", request.forward_activation_gibbs),
        ("reverse activation Gibbs energy", request.reverse_activation_gibbs),
        ("reaction Gibbs energy", request.reaction_gibbs),
    ):
        if evidence.basis == "hypothetical":
            problems.append(f"{label} is hypothetical")
    if forward_barrier < 0 or reverse_barrier < 0:
        problems.append("activation Gibbs energies cannot be negative")
    cycle_residual = forward_barrier - reverse_barrier - reaction_gibbs
    if abs(cycle_residual) > tolerance:
        problems.append(
            f"microscopic reversibility failed: ΔG‡f - ΔG‡r - ΔGr = {cycle_residual:.6g} kJ/mol"
        )

    forward_pressure, reverse_pressure, pressure_problems = _pressure_multipliers(
        request.pressure_dependence, pressure_Pa
    )
    problems.extend(pressure_problems)
    prefactor = request.transmission_coefficient * K_B_J_K * temperature_K / H_J_S
    forward_rate = (
        prefactor
        * math.exp(-forward_barrier * 1000.0 / (R_J_MOL_K * temperature_K))
        / standard_concentration ** (request.forward_order - 1)
        * forward_pressure
    )
    reverse_rate = (
        prefactor
        * math.exp(-reverse_barrier * 1000.0 / (R_J_MOL_K * temperature_K))
        / standard_concentration ** (request.reverse_order - 1)
        * reverse_pressure
    )
    equilibrium_constant = math.exp(-reaction_gibbs * 1000.0 / (R_J_MOL_K * temperature_K))
    pressure_cycle_residual = (
        -cycle_residual
        + R_J_MOL_K * temperature_K / 1000.0 * math.log(forward_pressure / reverse_pressure)
    )
    if abs(pressure_cycle_residual) > tolerance:
        problems.append(
            "pressure multipliers violate detailed balance at the requested pressure: "
            f"residual={pressure_cycle_residual:.6g} kJ/mol"
        )

    forward_uncertainty = _uncertainty_kj_mol(request.forward_activation_gibbs)
    reverse_uncertainty = _uncertainty_kj_mol(request.reverse_activation_gibbs)
    reaction_uncertainty = _uncertainty_kj_mol(request.reaction_gibbs)

    def rate_range(rate: float, uncertainty: float | None) -> list[float] | None:
        if uncertainty is None:
            return None
        factor = math.exp(uncertainty * 1000.0 / (R_J_MOL_K * temperature_K))
        return [rate / factor, rate * factor]

    uncertainty_values = [value for value in (forward_uncertainty, reverse_uncertainty, reaction_uncertainty) if value is not None]
    largest_factor = max(
        (math.exp(value * 1000.0 / (R_J_MOL_K * temperature_K)) for value in uncertainty_values),
        default=None,
    )
    active_learning_score = None if largest_factor is None else request.impact_weight * math.log(largest_factor)
    certified = not problems
    return {
        "schema_version": "1.0",
        "reaction_id": request.reaction_id,
        "certified": certified,
        "operation_ready": certified and len(uncertainty_values) == 3,
        "temperature_K": temperature_K,
        "pressure_Pa": pressure_Pa,
        "forward": {
            "order": request.forward_order,
            "rate_constant": forward_rate if certified else None,
            "unit": _rate_unit(request.forward_order),
            "uncertainty_range": rate_range(forward_rate, forward_uncertainty) if certified else None,
        },
        "reverse": {
            "order": request.reverse_order,
            "rate_constant": reverse_rate if certified else None,
            "unit": _rate_unit(request.reverse_order),
            "uncertainty_range": rate_range(reverse_rate, reverse_uncertainty) if certified else None,
        },
        "thermodynamics": {
            "reaction_gibbs_kj_mol": reaction_gibbs,
            "equilibrium_constant_dimensionless": equilibrium_constant,
            "cycle_residual_kj_mol": cycle_residual,
            "pressure_adjusted_cycle_residual_kj_mol": pressure_cycle_residual,
        },
        "pressure_dependence": request.pressure_dependence.model,
        "uncertainty": {
            "status": "complete" if len(uncertainty_values) == 3 else "incomplete",
            "largest_multiplicative_rate_factor": largest_factor,
        },
        "active_learning": {
            "impact_weight": request.impact_weight,
            "priority_score": active_learning_score,
            "reason": "Barrier uncertainty is prioritized by declared pathway impact; this is a queue score, not confidence.",
        },
        "problems": problems,
        "evidence": request.model_dump(mode="json"),
        "claim": "Detailed-balance-consistent Eyring rate artifact." if certified else "No certified rate constants are claimed.",
    }
