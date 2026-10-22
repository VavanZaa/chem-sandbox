from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from .schema import CatalysisSimulationRequest


R_J_MOL_K = 8.31446261815324
K_B_J_K = 1.380649e-23
H_J_S = 6.62607015e-34


def _evidence_value(evidence, unit: str) -> float:
    value = evidence.as_unit(unit)
    if evidence.basis == "hypothetical":
        raise ValueError(f"Operation-ready catalysis cannot use hypothetical evidence from {evidence.source}")
    return value


def simulate_catalysis(request: CatalysisSimulationRequest) -> dict[str, Any]:
    temperature = request.temperature.as_unit("K")
    pressure = request.pressure.as_unit("Pa")
    duration = request.duration.as_unit("s")
    active_site_amount = request.active_site_amount.as_unit("mol")
    if min(temperature, pressure, duration, active_site_amount) <= 0:
        raise ValueError("Temperature, pressure, duration, and active-site amount must be positive")
    tolerance = request.detailed_balance_tolerance.as_unit("kJ/mol")
    if tolerance < 0:
        raise ValueError("Detailed-balance tolerance cannot be negative")

    lookup = {item.id: item for item in request.species}
    surface_ids = [item.id for item in request.species if item.phase == "surface"]
    gas_ids = [item.id for item in request.species if item.phase == "gas"]
    vacant = next(item.id for item in request.species if item.vacant_site)
    surface_index = {key: index for index, key in enumerate(surface_ids)}
    gas_index = {key: index for index, key in enumerate(gas_ids)}

    initial_surface = np.zeros(len(surface_ids), dtype=float)
    for key, value in request.initial_surface_amounts_per_site.items():
        initial_surface[surface_index[key]] = value
    occupied = sum(lookup[key].site_occupancy * initial_surface[index] for key, index in surface_index.items())
    if vacant not in request.initial_surface_amounts_per_site:
        initial_surface[surface_index[vacant]] = 1.0 - occupied
    site_total = sum(lookup[key].site_occupancy * initial_surface[index] for key, index in surface_index.items())
    if abs(site_total - 1.0) > 1e-9:
        raise ValueError(f"Initial surface site balance must equal one; received {site_total:.12g}")

    prefactor = K_B_J_K * temperature / H_J_S
    kinetics = []
    for reaction in request.reactions:
        forward_barrier = _evidence_value(reaction.forward_activation_gibbs, "kJ/mol")
        reverse_barrier = _evidence_value(reaction.reverse_activation_gibbs, "kJ/mol")
        reaction_gibbs = _evidence_value(reaction.reaction_gibbs, "kJ/mol")
        if forward_barrier < 0 or reverse_barrier < 0:
            raise ValueError(f"Reaction {reaction.id} has a negative activation free energy")
        cycle_residual = forward_barrier - reverse_barrier - reaction_gibbs
        if abs(cycle_residual) > tolerance:
            raise ValueError(
                f"Reaction {reaction.id} violates microscopic reversibility by {cycle_residual:.6g} kJ/mol"
            )
        k_forward = reaction.transmission_coefficient * prefactor * math.exp(
            -forward_barrier * 1000.0 / (R_J_MOL_K * temperature)
        )
        k_reverse = reaction.transmission_coefficient * prefactor * math.exp(
            -reverse_barrier * 1000.0 / (R_J_MOL_K * temperature)
        )
        kinetics.append((reaction, k_forward, k_reverse, cycle_residual))

    deactivation_rate = 0.0
    if request.deactivation.model == "first_order":
        deactivation_rate = _evidence_value(request.deactivation.rate_constant, "1/s")
        if deactivation_rate < 0:
            raise ValueError("Deactivation rate cannot be negative")
    effectiveness = 1.0 if request.transport_effectiveness is None else request.transport_effectiveness.value
    if request.transport_effectiveness and request.transport_effectiveness.basis == "hypothetical":
        raise ValueError("Operation-ready catalysis cannot use a hypothetical transport effectiveness")

    y0 = np.concatenate([initial_surface, np.zeros(len(gas_ids), dtype=float)])

    def activities(state: np.ndarray) -> dict[str, float]:
        values = {key: max(float(state[index]), 0.0) for key, index in surface_index.items()}
        values.update({key: request.gas_activities.get(key, 0.0) for key in gas_ids})
        return values

    def rates_at(time_s: float, state: np.ndarray) -> np.ndarray:
        act = activities(state)
        active_fraction = math.exp(-deactivation_rate * time_s)
        rates = []
        for reaction, k_forward, k_reverse, _residual in kinetics:
            forward_activity = 1.0
            reverse_activity = 1.0
            for key, coefficient in reaction.stoichiometry.items():
                if coefficient < 0:
                    forward_activity *= act[key] ** (-coefficient)
                else:
                    reverse_activity *= act[key] ** coefficient
            rates.append(effectiveness * active_fraction * (k_forward * forward_activity - k_reverse * reverse_activity))
        return np.asarray(rates, dtype=float)

    def rhs(time_s: float, state: np.ndarray) -> np.ndarray:
        rates = rates_at(time_s, state)
        derivative = np.zeros_like(state)
        for rate, (reaction, _kf, _kr, _residual) in zip(rates, kinetics, strict=True):
            for key, coefficient in reaction.stoichiometry.items():
                if key in surface_index:
                    derivative[surface_index[key]] += coefficient * rate
                else:
                    derivative[len(surface_ids) + gas_index[key]] += coefficient * rate
        return derivative

    times = np.linspace(0.0, duration, request.points)
    solution = solve_ivp(
        rhs,
        (0.0, duration),
        y0,
        method="BDF",
        t_eval=times,
        rtol=1e-9,
        atol=1e-12,
    )
    if not solution.success or not np.isfinite(solution.y).all():
        raise ValueError(f"Microkinetic integration failed: {solution.message}")
    if float(np.min(solution.y[: len(surface_ids)])) < -1e-7:
        raise ValueError("Microkinetic solution produced a materially negative surface population")

    samples = []
    max_site_residual = 0.0
    for column, time_s in enumerate(solution.t):
        state = solution.y[:, column]
        surface = {key: max(float(state[index]), 0.0) for key, index in surface_index.items()}
        site_balance = sum(lookup[key].site_occupancy * value for key, value in surface.items())
        max_site_residual = max(max_site_residual, abs(site_balance - 1.0))
        rates = rates_at(float(time_s), state)
        samples.append({
            "time_s": float(time_s),
            "active_fraction": math.exp(-deactivation_rate * float(time_s)),
            "surface_amounts_per_site": surface,
            "site_balance": site_balance,
            "gas_net_change_mol": {
                key: float(state[len(surface_ids) + gas_index[key]]) * active_site_amount for key in gas_ids
            },
            "net_rates_per_site_s": {
                reaction.id: float(rate) for rate, (reaction, _kf, _kr, _residual) in zip(rates, kinetics, strict=True)
            },
        })

    final_state = solution.y[:, -1]
    final_rates = rates_at(duration, final_state)
    target_tof = {}
    for target in request.target_gas_species:
        production = 0.0
        for rate, (reaction, _kf, _kr, _residual) in zip(final_rates, kinetics, strict=True):
            production += reaction.stoichiometry.get(target, 0.0) * float(rate)
        target_tof[target] = production
    return {
        "status": "completed",
        "accepted": max_site_residual <= 1e-7,
        "model": "mean-field reversible microkinetics",
        "temperature_K": temperature,
        "pressure_Pa": pressure,
        "active_site_amount_mol": active_site_amount,
        "transport_effectiveness": effectiveness,
        "samples": samples,
        "target_turnover_frequencies_s_inv": target_tof,
        "site_balance": {"closed": max_site_residual <= 1e-7, "max_absolute_residual": max_site_residual},
        "kinetics": [
            {
                "reaction_id": reaction.id,
                "forward_rate_constant_s_inv": k_forward,
                "reverse_rate_constant_s_inv": k_reverse,
                "cycle_residual_kj_mol": residual,
            }
            for reaction, k_forward, k_reverse, residual in kinetics
        ],
        "boundary": "Post-discovery catalyst model; all elementary steps and energetic evidence were supplied explicitly.",
        "limitations": [
            "Gas activities are fixed reservoirs; gas depletion and reactor transport require coupling to a process model.",
            "Mean-field activities omit lateral interactions unless they are already represented in supplied free energies.",
            "A catalyst result is accepted only with explicit bound intermediates and a closed site balance.",
        ],
    }


def export_catmap_setup(request: CatalysisSimulationRequest) -> dict[str, Any]:
    lookup = {item.id: item for item in request.species}
    vacant = next(item.id for item in request.species if item.vacant_site)

    def catmap_name(key: str) -> str:
        item = lookup[key]
        if key == vacant:
            return "*"
        return f"{key}_{'g' if item.phase == 'gas' else 's'}"

    def side(stoichiometry: dict[str, float], sign: int) -> str:
        tokens = []
        for key, coefficient in stoichiometry.items():
            if coefficient * sign <= 0:
                continue
            magnitude = abs(coefficient)
            if abs(magnitude - round(magnitude)) > 1e-10:
                raise ValueError("CatMAP export currently requires integer stoichiometric coefficients")
            count = int(round(magnitude))
            tokens.extend([catmap_name(key)] * count)
        return " + ".join(tokens)

    expressions = [
        f"{side(reaction.stoichiometry, -1)} <-> {side(reaction.stoichiometry, 1)}"
        for reaction in request.reactions
    ]
    definitions = {"s": {"site_names": ["submitted_surface"], "total": 1.0}}
    for item in request.species:
        if item.phase == "gas":
            definitions[catmap_name(item.id)] = {"pressure": request.gas_activities.get(item.id, 0.0)}
        elif not item.vacant_site and item.site_occupancy != 1:
            definitions[catmap_name(item.id)] = {"n_sites": item.site_occupancy}
    lines = ["rxn_expressions = ["]
    lines.extend(f"    {expression!r}," for expression in expressions)
    lines.append("]")
    lines.append(f"species_definitions = {definitions!r}")
    lines.append(f"temperature = {request.temperature.as_unit('K')!r}")
    return {
        "format": "CatMAP setup fragment",
        "setup_text": "\n".join(lines) + "\n",
        "reaction_expressions": expressions,
        "note": "Supply a CatMAP energetic input table separately; this export does not invent formation energies or scaling relations.",
    }
