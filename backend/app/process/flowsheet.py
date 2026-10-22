from __future__ import annotations

from typing import Any

from .idaes_bridge import solve_idaes_algebraic
from .schema import ProcessFlowsheetRequest


def _fraction(evidence, label: str) -> tuple[float, float | None]:
    if evidence.basis == "hypothetical":
        raise ValueError(f"Operation-ready flowsheets cannot use hypothetical {label}")
    value = evidence.as_unit("dimensionless")
    if not 0 <= value <= 1:
        raise ValueError(f"{label} must lie between zero and one")
    uncertainty = None
    if evidence.uncertainty is not None:
        uncertainty = evidence.uncertainty.as_unit("dimensionless")
        if uncertainty < 0:
            raise ValueError(f"{label} uncertainty cannot be negative")
    return value, uncertainty


def _solve_once(
    request: ProcessFlowsheetRequest,
    conversion: float,
    recoveries: dict[str, float],
    recycle_fraction: float,
) -> dict[str, Any]:
    components = [item.id for item in request.components]
    fresh = {key: request.fresh_feed.get(key, None) for key in components}
    fresh_si = {key: 0.0 if quantity is None else quantity.as_unit("mol/s") for key, quantity in fresh.items()}
    if any(value < 0 for value in fresh_si.values()) or sum(fresh_si.values()) <= 0:
        raise ValueError("Fresh feed requires nonnegative component flows and at least one positive flow")
    stoichiometry = request.reaction.stoichiometry
    limiting = request.reaction.limiting_species
    limiting_coefficient = -stoichiometry[limiting]
    recycle = {key: 0.0 for key in components}
    converged = False
    iteration = 0
    for iteration in range(1, request.max_iterations + 1):
        mixed = {key: fresh_si[key] + recycle[key] for key in components}
        extent = conversion * mixed[limiting] / limiting_coefficient
        reactor_outlet = {key: mixed[key] + stoichiometry.get(key, 0.0) * extent for key in components}
        if min(reactor_outlet.values()) < -request.convergence_tolerance:
            raise ValueError(
                "Declared limiting species/conversion consumes another reactant below zero; "
                "reconcile feed stoichiometry or choose the actual limiting species"
            )
        reactor_outlet = {key: max(value, 0.0) for key, value in reactor_outlet.items()}
        product = {key: recoveries[key] * reactor_outlet[key] for key in components}
        next_recycle = {
            key: recycle_fraction * (reactor_outlet[key] - product[key]) for key in components
        }
        scale = max(sum(abs(value) for value in next_recycle.values()), 1.0)
        residual = max(abs(next_recycle[key] - recycle[key]) for key in components) / scale
        recycle = next_recycle
        if residual <= request.convergence_tolerance:
            converged = True
            break
    if not converged:
        raise ValueError("Recycle flowsheet did not converge within max_iterations")
    mixed = {key: fresh_si[key] + recycle[key] for key in components}
    extent = conversion * mixed[limiting] / limiting_coefficient
    reactor_outlet = {key: max(mixed[key] + stoichiometry.get(key, 0.0) * extent, 0.0) for key in components}
    product = {key: recoveries[key] * reactor_outlet[key] for key in components}
    recycle = {key: recycle_fraction * (reactor_outlet[key] - product[key]) for key in components}
    purge = {key: (1.0 - recycle_fraction) * (reactor_outlet[key] - product[key]) for key in components}
    internal_flow = sum(mixed.values())
    max_internal = None
    if request.max_internal_molar_flow is not None:
        max_internal = request.max_internal_molar_flow.as_unit("mol/s")
        if max_internal <= 0:
            raise ValueError("Maximum internal flow must be positive")
    feasible = max_internal is None or internal_flow <= max_internal * (1 + 1e-10)

    lookup = {item.id: item for item in request.components}
    element_residual: dict[str, float] = {}
    for key in components:
        net = product[key] + purge[key] - fresh_si[key]
        for element, count in lookup[key].elements.items():
            element_residual[element] = element_residual.get(element, 0.0) + net * count
    scale = max(sum(fresh_si.values()), 1.0)
    max_element_residual = max((abs(value) / scale for value in element_residual.values()), default=0.0)
    utility = None
    if request.utility_per_extent is not None:
        if request.utility_per_extent.basis == "hypothetical":
            raise ValueError("Operation-ready flowsheets cannot use hypothetical utility evidence")
        utility = extent * request.utility_per_extent.as_unit("J/mol")
    fresh_limiting = fresh_si[limiting]
    overall_conversion = None if fresh_limiting <= 0 else 1.0 - (product[limiting] + purge[limiting]) / fresh_limiting
    return {
        "converged": converged,
        "feasible": feasible,
        "iterations": iteration,
        "recycle_fraction": recycle_fraction,
        "single_pass_conversion": conversion,
        "fresh_feed_mol_s": fresh_si,
        "mixed_feed_mol_s": mixed,
        "reactor_outlet_mol_s": reactor_outlet,
        "product_mol_s": product,
        "recycle_mol_s": recycle,
        "purge_mol_s": purge,
        "reaction_extent_mol_s": extent,
        "internal_molar_flow_mol_s": internal_flow,
        "max_internal_molar_flow_mol_s": max_internal,
        "target_product_mol_s": product[request.target_species],
        "overall_limiting_species_conversion": overall_conversion,
        "utility_W": utility,
        "material_balance": {
            "closed": max_element_residual <= 1e-8,
            "max_relative_element_residual": max_element_residual,
            "element_residual_mol_s": element_residual,
        },
    }


def simulate_process_flowsheet(request: ProcessFlowsheetRequest) -> dict[str, Any]:
    conversion, conversion_uncertainty = _fraction(request.single_pass_conversion, "single-pass conversion")
    recoveries: dict[str, float] = {}
    recovery_uncertainties: dict[str, float | None] = {}
    for key, evidence in request.product_recoveries.items():
        recoveries[key], recovery_uncertainties[key] = _fraction(evidence, f"{key} product recovery")
    nominal_recycle, recycle_uncertainty = _fraction(request.recycle.fraction, "recycle fraction")

    candidates: list[dict[str, Any]] = []
    if request.recycle.optimize:
        for index in range(request.recycle.grid_points):
            fraction = request.recycle.lower_bound + (
                request.recycle.upper_bound - request.recycle.lower_bound
            ) * index / (request.recycle.grid_points - 1)
            candidate = _solve_once(request, conversion, recoveries, fraction)
            if candidate["feasible"]:
                candidates.append(candidate)
        if not candidates:
            raise ValueError("No recycle setting satisfies the declared internal-flow constraint")
        nominal = max(candidates, key=lambda item: item["target_product_mol_s"])
        optimization = {
            "enabled": True,
            "method": "bounded deterministic grid",
            "feasible_points": len(candidates),
            "objective": "maximize target product molar flow",
        }
    else:
        nominal = _solve_once(request, conversion, recoveries, nominal_recycle)
        optimization = {"enabled": False, "method": None, "feasible_points": int(nominal["feasible"])}
    if not nominal["feasible"]:
        raise ValueError("Nominal flowsheet exceeds max_internal_molar_flow")

    scenarios = []
    has_uncertainty = conversion_uncertainty is not None or recycle_uncertainty is not None or any(
        value is not None for value in recovery_uncertainties.values()
    )
    if has_uncertainty:
        for label, sign in (("low_parameters", -1.0), ("high_parameters", 1.0)):
            scenario_conversion = min(1.0, max(0.0, conversion + sign * (conversion_uncertainty or 0.0)))
            scenario_recycle = min(0.999999, max(0.0, nominal["recycle_fraction"] + sign * (recycle_uncertainty or 0.0)))
            scenario_recoveries = {
                key: min(1.0, max(0.0, value + sign * (recovery_uncertainties[key] or 0.0)))
                for key, value in recoveries.items()
            }
            scenario = _solve_once(request, scenario_conversion, scenario_recoveries, scenario_recycle)
            scenarios.append({
                "scenario": label,
                "feasible": scenario["feasible"],
                "target_product_mol_s": scenario["target_product_mol_s"],
                "overall_conversion": scenario["overall_limiting_species_conversion"],
                "internal_molar_flow_mol_s": scenario["internal_molar_flow_mol_s"],
            })

    idaes = solve_idaes_algebraic(
        list(nominal["fresh_feed_mol_s"]),
        nominal["fresh_feed_mol_s"],
        request.reaction.stoichiometry,
        request.reaction.limiting_species,
        nominal["single_pass_conversion"],
        recoveries,
        nominal["recycle_fraction"],
    )
    agreement = None
    if idaes.get("status") == "completed":
        fields = [
            abs(idaes["product_mol_s"][key] - nominal["product_mol_s"][key])
            for key in nominal["product_mol_s"]
        ]
        scale = max(sum(nominal["product_mol_s"].values()), 1.0)
        agreement = max(fields, default=0.0) / scale
    crosscheck_required = bool(
        idaes.get("capabilities", {}).get("idaes_available")
        and idaes.get("capabilities", {}).get("ipopt_available")
    )
    crosscheck_ok = (
        not crosscheck_required
        or (idaes.get("status") == "completed" and agreement is not None and agreement <= 1e-7)
    )
    accepted = bool(
        nominal["converged"]
        and nominal["feasible"]
        and nominal["material_balance"]["closed"]
        and crosscheck_ok
    )
    return {
        "status": "completed" if accepted else "crosscheck_failure",
        "accepted": accepted,
        "flowsheet": nominal,
        "optimization": optimization,
        "uncertainty_scenarios": scenarios,
        "idaes_crosscheck": idaes,
        "cross_engine_max_relative_product_difference": agreement,
        "boundary": "Post-discovery process model using only the submitted reaction, conversion, separation, and utility evidence.",
        "limitations": [
            "This foundation models one reactor, component separator, purge, and recycle loop at steady state.",
            "Thermodynamic property packages, phase equilibrium, equipment sizing, pressure drop, and safety relief require later unit-specific models.",
            "Optimization never changes chemistry; it only searches the declared recycle bounds under the declared internal-flow constraint.",
        ],
    }
