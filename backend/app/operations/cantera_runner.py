from __future__ import annotations

import hashlib
import threading
import warnings
from pathlib import Path
from typing import Any

import cantera as ct
import numpy as np
from scipy.optimize import OptimizeWarning, linprog

from .cantera_schema import CanteraSimulationRequest
from .validation import validate_operation


_ceiling_lock = threading.Lock()


def _load_solution(request: CanteraSimulationRequest) -> tuple[ct.Solution, dict[str, Any]]:
    source = request.cantera
    if source.yaml_text is not None:
        gas = ct.Solution(yaml=source.yaml_text, name=source.phase_name) if source.phase_name else ct.Solution(yaml=source.yaml_text)
        digest = hashlib.sha256(source.yaml_text.encode("utf-8")).hexdigest()
        identity = "submitted-yaml"
    else:
        gas = ct.Solution(source.builtin, source.phase_name) if source.phase_name else ct.Solution(source.builtin)
        model_path = next(
            (Path(directory) / source.builtin for directory in ct.get_data_directories()
             if (Path(directory) / source.builtin).is_file()),
            None,
        )
        digest = hashlib.sha256(model_path.read_bytes()).hexdigest() if model_path else None
        identity = source.builtin
    return gas, {
        "cantera_version": ct.__version__,
        "mechanism": identity,
        "mechanism_sha256": digest,
        "source": source.source,
        "phase_name": gas.name,
    }


def _phase_and_composition(request: CanteraSimulationRequest) -> tuple[dict[str, float], dict[str, float]]:
    operation = request.definition.operation
    gas_phases = [phase for phase in operation.phases if phase.kind == "gas"]
    if len(operation.phases) != 1 or len(gas_phases) != 1:
        raise ValueError("The native Cantera runner currently requires exactly one ideal-gas phase")
    if operation.gas_transfer.model != "none":
        raise ValueError("Gas transfer belongs to the multiphase runner, not a homogeneous Cantera reactor")
    if operation.catalyst.mode != "none":
        raise ValueError("Catalytic Cantera runs require the phase-6 surface mechanism adapter")
    inventory: dict[str, float] = {}
    for item in operation.initial_inventory:
        inventory[item.species_id] = inventory.get(item.species_id, 0.0) + item.amount.as_unit("mol")
    feed: dict[str, float] = {}
    for stream in operation.feeds:
        for item in stream.species_flows:
            feed[item.species_id] = feed.get(item.species_id, 0.0) + item.flow.as_unit("mol/s")
    return inventory, feed


def _species_moles(reactor) -> dict[str, float]:
    phase = reactor.phase
    moles = reactor.mass * phase.Y / phase.molecular_weights * 1000.0
    return {name: float(value) for name, value in zip(phase.species_names, moles, strict=True)}


def _species_molar_flows(phase: ct.Solution, mass_flow_kg_s: float) -> dict[str, float]:
    flows = mass_flow_kg_s * phase.Y / phase.molecular_weights * 1000.0
    return {name: float(value) for name, value in zip(phase.species_names, flows, strict=True)}


def _snapshot(reactor, coordinate: float, coordinate_name: str, flow_kg_s: float | None = None) -> dict[str, Any]:
    phase = reactor.phase
    payload = {
        coordinate_name: float(coordinate),
        "temperature_K": float(reactor.T),
        "pressure_Pa": float(phase.P),
        "volume_m3": float(reactor.volume),
        "mole_fractions": {
            name: float(value) for name, value in zip(phase.species_names, phase.X, strict=True) if value > 1e-16
        },
    }
    if flow_kg_s is None:
        payload["species_moles"] = _species_moles(reactor)
    else:
        payload["species_molar_flows_mol_s"] = _species_molar_flows(phase, flow_kg_s)
    return payload


def _element_inventory(gas: ct.Solution, amounts: dict[str, float]) -> dict[str, float]:
    result = {element: 0.0 for element in gas.element_names}
    for species, amount in amounts.items():
        if species not in gas.species_names:
            continue
        for element in gas.element_names:
            result[element] += amount * gas.n_atoms(species, element)
    return {key: value for key, value in result.items() if abs(value) > 1e-14}


def _balance(initial: dict[str, float], final: dict[str, float], gas: ct.Solution, basis: str) -> dict[str, Any]:
    before = _element_inventory(gas, initial)
    after = _element_inventory(gas, final)
    residual = {element: after.get(element, 0.0) - before.get(element, 0.0) for element in set(before) | set(after)}
    scale = max([abs(value) for value in before.values()] + [1.0])
    relative = max((abs(value) / scale for value in residual.values()), default=0.0)
    return {
        "basis": basis,
        "closed": relative <= 1e-6,
        "max_relative_element_residual": relative,
        "element_residual": {key: value for key, value in residual.items() if abs(value) > scale * 1e-10},
    }


def _theoretical_ceiling(definition, initial: dict[str, float], target: str) -> float | None:
    species_ids = [item.id for item in definition.mechanism.species]
    if target not in species_ids:
        return None
    index = {name: i for i, name in enumerate(species_ids)}
    matrix = np.zeros((len(species_ids), len(definition.mechanism.reactions)))
    for column, reaction in enumerate(definition.mechanism.reactions):
        for species, coefficient in reaction.stoichiometry.items():
            matrix[index[species], column] = coefficient
    start = np.asarray([initial.get(name, 0.0) for name in species_ids], dtype=float)
    with _ceiling_lock, warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Unrecognized options detected.*threads.*", category=OptimizeWarning)
        result = linprog(
            -matrix[index[target]], A_ub=-matrix, b_ub=start,
            bounds=(0, None), method="highs", options={"threads": 1},
        )
    return max(0.0, -float(result.fun)) if result.success else None


def _metrics(request: CanteraSimulationRequest, initial: dict[str, float], final: dict[str, float]) -> dict[str, Any]:
    metrics: dict[str, Any] = {}
    reference = request.reference_species
    if reference:
        if reference not in initial or initial[reference] <= 0:
            raise ValueError("Reference species requires a positive inlet/initial amount")
        metrics["reference_species"] = reference
        metrics["conversion"] = (initial[reference] - final.get(reference, 0.0)) / initial[reference]
    targets = {}
    for target in request.target_species:
        formed = final.get(target, 0.0) - initial.get(target, 0.0)
        ceiling = _theoretical_ceiling(request.definition, initial, target)
        targets[target] = {
            "net_formation": formed,
            "stoichiometric_ceiling": ceiling,
            "yield_fraction": None if ceiling is None or ceiling <= 1e-15 else formed / ceiling,
        }
    metrics["targets"] = targets
    return metrics


def _energy_balance(
    initial_energy: float,
    final_energy: float,
    basis: str,
    isothermal: bool,
    solver_relative_tolerance: float,
) -> dict[str, Any]:
    change = final_energy - initial_energy
    is_flow = "flow" in basis
    suffix = "W" if is_flow else "J"
    if isothermal:
        return {
            "closed": True,
            "basis": basis,
            f"state_function_change_{suffix}": change,
            f"required_external_heat_{suffix}": change,
            "note": "The reported external heat is the state-function duty needed to maintain the isothermal constraint.",
        }
    scale = max(abs(initial_energy), 1.0)
    # ReactorNet controls local integration error, while the balance compares two
    # independently reconstructed extensive states. Reactive PFR trajectories can
    # accumulate a few hundred local tolerances without violating conservation.
    # Keep the acceptance band bounded and always report it for auditability.
    relative_tolerance = max(1e-6, min(1e-5, 500.0 * solver_relative_tolerance))
    return {
        "closed": abs(change) / scale <= relative_tolerance,
        "basis": basis,
        f"residual_{suffix}": change,
        "relative_residual": abs(change) / scale,
        "relative_tolerance": relative_tolerance,
    }


def simulate_cantera(request: CanteraSimulationRequest) -> dict[str, Any]:
    validation = validate_operation(request.definition)
    if not validation["ready_for_solver"]:
        raise ValueError("Operation mechanism failed conservation checks")
    operation = request.definition.operation
    if operation.reactor_type == "semibatch":
        raise ValueError("Semibatch Cantera integration is not yet implemented; no batch approximation is substituted")
    if operation.heat_transfer.mode == "overall_UA":
        raise ValueError("overall_UA trajectories require integrated wall heat accounting and are not yet enabled")
    gas, provenance = _load_solution(request)
    if gas.thermo_model != "ideal-gas":
        raise ValueError("This runner currently supports Cantera ideal-gas phases only")
    inventory, feed = _phase_and_composition(request)
    declared_species = {item.id for item in request.definition.mechanism.species}
    unknown = declared_species - set(gas.species_names)
    if unknown:
        raise ValueError(f"Declared species are absent from the Cantera phase: {sorted(unknown)}")
    selected = set(request.target_species) | ({request.reference_species} if request.reference_species else set())
    if selected - set(gas.species_names):
        raise ValueError(f"Metric species are absent from the Cantera phase: {sorted(selected - set(gas.species_names))}")

    energy_mode = "on" if operation.heat_transfer.mode == "adiabatic" else "off"
    temperature_K = validation["normalized_si"]["temperature_K"]
    pressure_Pa = validation["normalized_si"]["pressure_Pa"]
    volume_m3 = validation["normalized_si"]["volume_m3"]
    network = None
    samples: list[dict[str, Any]] = []

    if operation.reactor_type == "batch":
        if not inventory or sum(inventory.values()) <= 0:
            raise ValueError("Batch Cantera simulation requires positive gas inventory")
        gas.TPX = temperature_K, pressure_Pa, inventory
        reactor_class = ct.IdealGasMoleReactor if request.batch_constraint == "constant_volume" else ct.IdealGasConstPressureMoleReactor
        reactor = reactor_class(gas, energy=energy_mode, clone=True, volume=volume_m3)
        actual_initial = _species_moles(reactor)
        declared_total, actual_total = sum(inventory.values()), sum(actual_initial.values())
        mismatch = abs(actual_total - declared_total) / declared_total
        if mismatch > request.inventory_state_tolerance:
            raise ValueError(
                f"Initial amount/temperature/pressure/volume are inconsistent by {mismatch:.1%}; reconcile the gas state"
            )
        network = ct.ReactorNet([reactor])
        network.rtol, network.atol = request.relative_tolerance, request.absolute_tolerance
        duration = operation.duration.as_unit("s")
        coordinate = np.linspace(0.0, duration, request.points)
        initial_energy = reactor.mass * (
            reactor.phase.int_energy_mass if request.batch_constraint == "constant_volume" else reactor.phase.enthalpy_mass
        )
        for time_s in coordinate:
            network.advance(float(time_s))
            samples.append(_snapshot(reactor, time_s, "time_s"))
        final = samples[-1]["species_moles"]
        material_balance = _balance(actual_initial, final, gas, "closed-reactor moles")
        final_energy = reactor.mass * (
            reactor.phase.int_energy_mass if request.batch_constraint == "constant_volume" else reactor.phase.enthalpy_mass
        )
        energy_balance = _energy_balance(
            initial_energy,
            final_energy,
            "internal energy" if request.batch_constraint == "constant_volume" else "enthalpy",
            operation.heat_transfer.mode == "isothermal",
            request.relative_tolerance,
        )
        metrics = _metrics(request, actual_initial, final)
    else:
        if len(operation.feeds) != 1 or not feed or sum(feed.values()) <= 0:
            raise ValueError("CSTR/PFR simulation currently requires exactly one positive gas feed")
        stream = operation.feeds[0]
        feed_temperature = stream.temperature.as_unit("K")
        feed_pressure = stream.pressure.as_unit("Pa")
        temperature_mismatch = abs(feed_temperature - temperature_K) / max(temperature_K, 1.0)
        pressure_mismatch = abs(feed_pressure - pressure_Pa) / max(pressure_Pa, 1.0)
        if max(temperature_mismatch, pressure_mismatch) > request.inventory_state_tolerance:
            raise ValueError(
                "The current CSTR/PFR runner requires feed and reactor temperature/pressure to match; "
                "a separate preheater, compressor, or mixer duty is not silently assumed"
            )
        gas.TPX = feed_temperature, feed_pressure, feed
        molecular_weights = dict(zip(gas.species_names, gas.molecular_weights, strict=True))
        mass_flow = sum(amount * molecular_weights[name] / 1000.0 for name, amount in feed.items())
        inlet_energy_flow = mass_flow * gas.enthalpy_mass
        residence_time = operation.residence_time.as_unit("s")
        if operation.reactor_type == "pfr":
            reactor = ct.IdealGasConstPressureMoleReactor(gas, energy=energy_mode, clone=True)
            reactor.volume = mass_flow * residence_time / gas.density
            network = ct.ReactorNet([reactor])
            network.rtol, network.atol = request.relative_tolerance, request.absolute_tolerance
            coordinate = np.linspace(0.0, residence_time, request.points)
            for time_s in coordinate:
                network.advance(float(time_s))
                samples.append(_snapshot(reactor, time_s / residence_time, "residence_fraction", mass_flow))
            final = samples[-1]["species_molar_flows_mol_s"]
            material_balance = _balance(feed, final, gas, "steady molar flow")
            outlet_energy_flow = mass_flow * reactor.phase.enthalpy_mass
            energy_balance = _energy_balance(
                inlet_energy_flow,
                outlet_energy_flow,
                "enthalpy flow",
                operation.heat_transfer.mode == "isothermal",
                request.relative_tolerance,
            )
            metrics = _metrics(request, feed, final)
        else:
            reactor = ct.IdealGasMoleReactor(gas, energy=energy_mode, clone=True, volume=volume_m3)
            inlet_residence_time = reactor.mass / mass_flow
            if abs(inlet_residence_time - residence_time) / residence_time > request.inventory_state_tolerance:
                raise ValueError(
                    "CSTR feed, volume, pressure, temperature, and residence time are inconsistent; "
                    f"inlet-state residence is {inlet_residence_time:.6g} s"
                )
            upstream = ct.Reservoir(gas, clone=True)
            downstream = ct.Reservoir(gas, clone=True)
            inlet = ct.MassFlowController(upstream, reactor, mdot=mass_flow)
            ct.PressureController(reactor, downstream, primary=inlet, K=1e-3)
            network = ct.ReactorNet([reactor])
            network.rtol, network.atol = request.relative_tolerance, request.absolute_tolerance
            network.advance_to_steady_state()
            steady_mass_flow = mass_flow
            actual_residence_time = reactor.mass / steady_mass_flow
            samples.append(_snapshot(reactor, actual_residence_time, "residence_time_s", steady_mass_flow))
            final = samples[-1]["species_molar_flows_mol_s"]
            material_balance = _balance(feed, final, gas, "steady molar flow")
            outlet_energy_flow = steady_mass_flow * reactor.phase.enthalpy_mass
            energy_balance = _energy_balance(
                inlet_energy_flow,
                outlet_energy_flow,
                "steady enthalpy flow",
                operation.heat_transfer.mode == "isothermal",
                request.relative_tolerance,
            )
            metrics = _metrics(request, feed, final)

    solver_stats = network.solver_stats if network is not None else {}
    accepted = bool(material_balance["closed"] and energy_balance["closed"])
    return {
        "status": "completed" if accepted else "balance_failure",
        "accepted": accepted,
        "reactor_type": operation.reactor_type,
        "provenance": provenance,
        "operation_validation": validation,
        "samples": samples,
        "metrics": metrics,
        "material_balance": material_balance,
        "energy_balance": energy_balance,
        "solver": {"name": "Cantera ReactorNet / SUNDIALS", "statistics": dict(solver_stats)},
        "boundary": "Post-discovery operation simulation; mechanism data never feed blind candidate generation.",
        "warnings": [
            "Cantera solves the submitted mechanism; it does not establish that the mechanism is complete or experimentally accurate.",
            "The declared phase-1 mechanism may be a subset of the Cantera mechanism and is used for conservation/yield ceilings only.",
            "PFR uses Cantera's Lagrangian constant-pressure particle formulation parameterized by residence time.",
        ],
    }
