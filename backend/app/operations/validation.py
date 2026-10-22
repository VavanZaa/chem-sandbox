from __future__ import annotations

from typing import Any

from .schema import MechanismDefinition, OperationValidationRequest


BALANCE_TOLERANCE = 1e-10


def mechanism_balance_report(mechanism: MechanismDefinition) -> dict[str, Any]:
    species = {item.id: item for item in mechanism.species}
    rows = []
    all_balanced = True
    for reaction in mechanism.reactions:
        element_delta: dict[str, float] = {}
        charge_delta = 0.0
        site_delta = 0.0
        for species_id, coefficient in reaction.stoichiometry.items():
            item = species[species_id]
            for element, count in item.elements.items():
                element_delta[element] = element_delta.get(element, 0.0) + coefficient * count
            charge_delta += coefficient * item.charge
            site_delta += coefficient * item.surface_sites
        element_delta = {key: value for key, value in element_delta.items() if abs(value) > BALANCE_TOLERANCE}
        balanced = not element_delta and abs(charge_delta) <= BALANCE_TOLERANCE and abs(site_delta) <= BALANCE_TOLERANCE
        all_balanced &= balanced
        rows.append(
            {
                "reaction_id": reaction.id,
                "balanced": balanced,
                "element_residual": element_delta,
                "charge_residual": charge_delta,
                "surface_site_residual": site_delta,
            }
        )
    return {"all_balanced": all_balanced, "tolerance": BALANCE_TOLERANCE, "reactions": rows}


def validate_operation(request: OperationValidationRequest) -> dict[str, Any]:
    operation = request.operation
    balance = mechanism_balance_report(request.mechanism)
    volume_m3 = operation.volume.as_unit("m^3")
    temperature_K = operation.temperature.as_unit("K")
    pressure_Pa = operation.pressure.as_unit("Pa")
    duration_s = None if operation.duration is None else operation.duration.as_unit("s")
    residence_time_s = None if operation.residence_time is None else operation.residence_time.as_unit("s")
    positive = {
        "volume": volume_m3,
        "temperature": temperature_K,
        "pressure": pressure_Pa,
        "duration": duration_s,
        "residence_time": residence_time_s,
    }
    invalid_positive = [name for name, value in positive.items() if value is not None and value <= 0]
    if invalid_positive:
        raise ValueError(f"Operation quantities must be greater than zero: {', '.join(invalid_positive)}")
    inventory = []
    for item in operation.initial_inventory:
        amount_mol = item.amount.as_unit("mol")
        if amount_mol < 0:
            raise ValueError("Initial species amounts cannot be negative")
        inventory.append({"species_id": item.species_id, "phase_id": item.phase_id, "amount_mol": amount_mol})
    feeds = []
    for feed in operation.feeds:
        start_time_s = feed.start_time.as_unit("s")
        end_time_s = None if feed.end_time is None else feed.end_time.as_unit("s")
        feed_temperature_K = feed.temperature.as_unit("K")
        feed_pressure_Pa = feed.pressure.as_unit("Pa")
        if start_time_s < 0 or (end_time_s is not None and end_time_s <= start_time_s):
            raise ValueError(f"Feed {feed.id} requires start >= 0 and end > start")
        if feed_temperature_K <= 0 or feed_pressure_Pa <= 0:
            raise ValueError(f"Feed {feed.id} temperature and pressure must be greater than zero")
        flow_values = {item.species_id: item.flow.as_unit("mol/s") for item in feed.species_flows}
        if any(value < 0 for value in flow_values.values()) or not any(value > 0 for value in flow_values.values()):
            raise ValueError(f"Feed {feed.id} requires nonnegative flows and at least one positive flow")
        feeds.append(
            {
                "id": feed.id,
                "phase_id": feed.phase_id,
                "temperature_K": feed_temperature_K,
                "pressure_Pa": feed_pressure_Pa,
                "start_time_s": start_time_s,
                "end_time_s": end_time_s,
                "species_molar_flows_mol_s": flow_values,
            }
        )
    normalized = {
        "id": operation.id,
        "reactor_type": operation.reactor_type,
        "volume_m3": volume_m3,
        "temperature_K": temperature_K,
        "pressure_Pa": pressure_Pa,
        "duration_s": duration_s,
        "residence_time_s": residence_time_s,
        "initial_inventory_mol": inventory,
        "feeds": feeds,
    }
    if operation.heat_transfer.setpoint_temperature:
        normalized["heat_setpoint_K"] = operation.heat_transfer.setpoint_temperature.as_unit("K")
        if normalized["heat_setpoint_K"] <= 0:
            raise ValueError("Heat setpoint must be greater than 0 K")
    if operation.heat_transfer.ambient_temperature:
        normalized["ambient_temperature_K"] = operation.heat_transfer.ambient_temperature.as_unit("K")
        if normalized["ambient_temperature_K"] <= 0:
            raise ValueError("Ambient temperature must be greater than 0 K")
    if operation.heat_transfer.overall_UA:
        normalized["overall_UA_W_K"] = operation.heat_transfer.overall_UA.as_unit("W/K")
        if normalized["overall_UA_W_K"] < 0:
            raise ValueError("overall_UA cannot be negative")
    if operation.mixing.speed:
        normalized["mixing_speed_s_inv"] = operation.mixing.speed.as_unit("1/s")
        if normalized["mixing_speed_s_inv"] < 0:
            raise ValueError("Mixing speed cannot be negative")
    if operation.mixing.dispersion_coefficient:
        normalized["dispersion_m2_s"] = operation.mixing.dispersion_coefficient.as_unit("m^2/s")
        if normalized["dispersion_m2_s"] < 0:
            raise ValueError("Dispersion coefficient cannot be negative")
    if operation.gas_transfer.kla:
        normalized["kla_s_inv"] = operation.gas_transfer.kla.as_unit("1/s")
        if normalized["kla_s_inv"] < 0:
            raise ValueError("kLa cannot be negative")
    if operation.catalyst.loading:
        normalized["catalyst_loading_mol"] = operation.catalyst.loading.as_unit("mol")
        if normalized["catalyst_loading_mol"] < 0:
            raise ValueError("Catalyst loading cannot be negative")
    if operation.catalyst.active_site_amount:
        normalized["active_site_amount_mol"] = operation.catalyst.active_site_amount.as_unit("mol")
        if normalized["active_site_amount_mol"] < 0:
            raise ValueError("Active-site amount cannot be negative")
    return {
        "status": "valid" if balance["all_balanced"] else "invalid_mechanism_balance",
        "ready_for_solver": balance["all_balanced"],
        "normalized_si": normalized,
        "mechanism_balance": balance,
        "energy_balance": {
            "status": "requires_thermochemistry_and_reactor_solution",
            "note": "Input validation cannot claim an energy balance before species enthalpies and a converged reactor trajectory exist.",
        },
        "solver_note": "Validation does not run kinetics or predict products/yield.",
    }
