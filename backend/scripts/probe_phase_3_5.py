from __future__ import annotations

import json
import math

from app.multiphase.gas_liquid import simulate_gas_liquid_transfer
from app.multiphase.schema import GasLiquidTransferRequest
from app.operations.cantera_runner import simulate_cantera
from app.operations.cantera_schema import CanteraSimulationRequest
from app.pathways.rates import build_rate_artifact
from app.pathways.schema import RateArtifactRequest


def energy(value: float) -> dict:
    return {
        "value": value,
        "unit": "kJ/mol",
        "source": "phase-3 validation input generated before references are unlocked",
        "method": "declared thermochemical validation datum",
        "basis": "computed",
        "uncertainty": {"value": 1.0, "unit": "kJ/mol"},
    }


def cantera_payload(reactor_type: str) -> dict:
    temperature = 1200.0
    pressure = 101325.0
    total_moles = 4.0
    volume = total_moles * 8.31446261815324 * temperature / pressure
    operation = {
        "id": f"validation-{reactor_type}",
        "reactor_type": reactor_type,
        "volume": {"value": volume, "unit": "m^3"},
        "temperature": {"value": temperature, "unit": "K"},
        "pressure": {"value": pressure, "unit": "Pa"},
        "phases": [{
            "id": "gas",
            "kind": "gas",
            "species_ids": ["H2", "O2", "H2O", "AR"],
            "activity_model": "Cantera ideal gas",
        }],
        "heat_transfer": {"mode": "adiabatic"},
    }
    if reactor_type == "batch":
        operation["duration"] = {"value": 0.005, "unit": "s"}
        operation["initial_inventory"] = [
            {"species_id": "H2", "phase_id": "gas", "amount": {"value": 2.0, "unit": "mol"}},
            {"species_id": "O2", "phase_id": "gas", "amount": {"value": 1.0, "unit": "mol"}},
            {"species_id": "AR", "phase_id": "gas", "amount": {"value": 1.0, "unit": "mol"}},
        ]
    else:
        residence_time = 0.005
        scale = 1e-3
        operation["residence_time"] = {"value": residence_time, "unit": "s"}
        operation["feeds"] = [{
            "id": "feed",
            "phase_id": "gas",
            "temperature": {"value": temperature, "unit": "K"},
            "pressure": {"value": pressure, "unit": "Pa"},
            "species_flows": [
                {"species_id": "H2", "flow": {"value": 2 * scale, "unit": "mol/s"}},
                {"species_id": "O2", "flow": {"value": scale, "unit": "mol/s"}},
                {"species_id": "AR", "flow": {"value": scale, "unit": "mol/s"}},
            ],
        }]
        if reactor_type == "cstr":
            operation["volume"] = {
                "value": sum((2 * scale, scale, scale)) * residence_time * 8.31446261815324 * temperature / pressure,
                "unit": "m^3",
            }
    return {
        "definition": {
            "operation": operation,
            "mechanism": {
                "id": "hydrogen-oxidation-validation-subset",
                "species": [
                    {"id": "H2", "name": "Hydrogen", "elements": {"H": 2}},
                    {"id": "O2", "name": "Oxygen", "elements": {"O": 2}},
                    {"id": "H2O", "name": "Water", "elements": {"H": 2, "O": 1}},
                    {"id": "AR", "name": "Argon", "elements": {"Ar": 1}},
                ],
                "reactions": [{
                    "id": "overall_balance_only",
                    "stoichiometry": {"H2": -1, "O2": -0.5, "H2O": 1},
                    "evidence_source": "balance-only subset; bundled Cantera YAML supplies kinetics",
                }],
            },
        },
        "cantera": {
            "builtin": "h2o2.yaml",
            "source": "Cantera 3.2 bundled hydrogen mechanism",
            "allow_declared_subset_of_cantera_mechanism": True,
        },
        "points": 21,
        "reference_species": "H2",
        "target_species": ["H2O"],
    }


def main() -> None:
    rate = build_rate_artifact(RateArtifactRequest.model_validate({
        "reaction_id": "validation_A_to_B",
        "temperature": {"value": 500.0, "unit": "K"},
        "pressure": {"value": 1.0, "unit": "atm"},
        "forward_order": 1,
        "reverse_order": 1,
        "forward_activation_gibbs": energy(80.0),
        "reverse_activation_gibbs": energy(90.0),
        "reaction_gibbs": energy(-10.0),
    }))
    cantera_results = {}
    for reactor in ("batch", "pfr", "cstr"):
        try:
            cantera_results[reactor] = simulate_cantera(
                CanteraSimulationRequest.model_validate(cantera_payload(reactor))
            )
        except Exception as error:
            cantera_results[reactor] = {"error": f"{type(error).__name__}: {error}"}
    gas_payload = {
        "solute_id": "validation-solute",
        "mode": "closed_headspace",
        "temperature": {"value": 298.15, "unit": "K"},
        "total_pressure": {"value": 2.0, "unit": "atm"},
        "gas_volume": {"value": 1.0, "unit": "liter"},
        "liquid_volume": {"value": 1.0, "unit": "liter"},
        "duration": {"value": 1000.0, "unit": "s"},
        "points": 11,
        "initial_gas_amount": {"value": 0.001, "unit": "mol"},
        "initial_dissolved_amount": {"value": 0.0, "unit": "mol"},
        "henry_hcp": {
            "value": 3.3e-4,
            "unit": "mol/(m^3*Pa)",
            "source": "validation-only measured coefficient fixture",
            "method": "declared Hcp convention",
            "basis": "experimental",
        },
        "kla": {
            "value": 0.1,
            "unit": "1/s",
            "source": "validation-only transfer coefficient fixture",
            "method": "declared dynamic transfer coefficient",
            "basis": "experimental",
        },
    }
    gas_liquid = simulate_gas_liquid_transfer(GasLiquidTransferRequest.model_validate(gas_payload))
    summary = {
        "rate": {
            "certified": rate["certified"],
            "kf": rate["forward"]["rate_constant"],
            "kr": rate["reverse"]["rate_constant"],
            "ratio": rate["forward"]["rate_constant"] / rate["reverse"]["rate_constant"],
            "keq": rate["thermodynamics"]["equilibrium_constant_dimensionless"],
        },
        "cantera": {
            name: ({
                "accepted": result["accepted"],
                "conversion": result["metrics"]["conversion"],
                "water_yield": result["metrics"]["targets"]["H2O"]["yield_fraction"],
                "final_temperature_K": result["samples"][-1]["temperature_K"],
                "material_residual": result["material_balance"]["max_relative_element_residual"],
                "energy_balance": result["energy_balance"],
            } if "error" not in result else result)
            for name, result in cantera_results.items()
        },
        "gas_liquid": {
            "accepted": gas_liquid["accepted"],
            "final": gas_liquid["samples"][-1],
            "equilibrium_residual": gas_liquid["final_equilibrium_residual_mol_m3"],
        },
    }
    print(json.dumps(summary, indent=2, default=lambda value: value.item()))


if __name__ == "__main__":
    main()
