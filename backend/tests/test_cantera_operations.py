import math

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.operations.cantera_runner import simulate_cantera
from app.operations.cantera_schema import CanteraSimulationRequest


def cantera_payload(reactor_type: str = "batch") -> dict:
    temperature = 300.0
    pressure = 101325.0
    inventory_mol = 7.0
    volume = inventory_mol * 8.31446261815324 * temperature / pressure
    operation = {
        "id": f"{reactor_type}-cantera-test",
        "reactor_type": reactor_type,
        "volume": {"value": volume, "unit": "m^3"},
        "temperature": {"value": temperature, "unit": "K"},
        "pressure": {"value": pressure, "unit": "Pa"},
        "phases": [
            {
                "id": "gas",
                "kind": "gas",
                "species_ids": ["H2", "O2", "H2O", "N2"],
                "activity_model": "Cantera ideal gas",
            }
        ],
        "heat_transfer": {"mode": "adiabatic"},
    }
    if reactor_type == "batch":
        operation["duration"] = {"value": 1e-6, "unit": "s"}
        operation["initial_inventory"] = [
            {"species_id": "H2", "phase_id": "gas", "amount": {"value": 2, "unit": "mol"}},
            {"species_id": "O2", "phase_id": "gas", "amount": {"value": 1, "unit": "mol"}},
            {"species_id": "N2", "phase_id": "gas", "amount": {"value": 4, "unit": "mol"}},
        ]
    else:
        residence_time = 1.0 if reactor_type == "cstr" else 1e-6
        operation["residence_time"] = {"value": residence_time, "unit": "s"}
        if reactor_type == "cstr":
            operation["volume"] = {"value": 0.007 * 8.31446261815324 * temperature / pressure, "unit": "m^3"}
        operation["feeds"] = [
            {
                "id": "feed",
                "phase_id": "gas",
                "temperature": {"value": temperature, "unit": "K"},
                "pressure": {"value": pressure, "unit": "Pa"},
                "species_flows": [
                    {"species_id": "H2", "flow": {"value": 0.002, "unit": "mol/s"}},
                    {"species_id": "O2", "flow": {"value": 0.001, "unit": "mol/s"}},
                    {"species_id": "N2", "flow": {"value": 0.004, "unit": "mol/s"}},
                ],
            }
        ]
    return {
        "definition": {
            "operation": operation,
            "mechanism": {
                "id": "declared-hydrogen-oxidation-subset",
                "species": [
                    {"id": "H2", "name": "Hydrogen", "elements": {"H": 2}},
                    {"id": "O2", "name": "Oxygen", "elements": {"O": 2}},
                    {"id": "H2O", "name": "Water", "elements": {"H": 2, "O": 1}},
                    {"id": "N2", "name": "Nitrogen", "elements": {"N": 2}},
                ],
                "reactions": [
                    {
                        "id": "overall_test_balance_only",
                        "stoichiometry": {"H2": -1, "O2": -0.5, "H2O": 1},
                        "evidence_source": "Declared balance subset; Cantera YAML is authoritative for kinetics",
                    }
                ],
            },
        },
        "cantera": {
            "builtin": "h2o2.yaml",
            "source": "Cantera 3.2 bundled hydrogen mechanism",
            "allow_declared_subset_of_cantera_mechanism": True,
        },
        "points": 3,
        "reference_species": "H2",
        "target_species": ["H2O"],
    }


def test_cantera_batch_closes_material_and_energy_balances() -> None:
    result = simulate_cantera(CanteraSimulationRequest.model_validate(cantera_payload()))
    assert result["status"] == "completed"
    assert result["accepted"] is True
    assert result["material_balance"]["closed"] is True
    assert result["energy_balance"]["closed"] is True
    assert len(result["samples"]) == 3
    assert len(result["provenance"]["mechanism_sha256"]) == 64
    assert math.isfinite(result["metrics"]["conversion"])


def test_cantera_pfr_lagrangian_path() -> None:
    result = simulate_cantera(CanteraSimulationRequest.model_validate(cantera_payload("pfr")))
    assert result["accepted"] is True
    assert result["samples"][-1]["residence_fraction"] == 1.0
    assert result["material_balance"]["closed"] is True


def test_cantera_cstr_steady_state_path() -> None:
    result = simulate_cantera(CanteraSimulationRequest.model_validate(cantera_payload("cstr")))
    assert result["accepted"] is True
    assert result["samples"][0]["residence_time_s"] == pytest.approx(1.0, rel=0.05)
    assert result["material_balance"]["closed"] is True


def test_cantera_operation_api() -> None:
    response = TestClient(app).post("/api/operations/cantera/simulate", json=cantera_payload())
    assert response.status_code == 200, response.text
    assert response.json()["boundary"].startswith("Post-discovery")
