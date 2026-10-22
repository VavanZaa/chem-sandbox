import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError

from app.main import app
from app.operations.schema import OperationValidationRequest
from app.operations.validation import validate_operation


def operation_payload() -> dict:
    return {
        "operation": {
            "id": "batch-validation",
            "reactor_type": "batch",
            "volume": {"value": 1, "unit": "liter"},
            "temperature": {"value": 25, "unit": "degC"},
            "pressure": {"value": 1, "unit": "atm"},
            "duration": {"value": 2, "unit": "hour"},
            "phases": [
                {
                    "id": "liquid",
                    "kind": "ideal_solution",
                    "species_ids": ["A", "B"],
                    "activity_model": "ideal molarity",
                }
            ],
            "initial_inventory": [
                {"species_id": "A", "phase_id": "liquid", "amount": {"value": 0.1, "unit": "mol"}}
            ],
            "heat_transfer": {
                "mode": "isothermal",
                "setpoint_temperature": {"value": 25, "unit": "degC"},
            },
        },
        "mechanism": {
            "id": "isomerization",
            "species": [
                {"id": "A", "name": "A", "elements": {"C": 1, "H": 4}},
                {"id": "B", "name": "B", "elements": {"C": 1, "H": 4}},
            ],
            "reactions": [
                {
                    "id": "r1",
                    "stoichiometry": {"A": -1, "B": 1},
                    "evidence_source": "user supplied elementary mechanism",
                }
            ],
        },
    }


def test_operation_units_normalize_and_mechanism_balances() -> None:
    result = validate_operation(OperationValidationRequest.model_validate(operation_payload()))
    assert result["ready_for_solver"] is True
    assert result["normalized_si"]["volume_m3"] == pytest.approx(0.001)
    assert result["normalized_si"]["temperature_K"] == pytest.approx(298.15)
    assert result["normalized_si"]["duration_s"] == pytest.approx(7200)
    assert result["mechanism_balance"]["all_balanced"] is True
    assert result["energy_balance"]["status"] == "requires_thermochemistry_and_reactor_solution"


def test_unbalanced_mechanism_is_not_ready_for_solver() -> None:
    payload = operation_payload()
    payload["mechanism"]["species"][1]["elements"] = {"C": 1, "H": 2}
    result = validate_operation(OperationValidationRequest.model_validate(payload))
    assert result["ready_for_solver"] is False
    assert result["mechanism_balance"]["reactions"][0]["element_residual"] == {"H": -2.0}


def test_reactor_requirements_and_dimensions_are_strict() -> None:
    payload = operation_payload()
    payload["operation"]["reactor_type"] = "cstr"
    with pytest.raises(ValidationError, match="residence_time"):
        OperationValidationRequest.model_validate(payload)
    payload = operation_payload()
    payload["operation"]["volume"] = {"value": 1, "unit": "second"}
    response = TestClient(app).post("/api/operations/validate", json=payload)
    assert response.status_code == 422
    assert "not compatible with m^3" in response.text


def test_operation_validation_api() -> None:
    response = TestClient(app).post("/api/operations/validate", json=operation_payload())
    assert response.status_code == 200
    assert response.json()["status"] == "valid"


def test_nonphysical_operation_values_are_rejected() -> None:
    payload = operation_payload()
    payload["operation"]["volume"] = {"value": -1, "unit": "liter"}
    response = TestClient(app).post("/api/operations/validate", json=payload)
    assert response.status_code == 422
    assert "greater than zero" in response.text
