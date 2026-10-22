import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.multiphase.gas_liquid import simulate_gas_liquid_transfer
from app.multiphase import reaktoro_bridge
from app.multiphase.schema import GasLiquidTransferRequest


def gas_liquid_payload(mode: str = "closed_headspace") -> dict:
    payload = {
        "solute_id": "CO2",
        "mode": mode,
        "temperature": {"value": 298.15, "unit": "K"},
        "total_pressure": {"value": 1, "unit": "atm"},
        "gas_volume": {"value": 1, "unit": "liter"},
        "liquid_volume": {"value": 1, "unit": "liter"},
        "duration": {"value": 100, "unit": "s"},
        "points": 21,
        "initial_gas_amount": {"value": 0.001, "unit": "mol"},
        "initial_dissolved_amount": {"value": 0, "unit": "mol"},
        "henry_hcp": {
            "value": 3.3e-4,
            "unit": "mol/(m^3*Pa)",
            "source": "temperature-specific measured Henry coefficient fixture",
            "method": "declared Hcp convention",
            "basis": "experimental",
            "uncertainty": {"value": 0.2e-4, "unit": "mol/(m^3*Pa)"},
        },
        "kla": {
            "value": 0.1,
            "unit": "1/s",
            "source": "measured vessel transfer coefficient fixture",
            "method": "dynamic transfer fit",
            "basis": "experimental",
        },
    }
    if mode == "open_fixed_partial_pressure":
        payload["ambient_partial_pressure"] = {"value": 400, "unit": "Pa"}
    return payload


def test_closed_headspace_transfer_conserves_solute() -> None:
    result = simulate_gas_liquid_transfer(GasLiquidTransferRequest.model_validate(gas_liquid_payload()))
    assert result["accepted"] is True
    assert result["mass_balance"]["closed"] is True
    assert result["samples"][-1]["dissolved_mol"] > 0
    assert len(result["uncertainty_scenarios"]) == 2
    total = result["samples"][-1]["dissolved_mol"] + result["samples"][-1]["headspace_mol"]
    assert total == pytest.approx(0.001)


def test_open_transfer_accounts_for_external_exchange() -> None:
    payload = gas_liquid_payload("open_fixed_partial_pressure")
    payload["initial_gas_amount"] = {"value": 0, "unit": "mol"}
    result = simulate_gas_liquid_transfer(GasLiquidTransferRequest.model_validate(payload))
    assert result["mass_balance"]["closed"] is True
    assert result["samples"][-1]["external_exchange_mol"] > 0


def test_multiphase_api_exposes_honest_reaktoro_unavailability(monkeypatch) -> None:
    monkeypatch.delenv("CHEM_SANDBOX_REAKTORO_URL", raising=False)
    client = TestClient(app)
    capability = client.get("/api/multiphase/capabilities")
    assert capability.status_code == 200
    assert capability.json()["native_gas_liquid_transfer"] is True
    assert capability.json()["electrolyte_equilibrium"]["available"] is False
    request = {
        "database_kind": "phreeqc",
        "database_name": "pitzer.dat",
        "aqueous_species": ["H2O", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2", "CO3-2"],
        "activity_model": "pitzer",
        "temperature": {"value": 298.15, "unit": "K"},
        "pressure": {"value": 1, "unit": "atm"},
        "components": [{"species": "H2O", "amount": {"value": 1, "unit": "kg"}}],
        "source": "explicit post-discovery electrolyte request",
    }
    response = client.post("/api/multiphase/electrolyte/equilibrate", json=request)
    assert response.status_code == 503
    assert "not configured" in response.text


def test_gas_liquid_api() -> None:
    response = TestClient(app).post("/api/multiphase/gas-liquid", json=gas_liquid_payload())
    assert response.status_code == 200, response.text
    assert response.json()["accepted"] is True


def test_remote_reaktoro_validation_error_remains_client_error(monkeypatch) -> None:
    request = {
        "database_kind": "phreeqc",
        "database_name": "pitzer.dat",
        "aqueous_species": ["H2O", "H+", "OH-"],
        "activity_model": "pitzer",
        "temperature": {"value": 298.15, "unit": "K"},
        "pressure": {"value": 1, "unit": "atm"},
        "components": [{"species": "H2O", "amount": {"value": 1, "unit": "kg"}}],
        "source": "invalid-species response propagation fixture",
    }

    class RejectedResponse:
        status_code = 422
        text = '{"detail":"unknown species"}'

        @staticmethod
        def json():
            return {"detail": "unknown species"}

    monkeypatch.setattr(
        reaktoro_bridge,
        "reaktoro_status",
        lambda: {"available": True, "url": "http://127.0.0.1:8011"},
    )
    monkeypatch.setattr(reaktoro_bridge.httpx, "post", lambda *args, **kwargs: RejectedResponse())
    response = TestClient(app).post("/api/multiphase/electrolyte/equilibrate", json=request)
    assert response.status_code == 422
    assert "unknown species" in response.text
