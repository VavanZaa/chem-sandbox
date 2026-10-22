import pytest
from fastapi.testclient import TestClient

from app.catalysis.microkinetics import export_catmap_setup, simulate_catalysis
from app.catalysis.schema import CatalysisSimulationRequest
from app.main import app


def energy(value: float) -> dict:
    return {
        "value": value,
        "unit": "kJ/mol",
        "source": "computed surface pathway artifact sha256:test",
        "method": "consistent free-energy protocol",
        "basis": "computed",
        "uncertainty": {"value": 2.0, "unit": "kJ/mol"},
    }


def catalytic_payload() -> dict:
    return {
        "id": "generic-catalytic-cycle",
        "temperature": {"value": 500, "unit": "K"},
        "pressure": {"value": 1, "unit": "atm"},
        "duration": {"value": 0.1, "unit": "s"},
        "points": 21,
        "active_site_amount": {"value": 1e-6, "unit": "mol"},
        "species": [
            {"id": "A", "name": "Gas A", "phase": "gas", "elements": {"C": 1}},
            {"id": "B", "name": "Gas B", "phase": "gas", "elements": {"C": 1}},
            {"id": "star", "name": "Vacant site", "phase": "surface", "site_occupancy": 1, "vacant_site": True},
            {"id": "Astar", "name": "Adsorbed A", "phase": "surface", "elements": {"C": 1}, "site_occupancy": 1},
            {"id": "Bstar", "name": "Adsorbed B", "phase": "surface", "elements": {"C": 1}, "site_occupancy": 1},
        ],
        "reactions": [
            {
                "id": "adsorb_A",
                "process": "adsorption",
                "stoichiometry": {"A": -1, "star": -1, "Astar": 1},
                "forward_activation_gibbs": energy(90),
                "reverse_activation_gibbs": energy(100),
                "reaction_gibbs": energy(-10),
            },
            {
                "id": "convert_surface",
                "process": "surface_reaction",
                "stoichiometry": {"Astar": -1, "Bstar": 1},
                "forward_activation_gibbs": energy(100),
                "reverse_activation_gibbs": energy(100),
                "reaction_gibbs": energy(0),
            },
            {
                "id": "desorb_B",
                "process": "desorption",
                "stoichiometry": {"Bstar": -1, "B": 1, "star": 1},
                "forward_activation_gibbs": energy(100),
                "reverse_activation_gibbs": energy(90),
                "reaction_gibbs": energy(10),
            },
        ],
        "initial_surface_amounts_per_site": {"star": 1.0},
        "gas_activities": {"A": 1.0, "B": 0.0},
        "target_gas_species": ["B"],
        "transport_effectiveness": {
            "value": 0.9,
            "source": "computed effectiveness-factor validation fixture",
            "basis": "computed",
        },
    }


def test_microkinetics_closes_site_balance_and_produces_target() -> None:
    result = simulate_catalysis(CatalysisSimulationRequest.model_validate(catalytic_payload()))
    assert result["accepted"] is True
    assert result["site_balance"]["max_absolute_residual"] < 1e-8
    assert result["samples"][-1]["gas_net_change_mol"]["B"] > 0
    assert result["target_turnover_frequencies_s_inv"]["B"] > 0


def test_catalyst_free_endpoint_cancellation_is_rejected() -> None:
    payload = catalytic_payload()
    payload["reactions"] = [{
        "id": "fake_catalyst_label",
        "process": "surface_reaction",
        "stoichiometry": {"A": -1, "B": 1},
        "forward_activation_gibbs": energy(100),
        "reverse_activation_gibbs": energy(100),
        "reaction_gibbs": energy(0),
    }]
    with pytest.raises(ValueError, match="no bound intermediate"):
        CatalysisSimulationRequest.model_validate(payload)


def test_catmap_export_and_api() -> None:
    request = CatalysisSimulationRequest.model_validate(catalytic_payload())
    exported = export_catmap_setup(request)
    assert "rxn_expressions" in exported["setup_text"]
    assert "A_g + * <-> Astar_s" in exported["reaction_expressions"]
    response = TestClient(app).post("/api/catalysis/simulate", json=catalytic_payload())
    assert response.status_code == 200, response.text
    assert response.json()["accepted"] is True

