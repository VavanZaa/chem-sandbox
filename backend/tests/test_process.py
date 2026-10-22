import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.process.flowsheet import simulate_process_flowsheet
from app.process.schema import ProcessFlowsheetRequest


def fraction(value: float, uncertainty: float | None = None) -> dict:
    payload = {
        "value": value,
        "unit": "dimensionless",
        "source": "post-discovery process measurement fixture",
        "method": "steady-state process measurement",
        "basis": "experimental",
    }
    if uncertainty is not None:
        payload["uncertainty"] = {"value": uncertainty, "unit": "dimensionless"}
    return payload


def process_payload(optimize: bool = False) -> dict:
    return {
        "id": "single-recycle-validation",
        "components": [
            {"id": "A", "name": "Reactant A", "elements": {"C": 1}},
            {"id": "B", "name": "Product B", "elements": {"C": 1}},
        ],
        "fresh_feed": {
            "A": {"value": 1.0, "unit": "mol/s"},
            "B": {"value": 0.0, "unit": "mol/s"},
        },
        "reaction": {
            "id": "A_to_B",
            "stoichiometry": {"A": -1, "B": 1},
            "limiting_species": "A",
            "evidence_source": "balanced post-discovery reaction artifact",
        },
        "single_pass_conversion": fraction(0.5, 0.02),
        "product_recoveries": {"A": fraction(0.0), "B": fraction(1.0, 0.01)},
        "recycle": {
            "fraction": fraction(0.5, 0.02),
            "optimize": optimize,
            "lower_bound": 0.0,
            "upper_bound": 0.9,
            "grid_points": 31,
        },
        "target_species": "B",
        "utility_per_extent": {
            "value": 10.0,
            "unit": "kJ/mol",
            "source": "measured specific reactor utility fixture",
            "method": "steady-state energy balance",
            "basis": "experimental",
        },
        "max_internal_molar_flow": {"value": 2.0, "unit": "mol/s"},
    }


def test_recycle_flowsheet_matches_analytic_solution_and_closes_balance() -> None:
    result = simulate_process_flowsheet(ProcessFlowsheetRequest.model_validate(process_payload()))
    flowsheet = result["flowsheet"]
    assert result["accepted"] is True
    assert flowsheet["mixed_feed_mol_s"]["A"] == pytest.approx(4 / 3, rel=1e-8)
    assert flowsheet["product_mol_s"]["B"] == pytest.approx(2 / 3, rel=1e-8)
    assert flowsheet["material_balance"]["closed"] is True
    assert len(result["uncertainty_scenarios"]) == 2
    if result["idaes_crosscheck"]["capabilities"]["ipopt_available"]:
        assert result["idaes_crosscheck"]["status"] == "completed", result["idaes_crosscheck"]
        assert result["cross_engine_max_relative_product_difference"] < 1e-7


def test_recycle_optimization_respects_internal_flow_constraint() -> None:
    result = simulate_process_flowsheet(ProcessFlowsheetRequest.model_validate(process_payload(True)))
    assert result["accepted"] is True
    assert result["optimization"]["enabled"] is True
    assert result["flowsheet"]["internal_molar_flow_mol_s"] <= 2.0 * (1 + 1e-10)
    assert result["flowsheet"]["target_product_mol_s"] > 2 / 3


def test_process_api_exposes_capability_and_simulation() -> None:
    client = TestClient(app)
    capability = client.get("/api/process/capabilities")
    assert capability.status_code == 200
    assert "idaes_available" in capability.json()
    response = client.post("/api/process/simulate", json=process_payload())
    assert response.status_code == 200, response.text
    assert response.json()["accepted"] is True
