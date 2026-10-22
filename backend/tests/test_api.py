from fastapi.testclient import TestClient
from ase import Atoms
from ase.io import write

from app.main import app
from app.api import routes
from app.experiment.schema import default_experiment


def test_health_and_default_definition() -> None:
    client = TestClient(app)
    assert client.get("/api/health").json()["status"] == "ok"
    definition = client.get("/api/experiments/default")
    assert definition.status_code == 200
    assert len(definition.json()["substances"]) == 3
    assert definition.json()["knowledge_policy"]["mode"] == "blind_physics_only"


def test_scientific_policy_and_operation_capabilities() -> None:
    client = TestClient(app)
    policy = client.get("/api/scientific-policy")
    assert policy.status_code == 200
    assert policy.json()["knowledge_assisted_mode_available"] is False
    capabilities = client.get("/api/operations/capabilities")
    assert capabilities.status_code == 200
    installed = {item["id"]: item for item in capabilities.json()["installed"]}
    assert installed["cantera"]["available"] is True
    assert installed["idaes"]["available"] is True
    assert installed["libsbml"]["available"] is True
    assert installed["roadrunner"]["available"] is True


def test_saved_structure_endpoint_returns_atoms_and_detected_bonds(tmp_path, monkeypatch) -> None:
    runs = tmp_path / "runs"
    structures = runs / "test_job" / "structures"
    structures.mkdir(parents=True)
    write(
        structures / "water.xyz",
        Atoms("OH2", positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]),
    )
    monkeypatch.setattr(routes, "RUNS_DIR", runs)
    client = TestClient(app)
    response = client.get("/api/experiments/test_job/structures/water.xyz")
    assert response.status_code == 200
    payload = response.json()
    assert len(payload["atoms"]) == 3
    assert len(payload["bonds"]) == 2
    assert payload["species"][0]["name"] == "Water"


def test_stoichiometry_api_without_model():
    response = TestClient(app).post("/api/estimates/stoichiometry", json=default_experiment().model_dump(mode="json"))
    assert response.status_code == 200
    assert response.json()["co2_mol"] == .05


def test_persisted_result_and_kinetics_validation(tmp_path, monkeypatch):
    import json
    runs = tmp_path / "runs"
    job = runs / "saved"
    job.mkdir(parents=True)
    result = {"graph": {"nodes": [{"id": "A"}, {"id": "B"}], "edges": [{"id": "ab", "source": "A", "target": "B"}]}}
    (job / "result.json").write_text(json.dumps(result), encoding="utf-8")
    monkeypatch.setattr(routes, "RUNS_DIR", runs)
    client = TestClient(app)
    assert client.get("/api/experiments/saved/result").json() == result
    request = {"initial_state_id": "A", "duration_s": 1, "acknowledge_conditional_model": True, "rates": []}
    assert client.post("/api/experiments/saved/kinetics", json=request).status_code == 422
    request["rates"] = [{"edge_id": "ab", "rate_s_inv": 1, "temperature_K": 298.15, "source": "hypothetical"}]
    response = client.post("/api/experiments/saved/kinetics", json=request)
    assert response.status_code == 200
    assert response.json()["actual_yield"] is None
    assert client.get("/api/experiments/missing/result").status_code == 404
