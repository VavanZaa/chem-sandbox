import math
import pytest
from fastapi.testclient import TestClient
from pydantic import ValidationError
from app.experiment.efficiency import EfficiencyRequest, simulate_efficiency


def barrier(rate, temperature=298.15):
    return -8.31446261815324 * temperature * math.log(rate * 6.62607015e-34 / (1.380649e-23 * temperature)) / 1000


def request():
    return EfficiencyRequest.model_validate({
        "species": [{"id": "A", "name": "Ethanol", "smiles": "CCO", "initial_mol": .05},
                    {"id": "B", "name": "Ether", "smiles": "COC", "initial_mol": 0},
                    {"id": "C", "name": "Ethene", "smiles": "C=C", "initial_mol": 0},
                    {"id": "W", "name": "Water", "smiles": "O", "initial_mol": 0}],
        "models": [{"id": "test", "basis": "hypothetical", "channels": [
            {"id": "ab", "reactants": {"A": 1}, "products": {"B": 1}, "activation_free_energy_kj_mol": barrier(2), "source": "analytic test"},
            {"id": "ac", "reactants": {"A": 1}, "products": {"C": 1, "W": 1}, "activation_free_energy_kj_mol": barrier(1), "source": "analytic test"}]}],
        "reference_reactant": "A", "target_product": "B", "product_pool": ["B", "C"],
        "duration_s": 1, "sensitivity_kj_mol": 0, "acknowledge_assumptions": True,
    })


def test_competing_yields_match_analytic_and_conserve_atoms():
    result = simulate_efficiency(request())
    model = result["models"][0]
    nominal = model["nominal"]
    assert nominal["conversion"] == pytest.approx(1 - math.exp(-3), rel=1e-6)
    assert nominal["target_yield"] == pytest.approx(2 / 3 * (1 - math.exp(-3)), rel=1e-6)
    assert nominal["target_product_share"] == pytest.approx(2 / 3)
    assert model["stoichiometric_ceiling_mol"] == pytest.approx(.05)
    assert result["calibrated_product_probabilities"] is None
    for sample in model["samples"]:
        n = sample["moles"]
        assert n["A"] + n["B"] + n["C"] == pytest.approx(.05)
        assert n["W"] == pytest.approx(n["C"])


def test_time_temperature_barrier_and_sensitivity_change_conversion():
    req = request()
    nominal = simulate_efficiency(req)["models"][0]["nominal"]["conversion"]
    req.temperature_K = 330
    assert simulate_efficiency(req)["models"][0]["nominal"]["conversion"] > nominal
    req.temperature_K = 298.15
    req.duration_s = .1
    assert simulate_efficiency(req)["models"][0]["nominal"]["conversion"] < nominal
    req.sensitivity_kj_mol = 5
    model = simulate_efficiency(req)["models"][0]
    assert len(model["scenarios"]) == 5
    assert model["yield_sensitivity_range"][0] < model["nominal"]["target_yield"] < model["yield_sensitivity_range"][1]


def test_bimolecular_concentration_and_gas_pressure_are_used():
    req = request()
    req.species[0].initial_mol = 0
    req.species[2].initial_mol = req.species[3].initial_mol = .05
    req.reference_reactant, req.target_product, req.product_pool = "C", "A", ["A"]
    req.models[0].channels = [req.models[0].channels[1]]
    channel = req.models[0].channels[0]
    channel.reactants, channel.products = channel.products, channel.reactants
    req.volume_L = .05
    # C + W -> A, k=1 L/mol/s: equimolar half-conversion at C0=1 M, t=1 s.
    assert simulate_efficiency(req)["models"][0]["nominal"]["conversion"] == pytest.approx(.5, rel=1e-6)
    req.volume_L = .5
    assert simulate_efficiency(req)["models"][0]["nominal"]["conversion"] == pytest.approx(1 / 11, rel=1e-6)
    req.phase = "ideal_gas"
    low = simulate_efficiency(req)["models"][0]["nominal"]["conversion"]
    req.pressure_atm = 5
    assert simulate_efficiency(req)["models"][0]["nominal"]["conversion"] > low


def test_bad_channels_and_incomplete_engine_comparison_rejected():
    req = request()
    req.models[0].channels[1].products.pop("W")
    with pytest.raises(ValueError, match="conserve"):
        simulate_efficiency(req)
    req = request()
    other = req.models[0].model_copy(deep=True)
    other.id = "incomplete"
    other.channels.pop()
    req.models.append(other)
    with pytest.raises(ValueError, match="identical reaction channels"):
        simulate_efficiency(req)
    req = request()
    req.models[0].channels[0].products = {"UNKNOWN": 1}
    with pytest.raises(ValueError, match="Unknown species"):
        simulate_efficiency(req)


def test_no_artificial_efficiency_cap_and_api_validation():
    req = request()
    req.models[0].channels = req.models[0].channels[:1]
    req.duration_s = 60
    result = simulate_efficiency(req)
    assert result["models"][0]["nominal"]["target_yield"] == pytest.approx(1, abs=1e-6)
    from app.main import app
    client = TestClient(app)
    assert client.post("/api/efficiency", json=req.model_dump()).status_code == 200
    data = req.model_dump()
    data["acknowledge_assumptions"] = False
    assert client.post("/api/efficiency", json=data).status_code == 422
    data["temperature_K"] = float("nan")
    with pytest.raises(ValidationError):
        EfficiencyRequest.model_validate(data)
