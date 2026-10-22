import math
import pytest
from pydantic import ValidationError
from app.experiment.kinetics import KineticsRequest, RateEvidence, simulate_kinetics
from app.experiment.schema import default_experiment
from app.experiment.stoichiometry import theoretical_yield
from app.graph.graph_statistics import build_statistics


def test_stoichiometry_defaults_and_unit_equivalence():
    experiment = default_experiment()
    first = theoretical_yield(experiment)
    assert first["limiting_reagent"] == "Acetic acid"
    assert first["co2_mol"] == .05
    assert first["co2_mass_g"] == pytest.approx(2.200475)
    assert first["excess_bicarbonate_mol"] == pytest.approx(.00952)
    assert first["co2_ideal_gas_L"] == pytest.approx(1.22327, rel=1e-5)
    experiment.substances[0].amount *= 84.0066
    experiment.substances[0].unit = "g"
    experiment.substances[1].amount *= 60.052
    experiment.substances[1].unit = "g"
    second = theoretical_yield(experiment)
    assert second["co2_mol"] == pytest.approx(first["co2_mol"])
    experiment.substances[0].amount = .84
    assert theoretical_yield(experiment)["limiting_reagent"] == "Baking soda"


GRAPH = {"nodes": [{"id": i} for i in "ABC"], "edges": [
    {"id": "ab", "source": "A", "target": "B"}, {"id": "ac", "source": "A", "target": "C"}]}


def request():
    return KineticsRequest(initial_state_id="A", duration_s=1, acknowledge_conditional_model=True,
                           rates=[RateEvidence(edge_id="ab", rate_s_inv=2, temperature_K=298.15, source="test assumption"),
                                  RateEvidence(edge_id="ac", rate_s_inv=1, temperature_K=298.15, source="test assumption")])


def test_competing_products_analytic_solution_and_conservation():
    result = simulate_kinetics(GRAPH, request())
    final = result["samples"][-1]["fractions"]
    assert final["A"] == pytest.approx(math.exp(-3))
    assert final["B"] == pytest.approx(2 / 3 * (1 - math.exp(-3)))
    assert final["C"] == pytest.approx(1 / 3 * (1 - math.exp(-3)))
    assert all(sum(sample["fractions"].values()) == pytest.approx(1) for sample in result["samples"])
    assert result["actual_yield"] is None


def test_reversible_kinetics_requires_independent_reverse_rate():
    graph = {"nodes": [{"id": "A"}, {"id": "B"}], "edges": [GRAPH["edges"][0], {"id": "ba", "source": "B", "target": "A"}]}
    req = request()
    req.rates[1].edge_id = "ba"
    req.duration_s = 10
    result = simulate_kinetics(graph, req)
    assert result["samples"][-1]["fractions"]["B"] == pytest.approx(2 / 3)


def test_missing_duplicate_or_wrong_temperature_rates_rejected():
    req = request()
    req.rates.pop()
    with pytest.raises(ValueError, match="every directed"):
        simulate_kinetics(GRAPH, req)
    req = request()
    req.rates[1].temperature_K = 300
    with pytest.raises(ValueError, match="temperatures"):
        simulate_kinetics(GRAPH, req)
    req = request()
    req.rates.append(req.rates[0])
    with pytest.raises(ValueError, match="Duplicate"):
        simulate_kinetics(GRAPH, req)


def test_eyring_and_invalid_parameters():
    req = request()
    rate = req.rates[0]
    rate.rate_s_inv = None
    rate.activation_free_energy_ev = .8
    result = simulate_kinetics(GRAPH, req)
    expected = 1.380649e-23 * 298.15 / 6.62607015e-34 * math.exp(-.8 / (8.617333262145e-5 * 298.15))
    assert result["rates"][0]["computed_rate_s_inv"] == pytest.approx(expected)
    for bad in [float("nan"), float("inf"), -1]:
        with pytest.raises(ValidationError):
            RateEvidence(edge_id="ab", rate_s_inv=bad, temperature_K=298.15, source="test")


def test_co2_runs_use_union_not_sum_or_max():
    nodes = [{"formula": "CO2", "run_indices": [0, 1], "endpoint_run_indices": [0]},
             {"formula": "CO2 + H2O", "run_indices": [1, 2], "endpoint_run_indices": [0, 2]}]
    stats = build_statistics({"nodes": nodes, "edges": []}, total_runs=3, valid_states=4, failed_optimizations=0)
    assert stats["co2_containing_runs"] == 3
    assert stats["co2_terminal_runs"] == 2
