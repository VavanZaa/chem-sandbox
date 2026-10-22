import pytest

from app.pathways.rates import build_rate_artifact
from app.pathways.schema import RateArtifactRequest


def energy(value: float, uncertainty: float = 2.0) -> dict:
    return {
        "value": value,
        "unit": "kJ/mol",
        "source": "certified thermochemistry artifact sha256:example",
        "method": "consistent electronic and thermal method",
        "basis": "computed",
        "uncertainty": {"value": uncertainty, "unit": "kJ/mol"},
    }


def rate_payload() -> dict:
    return {
        "reaction_id": "A_to_B",
        "temperature": {"value": 298.15, "unit": "K"},
        "pressure": {"value": 1, "unit": "atm"},
        "forward_order": 1,
        "reverse_order": 1,
        "forward_activation_gibbs": energy(60),
        "reverse_activation_gibbs": energy(50),
        "reaction_gibbs": energy(10),
        "impact_weight": 0.8,
    }


def test_rate_artifact_enforces_detailed_balance_and_uncertainty() -> None:
    result = build_rate_artifact(RateArtifactRequest.model_validate(rate_payload()))
    assert result["certified"] is True
    assert result["operation_ready"] is True
    assert result["forward"]["rate_constant"] > 0
    assert result["forward"]["uncertainty_range"][0] < result["forward"]["rate_constant"]
    assert result["thermodynamics"]["cycle_residual_kj_mol"] == pytest.approx(0)
    assert result["active_learning"]["priority_score"] > 0


def test_inconsistent_reverse_barrier_refuses_rate_claim() -> None:
    payload = rate_payload()
    payload["reverse_activation_gibbs"]["value"] = 30
    result = build_rate_artifact(RateArtifactRequest.model_validate(payload))
    assert result["certified"] is False
    assert result["forward"]["rate_constant"] is None
    assert any("microscopic reversibility" in item for item in result["problems"])


def test_pressure_interpolation_must_be_bracketed_and_reversible() -> None:
    payload = rate_payload()
    payload["pressure_dependence"] = {
        "model": "tabulated_log_interpolation",
        "source": "pressure-dependent master-equation calculation",
        "basis": "computed",
        "points": [
            {"pressure": {"value": 0.1, "unit": "atm"}, "forward_multiplier": 0.5, "reverse_multiplier": 0.5},
            {"pressure": {"value": 10, "unit": "atm"}, "forward_multiplier": 2.0, "reverse_multiplier": 2.0},
        ],
    }
    result = build_rate_artifact(RateArtifactRequest.model_validate(payload))
    assert result["certified"] is True

    payload["pressure"]["value"] = 100
    result = build_rate_artifact(RateArtifactRequest.model_validate(payload))
    assert result["certified"] is False
    assert any("outside" in item for item in result["problems"])
