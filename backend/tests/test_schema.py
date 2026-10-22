import pytest
from pydantic import ValidationError

from app.experiment.compiler import ExperimentCompiler
from app.experiment.schema import ExperimentDefinition, default_experiment


def test_default_experiment_round_trip() -> None:
    experiment = default_experiment()
    assert ExperimentDefinition.model_validate_json(experiment.model_dump_json()) == experiment
    assert experiment.environment.temperature_K == 298.15
    assert {substance.unit for substance in experiment.substances} == {"mol"}


def test_v01_rejects_missing_substance() -> None:
    payload = default_experiment().model_dump(mode="json")
    payload["substances"].pop()
    with pytest.raises(ValidationError, match="requires baking soda"):
        ExperimentDefinition.model_validate(payload)


def test_rejects_wrong_unit() -> None:
    payload = default_experiment().model_dump(mode="json")
    payload["substances"][0]["unit"] = "mL"
    with pytest.raises(ValidationError, match="supports mol, g"):
        ExperimentDefinition.model_validate(payload)


def test_accepts_editable_mass_and_molar_units() -> None:
    payload = default_experiment().model_dump(mode="json")
    payload["substances"][0].update(amount=0.1, unit="mol")
    payload["substances"][1].update(amount=3.0, unit="g")
    experiment = ExperimentDefinition.model_validate(payload)
    assert experiment.substances[0].unit == "mol"
    assert experiment.substances[1].unit == "g"


def test_relaxation_budget_is_explicit_not_silently_capped() -> None:
    experiment = default_experiment()
    experiment.search.effort = "high"
    experiment.search.runs = 12
    experiment.search.optimization_steps = 16
    plan = ExperimentCompiler().compile(experiment)
    assert plan.runs == 12
    assert plan.optimization_steps == 16
    assert plan.candidate_optimization_steps == 6
    assert any("user-configurable" in assumption for assumption in plan.assumptions)
