import ast
from pathlib import Path

import pytest
from pydantic import ValidationError

from app.experiment.schema import ExperimentDefinition, default_experiment
from app.science.blind_policy import build_blind_manifest


@pytest.mark.parametrize(
    "field",
    ["expected_product", "reaction_recipe", "known_yield", "reference_energy"],
)
def test_target_specific_knowledge_is_rejected_at_api_schema(field: str) -> None:
    payload = default_experiment().model_dump(mode="json")
    payload[field] = "must not enter discovery"
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        ExperimentDefinition.model_validate(payload)


def test_nested_unknown_knowledge_is_not_silently_ignored() -> None:
    payload = default_experiment().model_dump(mode="json")
    payload["search"]["expected_pathway"] = "answer leakage"
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        ExperimentDefinition.model_validate(payload)


def test_blind_manifest_is_stable_and_disclaims_target_knowledge() -> None:
    experiment = default_experiment()
    first = build_blind_manifest(experiment)
    second = build_blind_manifest(experiment)
    assert first == second
    assert len(first["input_sha256"]) == 64
    assert first["target_specific_knowledge_used"] is False
    assert first["reaction_recipe_used"] is False
    assert first["reference_access"] == "post_run_evaluation_only"


def test_runtime_discovery_does_not_import_benchmark_or_evaluation_answers() -> None:
    app_root = Path(__file__).resolve().parents[1] / "app"
    forbidden = ("app.benchmarks", "evaluation", "reference_only")
    violations: list[str] = []
    for path in app_root.rglob("*.py"):
        if "benchmarks" in path.relative_to(app_root).parts:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        imported = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.append(node.module)
        if any(module == token or module.startswith(f"{token}.") for module in imported for token in forbidden):
            violations.append(str(path.relative_to(app_root)))
    assert not violations, f"Reference data crossed into discovery runtime: {violations}"
