import json

import pytest

from evaluation.protocol import freeze_prediction_artifact, unlock_evaluation
from scripts.run_reaction_benchmarks import prediction_only_payload


def test_evaluation_references_unlock_only_after_prediction_freeze(tmp_path) -> None:
    frozen_path = tmp_path / "prediction.frozen.json"
    evaluation_path = tmp_path / "evaluation.json"
    frozen = freeze_prediction_artifact(
        {"run_id": "blind-1", "predicted_species": ["state-a"], "score": -1.25},
        frozen_path,
    )
    unlocked = unlock_evaluation(
        frozen_path,
        {"expected_product": "kept outside runtime", "literature_yield": 0.72},
        evaluation_path,
    )
    assert unlocked["prediction_sha256"] == frozen["prediction_sha256"]
    assert unlocked["status"] == "references_unlocked_post_run"
    assert json.loads(evaluation_path.read_text(encoding="utf-8"))["reference_sha256"]


def test_prediction_freeze_rejects_reference_leakage(tmp_path) -> None:
    with pytest.raises(ValueError, match="Reference fields are forbidden"):
        freeze_prediction_artifact(
            {"states": [], "expected_product": "leak"},
            tmp_path / "must-not-exist.json",
        )


def test_unlock_detects_modified_prediction(tmp_path) -> None:
    path = tmp_path / "prediction.frozen.json"
    freeze_prediction_artifact({"score": 2.0}, path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["prediction"]["score"] = 3.0
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="altered"):
        unlock_evaluation(path, {"answer": 2.5}, tmp_path / "evaluation.json")


def test_benchmark_projection_removes_reference_scoring_fields(tmp_path) -> None:
    payload = {
        "cases": [
            {
                "id": "case-1",
                "stoichiometry": {"A": -1, "B": 1},
                "reference_kj_mol": -10,
                "reference_note": "locked",
                "references": ["https://example.invalid"],
            }
        ],
        "results": [
            {
                "reaction_id": "case-1",
                "predicted_kj_mol": -12,
                "reference_kj_mol": -10,
                "raw_deviation_kj_mol": -2,
                "sign_agreement": True,
            }
        ],
    }
    projected = prediction_only_payload(payload)
    frozen = freeze_prediction_artifact(projected, tmp_path / "prediction.frozen.json")
    assert frozen["prediction"]["results"][0] == {"reaction_id": "case-1", "predicted_kj_mol": -12}
