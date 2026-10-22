import json

from app.experiment.schema import default_experiment
from app.search.search_runner import SearchRunner


def test_search_runner_serializes_logs(tmp_path, deterministic_calculator) -> None:
    experiment = default_experiment()
    experiment.search.runs = 2
    experiment.search.optimization_steps = 1
    runner = SearchRunner(deterministic_calculator, output_dir=tmp_path / "run")
    result = runner.run(experiment)
    assert result["status"] == "completed"
    assert result["statistics"]["search_runs"] == 2
    assert result["barriers_calculated"] is False
    assert (tmp_path / "run" / "experiment.json").exists()
    assert (tmp_path / "run" / "events.jsonl").exists()
    saved = json.loads((tmp_path / "run" / "result.json").read_text(encoding="utf-8"))
    def all_keys(value):
        if isinstance(value, dict):
            for key, nested in value.items():
                yield key
                yield from all_keys(nested)
        elif isinstance(value, list):
            for nested in value:
                yield from all_keys(nested)

    assert "reaction_probability" not in set(all_keys(saved))
    assert "Neither is a yield" in saved["statistics"]["disclaimer"]
    assert all(node["interpretation_quality"] == "resolved" for node in saved["graph"]["nodes"])
    assert len(saved["run_outcomes"]) == 2
    assert all(outcome["endpoints"] for outcome in saved["run_outcomes"])
    assert saved["statistics"]["accounted_runs"] == 2
    assert all(outcome["evaluations"] <= experiment.search.max_evaluations_per_run for outcome in saved["run_outcomes"])
    assert saved["statistics"]["rejected_candidate_states"] >= 0
    assert len(saved["graph"]["edges"]) >= 1


def test_initial_model_failures_are_all_accounted(tmp_path, deterministic_calculator, monkeypatch):
    experiment = default_experiment()
    experiment.search.runs = 3
    def fail(*args, **kwargs):
        raise RuntimeError("test model offline")
    monkeypatch.setattr(deterministic_calculator, "optimize", fail)
    result = SearchRunner(deterministic_calculator, output_dir=tmp_path / "failed").run(experiment)
    assert result["status"] == "model_failure"
    assert result["statistics"]["accounted_runs"] == 3
    assert result["statistics"]["run_status_counts"] == {"model_failure": 3}
    assert result["statistics"]["failed_optimizations"] == 3
    assert all(item["endpoints"][0]["state_id"] is None for item in result["run_outcomes"])


def test_unsupported_structures_are_saved_for_diagnostics(tmp_path, deterministic_calculator):
    from ase import Atoms
    runner = SearchRunner(deterministic_calculator, output_dir=tmp_path / "unsupported")
    runner.output_dir.mkdir()
    payload = runner._evaluate(Atoms("H", positions=[[0, 0, 0]]), run_index=0, name="isolated",
                               steps=1, fmax=.35, method="test")
    assert payload is None
    assert runner.last_failure == "unsupported_connectivity"
    assert runner.errors[0]["structure_reference"]
    assert (runner.output_dir / "structures" / "run_000_isolated.xyz").is_file()
