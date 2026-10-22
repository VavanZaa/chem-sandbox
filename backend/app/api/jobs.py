from __future__ import annotations

import threading
import traceback
import uuid
import json
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from app.experiment.schema import ExperimentDefinition
from app.models.registry import compatible_engines, make_engine
from app.search.search_runner import SearchRunner
from app.search.consensus import build_discovery_comparison
from app.science.blind_policy import build_blind_manifest
from app.science.provenance import build_runtime_provenance


@dataclass(slots=True)
class ExperimentJob:
    id: str
    status: str = "queued"
    stage: str = "queued"
    progress: int = 0
    message: str = "Waiting for compute slot…"
    created_at: str = field(default_factory=lambda: datetime.now(UTC).isoformat())
    updated_at: str = field(default_factory=lambda: datetime.now(UTC).isoformat())
    result: dict[str, Any] | None = None
    error: str | None = None

    def payload(self) -> dict[str, Any]:
        return asdict(self)


class JobManager:
    def __init__(self, runs_dir: Path) -> None:
        self.runs_dir = runs_dir
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self.jobs: dict[str, ExperimentJob] = {}
        self.lock = threading.RLock()
        self.compute_lock = threading.Lock()
        self.engines = {}

    def engine(self, engine_id: str, device: str):
        key = (engine_id, device)
        if key not in self.engines:
            self.engines[key] = make_engine(engine_id, device)
        return self.engines[key]

    def release_engines(self, engine_ids, device):
        for engine_id in engine_ids:
            self.engines.pop((engine_id, device), None)
        try:
            import torch
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        except ImportError:
            pass

    def create(self) -> ExperimentJob:
        job = ExperimentJob(id=uuid.uuid4().hex[:12])
        with self.lock:
            self.jobs[job.id] = job
        return job

    def get(self, job_id: str) -> ExperimentJob | None:
        with self.lock:
            return self.jobs.get(job_id)

    def _update(self, job: ExperimentJob, stage: str, progress: int, message: str) -> None:
        with self.lock:
            job.stage = stage
            job.progress = progress
            job.message = message
            job.updated_at = datetime.now(UTC).isoformat()

    def run(self, job_id: str, experiment: ExperimentDefinition) -> None:
        job = self.get(job_id)
        if job is None:
            return
        with self.compute_lock:
            try:
                blind_manifest = build_blind_manifest(experiment)
                with self.lock:
                    job.status = "running"
                if experiment.compute.engine == "all_compatible":
                    support = compatible_engines(experiment)
                    selected = [item["engine"] for item in support if item["supported"]]
                    excluded = [item for item in support if not item["supported"]]
                    if not selected:
                        raise ValueError("No installed engine supports this submitted atom/charge/spin inventory")
                    completed, failures = [], []
                    for index, engine_id in enumerate(selected):
                        concrete = experiment.model_copy(deep=True)
                        concrete.compute.engine = engine_id
                        result_job_id = job_id if index == 0 else f"{job_id}--{engine_id}"
                        start = int(index * 96 / len(selected))
                        span = max(1, int(96 / len(selected)))
                        try:
                            runner = SearchRunner(
                                self.engine(engine_id, experiment.compute.device),
                                output_dir=self.runs_dir / result_job_id,
                                progress=lambda stage, progress, message, i=engine_id, s=start, w=span:
                                    self._update(job, f"{i}:{stage}", min(96, s + int(w * progress / 100)),
                                                 f"{i} · {message}"),
                            )
                            engine_result = runner.run(concrete)
                            if engine_result["status"] != "completed":
                                raise RuntimeError("engine produced no valid graph")
                            completed.append((engine_id, engine_result, result_job_id))
                        except Exception as error:
                            failures.append({"engine": engine_id, "status": "failed",
                                             "reason": f"{type(error).__name__}: {error}"})
                        finally:
                            # Four resident models can consume substantial VRAM;
                            # all-compatible mode deliberately trades reload time
                            # for predictable memory use.
                            self.release_engines([engine_id], experiment.compute.device)
                    if not completed:
                        raise RuntimeError("All compatible engine searches failed")
                    # The first compatible engine remains the displayed graph so
                    # existing geometry/NEB routes retain stable references.
                    result = completed[0][1]
                    result["experiment"] = experiment.model_dump(mode="json")
                    result["discovery_comparison"] = build_discovery_comparison(completed, excluded, failures)
                    result["display_engine"] = completed[0][0]
                    (self.runs_dir / job_id / "result.json").write_text(
                        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
                    )
                else:
                    runner = SearchRunner(
                        self.engine(experiment.compute.engine, experiment.compute.device),
                        output_dir=self.runs_dir / job_id,
                        progress=lambda stage, progress, message: self._update(
                            job, stage, progress, message
                        ),
                    )
                    result = runner.run(experiment)
                result["blind_manifest"] = blind_manifest
                model_records = [item[1].get("model", {}) for item in completed] if experiment.compute.engine == "all_compatible" else [result.get("model", {})]
                result["provenance"] = build_runtime_provenance(experiment, model_records)
                result_path = self.runs_dir / job_id / "result.json"
                result_path.parent.mkdir(parents=True, exist_ok=True)
                result_path.write_text(
                    json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
                )
                with self.lock:
                    job.result = result
                    job.status = "completed" if result["status"] == "completed" else "failed"
                    job.stage = "complete" if job.status == "completed" else "model_failure"
                    job.progress = 100
                    job.message = (
                        (f"Compared {result['discovery_comparison']['successful_engine_count']} compatible engine searches"
                         if result.get("discovery_comparison") else "Reaction network ready")
                        if job.status == "completed"
                        else "The model produced no valid scientific result"
                    )
                    job.updated_at = datetime.now(UTC).isoformat()
            except Exception as error:
                with self.lock:
                    job.status = "failed"
                    job.stage = "failed"
                    job.error = f"{type(error).__name__}: {error}"
                    job.message = "Model failure / result uncertain"
                    job.updated_at = datetime.now(UTC).isoformat()
                failure_dir = self.runs_dir / job_id
                failure_dir.mkdir(parents=True, exist_ok=True)
                (failure_dir / "failure.txt").write_text(traceback.format_exc(), encoding="utf-8")

    def refine(self, job_id, parent_result, settings):
        from app.search.refinement import compare_engines
        job = self.get(job_id)
        if job is None:
            return
        with self.compute_lock:
            job.status = "running"
            try:
                job.result = compare_engines(parent_result, settings, self.engine, self.runs_dir,
                                            self.runs_dir / job_id,
                                            lambda stage, percent, message: self._update(job, stage, percent, message))
                job.status = "completed"
                self._update(job, "complete", 100, "Engine comparison complete; inspect convergence and domain warnings")
            except Exception as error:
                job.status = "failed"
                job.error = f"{type(error).__name__}: {error}"
                self._update(job, "failed", 100, "Path comparison failed")
            finally:
                self.release_engines(settings.engines, settings.device)
