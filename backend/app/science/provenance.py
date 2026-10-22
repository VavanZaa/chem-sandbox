from __future__ import annotations

import platform
import sys
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from typing import Any

from app.experiment.schema import ExperimentDefinition


TRACKED_DISTRIBUTIONS = (
    "ase",
    "aimnet",
    "Cantera",
    "idaes-pse",
    "mace-torch",
    "numpy",
    "rdkit",
    "scipy",
    "torch",
    "torchani",
)


def _package_versions() -> dict[str, str | None]:
    result: dict[str, str | None] = {}
    for distribution in TRACKED_DISTRIBUTIONS:
        try:
            result[distribution] = version(distribution)
        except PackageNotFoundError:
            result[distribution] = None
    return result


def build_runtime_provenance(
    experiment: ExperimentDefinition,
    model_records: list[dict[str, Any]],
) -> dict[str, Any]:
    hardware: dict[str, Any] = {"cuda_available": False, "gpu": None}
    try:
        import torch

        hardware["cuda_available"] = bool(torch.cuda.is_available())
        hardware["cuda_runtime"] = torch.version.cuda
        if torch.cuda.is_available():
            hardware["gpu"] = torch.cuda.get_device_name(0)
            hardware["gpu_capability"] = list(torch.cuda.get_device_capability(0))
    except ImportError:
        hardware["torch"] = "not installed"
    checkpoint_records = [
        {
            "model": record.get("model"),
            "file": record.get("checkpoint_file"),
            "sha256": record.get("checkpoint_sha256"),
        }
        for record in model_records
    ]
    all_checkpoint_hashes = bool(checkpoint_records) and all(item["sha256"] for item in checkpoint_records)
    return {
        "schema_version": "1.0",
        "captured_at": datetime.now(UTC).isoformat(),
        "application": {"name": "chem-sandbox", "version": "0.26.0"},
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": _package_versions(),
        "hardware": hardware,
        "models": model_records,
        "requested_engine": experiment.compute.engine,
        "requested_device": experiment.compute.device,
        "search_seed": experiment.search.seed,
        "search_settings": experiment.search.model_dump(mode="json"),
        "checkpoint_hashes": {
            "status": "complete" if all_checkpoint_hashes else "incomplete",
            "records": checkpoint_records,
            "requirement": "Every executed model must expose a checkpoint SHA-256 for a fully locked run.",
        },
    }
