"""Freeze the assay-free dry probe, then inspect the held-out paper ordering."""
from __future__ import annotations

import json

from scripts.run_dry_flow_hcn_probe import run_dry_flow_probe


def run_heldout_check() -> dict[str, object]:
    prediction = run_dry_flow_probe(source_forcing="evolving_eedf", substeps=16)
    # Import the reference-only module after the model has finished.
    from evaluation.reference_only.dry_hcn_2024 import evaluate_frozen_dry_flow_prediction

    return evaluate_frozen_dry_flow_prediction(prediction)


if __name__ == "__main__":
    print(json.dumps(run_heldout_check(), indent=2))
