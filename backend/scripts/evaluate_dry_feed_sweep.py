"""Compute the dry feed sweep, then read the held-out paper trend."""
from __future__ import annotations

import json

from scripts.audit_dry_feed_sweep import audit_dry_feed_sweep


def run_heldout_feed_check() -> dict[str, object]:
    sweep = audit_dry_feed_sweep()
    from evaluation.reference_only.dry_hcn_2024 import evaluate_frozen_dry_feed_sweep

    return evaluate_frozen_dry_feed_sweep(sweep)


if __name__ == "__main__":
    print(json.dumps(run_heldout_feed_check(), indent=2))
