"""Run four assay-free residence conditions before checking Figure 1D bands."""
from __future__ import annotations

import json
from typing import Any

from scripts.audit_dry_residence_sweep import audit_dry_residence_sweep


def run_heldout_residence_plot_check() -> dict[str, Any]:
    sweep = audit_dry_residence_sweep()
    from evaluation.reference_only.dry_conversion_plot_bands import (
        evaluate_frozen_dry_residence_plot_bands,
    )

    return evaluate_frozen_dry_residence_plot_bands(sweep)


if __name__ == "__main__":
    print(json.dumps(run_heldout_residence_plot_check(), indent=2))
