"""Compute seven assay-free feed predictions, then evaluate ethane bars."""
from __future__ import annotations

import json
from typing import Any

from scripts.audit_dry_feed_sweep import (
    FIGURE_1A_METHANE_FEED_FRACTIONS,
    audit_dry_feed_sweep,
)


def run_heldout_ethane_plot_check() -> dict[str, Any]:
    sweep = audit_dry_feed_sweep(
        methane_feed_fractions=FIGURE_1A_METHANE_FEED_FRACTIONS,
    )
    from evaluation.reference_only.dry_ethane_plot_bands import (
        evaluate_frozen_dry_ethane_plot_bands,
    )

    return evaluate_frozen_dry_ethane_plot_bands(sweep)


if __name__ == "__main__":
    print(json.dumps(run_heldout_ethane_plot_check(), indent=2))
