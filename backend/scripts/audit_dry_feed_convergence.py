"""Assay-free temporal-resolution audit of the dry-feed trend.

This compares discretizations of the same uniform-packet scenario. It does
not resolve individual streamer pulses or validate their space-time exposure.
"""
from __future__ import annotations

import json
from typing import Any

from scripts.audit_dry_feed_sweep import METHANE_FEED_FRACTIONS, audit_dry_feed_sweep


def audit_dry_feed_convergence(
    *, substep_counts: tuple[int, ...] = (16, 32, 64),
) -> dict[str, Any]:
    if len(substep_counts) < 2 or any(
        not isinstance(count, int) or count < 1 or count > 128
        for count in substep_counts
    ) or any(a >= b for a, b in zip(substep_counts, substep_counts[1:])):
        raise ValueError("Substep counts must be strictly increasing integers in [1, 128]")

    sweeps = [audit_dry_feed_sweep(substeps=count) for count in substep_counts]
    rows = []
    for index, methane_fraction in enumerate(METHANE_FEED_FRACTIONS):
        conversions = [
            sweep["rows"][index]["ch4_conversion_fraction"] for sweep in sweeps
        ]
        ethane_selectivities = [
            sweep["rows"][index]["c2h6_carbon_selectivity_fraction"]
            for sweep in sweeps
        ]
        hydrogen_selectivities = [
            sweep["rows"][index]["h2_hydrogen_selectivity_fraction"]
            for sweep in sweeps
        ]
        rows.append({
            "feed_CH4_mole_fraction": methane_fraction,
            "ch4_conversion_by_substeps": dict(zip(substep_counts, conversions)),
            "ethane_carbon_selectivity_by_substeps": dict(
                zip(substep_counts, ethane_selectivities)
            ),
            "hydrogen_selectivity_by_substeps": dict(
                zip(substep_counts, hydrogen_selectivities)
            ),
            "last_pair_ch4_conversion_absolute_change": abs(
                conversions[-1] - conversions[-2]
            ),
            "last_pair_ethane_selectivity_absolute_change": abs(
                ethane_selectivities[-1] - ethane_selectivities[-2]
            ),
            "last_pair_hydrogen_selectivity_absolute_change": abs(
                hydrogen_selectivities[-1] - hydrogen_selectivities[-2]
            ),
        })
    return {
        "status": "same_scenario_temporal_discretization_audit_not_experimental_validation",
        "substep_counts": substep_counts,
        "rows": rows,
        "ch4_conversion_feed_trend_by_substeps": {
            count: (
                sweep["rows"][-1]["ch4_conversion_fraction"]
                - sweep["rows"][0]["ch4_conversion_fraction"]
            )
            for count, sweep in zip(substep_counts, sweeps)
        },
        "all_material_ledgers_closed": all(
            row["material_ledger_closed"] for sweep in sweeps for row in sweep["rows"]
        ),
        "all_collision_energy_upper_bounds_closed": all(
            row["collision_energy_upper_bound_closed"]
            for sweep in sweeps for row in sweep["rows"]
        ),
        "all_n2_product_state_sources_excluded": all(
            row["n2_product_state_source_excluded"]
            for sweep in sweeps for row in sweep["rows"]
        ),
        "measured_product_assay_loaded": False,
        "streamer_pulses_resolved": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_feed_convergence(), indent=2))
