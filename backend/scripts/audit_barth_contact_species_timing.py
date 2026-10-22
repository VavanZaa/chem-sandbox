"""Diagnostic-only test of contact-species sensitivity in the Barth feed probe.

This deliberately removes selected Henry entries for one controlled comparison.
It does not define a new predictive property set, load an HCN assay, or
validate any product amount. The unmodified property data are restored after
each run because the override exists only inside this script's context.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
from unittest.mock import patch

from app.prebiotic import phase_partition
from scripts.probe_barth_feed_pathway import run_probe


def compare_contact_timing(*, excluded_species: tuple[str, ...],
                           cycles: int = 1, substeps: int = 48) -> dict:
    reference = phase_partition.load_henry_subset()
    unknown = set(excluded_species) - set(reference["entries"])
    if unknown:
        raise ValueError(f"Unknown Henry entry in diagnostic exclusion: {sorted(unknown)}")
    restricted = deepcopy(reference)
    for species in excluded_species:
        del restricted["entries"][species]
    reports = {}
    with patch.object(phase_partition, "load_henry_subset", return_value=restricted):
        for solver in ("post_step_analytic", "coupled_bdf"):
            reports[solver] = run_probe(
                duration_s=3600, cycles=cycles, substeps=substeps,
                gas_water_kla_per_s=1e-4, gas_water_solver=solver,
            )
    accepted = all(
        report["status"] == "feed_matched_short_pathway_diagnostic_not_yield_validation"
        and report["material_ledger_closed"] is True
        and report["charged_branch_material_ledger_closed"] is True
        and report["charged_branch_failed_substeps"] == 0
        for report in reports.values()
    )
    if not accepted:
        return {
            "status": "contact_species_timing_diagnostic_failed",
            "excluded_henry_species": list(excluded_species),
            "branch_status": {name: item["status"] for name, item in reports.items()},
            "numerical_amount_accepted": False,
            "experimental_hcn_amount_loaded": False,
        }
    post = reports["post_step_analytic"]["final_aqueous_HCN_plus_CN_minus_mol"]
    coupled = reports["coupled_bdf"]["final_aqueous_HCN_plus_CN_minus_mol"]
    return {
        "status": "diagnostic_only_not_yield_validation",
        "excluded_henry_species": list(excluded_species),
        "duration_s": 3600,
        "integration_steps": cycles * substeps,
        "hypothetical_kla_per_s": 1e-4,
        "post_step_aqueous_hcn_plus_cn_mol": post,
        "coupled_aqueous_hcn_plus_cn_mol": coupled,
        "post_step_to_coupled_ratio": post / coupled if coupled > 0 else None,
        "material_and_charged_controls_pass": True,
        "experimental_hcn_amount_loaded": False,
        "science_boundary": (
            "Removed Henry entries are a mechanism-timing diagnostic, not a "
            "defensible property model or an experimental HCN comparison."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exclude", nargs="+", required=True)
    parser.add_argument("--cycles", type=int, choices=(1, 2, 4), default=1)
    parser.add_argument("--substeps", type=int, choices=(48,), default=48)
    arguments = parser.parse_args()
    result = compare_contact_timing(
        excluded_species=tuple(arguments.exclude),
        cycles=arguments.cycles, substeps=arguments.substeps,
    )
    print(json.dumps(result, indent=2))
    if result["status"] != "diagnostic_only_not_yield_validation":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
