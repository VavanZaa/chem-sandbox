"""Feed-only numerical screen for simultaneous gas chemistry and water contact.

The 2024 HCN assay is not loaded. The contact coefficient and spark exposure
are hypothetical, so no experimental yield validation is possible here.
"""

from __future__ import annotations

import argparse
import json
import sys

from scripts.audit_barth_feed_convergence import CONFIGURATIONS, summarize_convergence
from scripts.probe_barth_feed_pathway import run_probe


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gas-water-kla-per-s", type=float, default=1e-4)
    parser.add_argument("--primary-source-forcing", choices=("time_averaged", "inventory_limited"),
                        default="time_averaged")
    parser.add_argument("--formaldehyde-hydration-model", choices=(
        "off", "winkelman_2002_source_bounded",
        "winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    ), default="off")
    arguments = parser.parse_args()
    reports = []
    failures = []
    for cycles, substeps in CONFIGURATIONS:
        print(f"Running coupled contact: {cycles * substeps} one-hour steps", file=sys.stderr, flush=True)
        try:
            report = run_probe(
                duration_s=3600, cycles=cycles, substeps=substeps,
                gas_water_kla_per_s=arguments.gas_water_kla_per_s,
                gas_water_solver="coupled_bdf",
                primary_source_forcing=arguments.primary_source_forcing,
                formaldehyde_hydration_model=arguments.formaldehyde_hydration_model,
            )
        except Exception as error:
            failures.append({
                "steps": cycles * substeps,
                "error": f"{type(error).__name__}: {error}",
                "numerical_result_accepted": False,
            })
            continue
        if (
            report["status"] != "feed_matched_short_pathway_diagnostic_not_yield_validation"
            or report["charged_branch_failed_substeps"] != 0
        ):
            failures.append({
                "steps": cycles * substeps,
                "status": report["status"],
                "charged_branch_failed_substeps": report["charged_branch_failed_substeps"],
                "charged_branch_failure_examples": report["charged_branch_failure_examples"],
                "numerical_result_accepted": False,
            })
        else:
            reports.append(report)
    if failures:
        print(json.dumps({
            "status": "coupled_contact_numerical_screen_failed",
            "failed_resolutions": failures,
            "accepted_resolution_count": len(reports),
            "experimental_hcn_amount_loaded": False,
            "science_boundary": "No HCN amount is accepted from a failed resolution screen.",
        }, indent=2))
        raise SystemExit(1)
    summary = summarize_convergence(reports)
    summary["gas_water_contact_solver"] = "coupled_bdf"
    summary["primary_source_forcing"] = arguments.primary_source_forcing
    summary["hypothetical_liquid_side_kla_per_s"] = arguments.gas_water_kla_per_s
    summary["formaldehyde_hydration_model"] = arguments.formaldehyde_hydration_model
    summary["experimental_hcn_amount_loaded"] = False
    print(json.dumps(summary, indent=2))
    if summary["status"] != "one_hour_feed_only_numerical_screen_passed_not_yield_validation":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
