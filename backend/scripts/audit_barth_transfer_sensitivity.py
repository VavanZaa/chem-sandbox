"""One-hour feed-only sensitivity to declared gas/water contact kinetics.

The Barth HCN assay is never opened here. All finite kLa values are
hypothetical; this audit is not a calibration or an experimental yield test.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from scripts.probe_barth_feed_pathway import run_probe


DEFAULT_KLA_PER_S = (0.0, 1e-5, 1e-4, 1e-3)


def summarize_transfer_sensitivity(
    reports: list[dict[str, Any]], failures: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    if len(reports) < 2:
        raise ValueError("Transfer sensitivity needs finite-rate and equilibrium reports")
    failures = failures or []
    reference = reports[0]
    same_conditions = all(
        item.get("duration_s") == 3600
        and item.get("spark_on") is True
        and item.get("initial_feed_mol") == reference.get("initial_feed_mol")
        and item.get("energy_basis") == "declared_electrical_estimate"
        and item.get("total_integration_steps") == reference.get("total_integration_steps")
        and item.get("cycles") == reference.get("cycles")
        and item.get("integration_substeps_per_cycle")
        == reference.get("integration_substeps_per_cycle")
        and item.get("gas_water_contact_solver") == "post_step_analytic"
        for item in reports
    ) and "HCN" not in reference.get("initial_feed_mol", {})
    contact_labels_valid = (
        sum(item.get("gas_water_contact_model") == "step_equilibrium" for item in reports) == 1
        and all(
            item.get("gas_water_contact_model") == (
                "step_equilibrium"
                if item.get("hypothetical_liquid_side_kla_per_s") is None
                else "finite_rate_contact"
            )
            for item in reports
        )
    )
    controls_pass = same_conditions and contact_labels_valid and not failures and all(
        item.get("status") == "feed_matched_short_pathway_diagnostic_not_yield_validation"
        and item.get("material_ledger_closed") is True
        and item.get("charged_branch_material_ledger_closed") is True
        and item.get("charged_branch_failed_substeps") == 0
        for item in reports
    )
    scenarios = [
        {
            "gas_water_contact_model": item["gas_water_contact_model"],
            "hypothetical_liquid_side_kla_per_s": item["hypothetical_liquid_side_kla_per_s"],
            "aqueous_hcn_plus_cn_mol": item["final_aqueous_HCN_plus_CN_minus_mol"],
            "gas_hcn_mol": item["final_gas_HCN_mol"],
            "numerical_reporting_floor_mol": item["HCN_numerical_reporting_floor_mol"],
            "aqueous_amount_above_reporting_floor": (
                item["final_aqueous_HCN_plus_CN_minus_mol"]
                > item["HCN_numerical_reporting_floor_mol"]
            ),
        }
        for item in reports
    ]
    return {
        "status": (
            "feed_only_transfer_sensitivity_controls_passed_not_yield_validation"
            if controls_pass else "feed_only_transfer_sensitivity_controls_failed"
        ),
        "duration_s": 3600,
        "integration_steps": reference.get("total_integration_steps"),
        "same_feed_duration_energy_and_resolution": same_conditions,
        "contact_labels_valid": contact_labels_valid,
        "material_and_charged_controls_pass": controls_pass,
        "scenarios": scenarios,
        "failed_scenarios": failures,
        "experimental_hcn_amount_loaded": False,
        "science_boundary": (
            "kLa values and spark exposure are hypothetical. This is a numerical "
            "sensitivity map, not a measured kLa interval, calibrated prediction, "
            "full-duration comparison, or HCN yield validation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cycles", type=int, choices=(1, 2, 4), default=1)
    parser.add_argument("--substeps", type=int, choices=(48,), default=48)
    arguments = parser.parse_args()
    reports: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    for value in (*DEFAULT_KLA_PER_S, None):
        label = "step_equilibrium" if value is None else f"kLa={value:g} 1/s"
        print(f"Running feed-only contact scenario {label}", file=sys.stderr, flush=True)
        try:
            reports.append(run_probe(
                duration_s=3600, cycles=arguments.cycles, substeps=arguments.substeps,
                gas_water_kla_per_s=value,
            ))
        except ValueError as error:
            failures.append({
                "gas_water_contact_model": (
                    "step_equilibrium" if value is None else "finite_rate_contact"
                ),
                "hypothetical_liquid_side_kla_per_s": value,
                "error": str(error),
                "numerical_result_accepted": False,
            })
    print(json.dumps(summarize_transfer_sensitivity(reports, failures), indent=2))


if __name__ == "__main__":
    main()
