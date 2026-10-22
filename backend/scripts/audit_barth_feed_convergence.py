"""Read-only substep screen for the feed-only Barth pathway diagnostic.

The output assesses numerical resolution, not apparatus fidelity or HCN yield.
No experimental HCN assay is loaded by this module.
"""

from __future__ import annotations

import json
from typing import Any

from scripts.probe_barth_feed_pathway import run_probe


CONFIGURATIONS = ((1, 48), (2, 48), (4, 48))


def summarize_convergence(reports: list[dict[str, Any]]) -> dict[str, Any]:
    if len(reports) != 3:
        raise ValueError("The convergence screen requires three ordered resolutions")
    expected_steps = [48, 96, 192]
    if [item["total_integration_steps"] for item in reports] != expected_steps:
        raise ValueError("Convergence reports must be ordered at 48, 96, 192 steps")
    same_conditions = (
        [(item.get("cycles"), item.get("integration_substeps_per_cycle")) for item in reports]
        == list(CONFIGURATIONS)
        and all(item.get("duration_s") == 3600 for item in reports)
        and all(item.get("spark_on") is True for item in reports)
        and all(item.get("energy_basis") == "declared_electrical_estimate" for item in reports)
        and all(item.get("initial_feed_mol") == reports[0].get("initial_feed_mol") for item in reports)
        and all(item.get("gas_water_contact_model") == reports[0].get("gas_water_contact_model") for item in reports)
        and all(item.get("gas_water_contact_solver") == reports[0].get("gas_water_contact_solver") for item in reports)
        and all(item.get("primary_source_forcing") == reports[0].get("primary_source_forcing") for item in reports)
        and all(item.get("formaldehyde_hydration_model") == reports[0].get("formaldehyde_hydration_model") for item in reports)
        and all(item.get("hypothetical_liquid_side_kla_per_s") == reports[0].get("hypothetical_liquid_side_kla_per_s") for item in reports)
        and "HCN" not in reports[0].get("initial_feed_mol", {})
    )
    controls_pass = all(
        item["status"] == "feed_matched_short_pathway_diagnostic_not_yield_validation"
        and item["material_ledger_closed"] is True
        and item["charged_branch_material_ledger_closed"] is True
        and item["charged_branch_failed_substeps"] == 0
        and item["charged_branch_unprojected_ion_roundoff_clip_mol"] <= 1e-15
        for item in reports
    ) and same_conditions
    endpoints = (
        "final_aqueous_HCN_plus_CN_minus_mol",
        "charged_final_aqueous_HCN_plus_CN_minus_mol",
    )
    comparisons = {}
    for endpoint in endpoints:
        coarse, medium, fine = (float(item[endpoint]) for item in reports)
        comparisons[endpoint] = {
            "48_steps_mol": coarse,
            "96_steps_mol": medium,
            "192_steps_mol": fine,
            "96_to_192_relative_change": (
                abs(fine - medium) / abs(fine) if fine > 0 else None
            ),
        }
    resolution_screen_pass = controls_pass and all(
        item["96_to_192_relative_change"] is not None
        and item["96_to_192_relative_change"] <= 0.05
        for item in comparisons.values()
    )
    return {
        "status": (
            "one_hour_feed_only_numerical_screen_passed_not_yield_validation"
            if resolution_screen_pass else "one_hour_feed_only_numerical_screen_failed"
        ),
        "criterion": "both primary and charged aqueous HCN+CN- change by <=5% from 96 to 192 steps, with material and charged controls passing",
        "comparison": comparisons,
        "material_and_charged_controls_pass": controls_pass,
        "same_feed_duration_and_electrical_basis": same_conditions,
        "primary_maximum_projection_total_absolute_mol": max(
            float(item["primary_negative_species_projection_total_absolute_mol"])
            for item in reports
        ),
        "primary_maximum_positivity_retry_substeps": max(
            int(item.get("primary_positivity_retry_substeps") or 0)
            for item in reports
        ),
        "primary_maximum_HCN_projection_adjustment_total_absolute_mol": max(
            float(item.get("primary_HCN_projection_adjustment_total_absolute_mol") or 0.0)
            for item in reports
        ),
        "charged_maximum_unprojected_roundoff_clip_mol": max(
            float(item["charged_branch_unprojected_ion_roundoff_clip_mol"])
            for item in reports
        ),
        "formaldehyde_hydration_model": reports[0].get("formaldehyde_hydration_model"),
        "primary_source_forcing": reports[0].get("primary_source_forcing"),
        "formaldehyde_hydration_pH_extrapolation_used": any(
            bool(item.get("formaldehyde_hydration_pH_extrapolation_used")) for item in reports
        ),
        "minimum_ideal_pH": min(
            (float(item["formaldehyde_hydration_minimum_ideal_pH"])
             for item in reports
             if item.get("formaldehyde_hydration_minimum_ideal_pH") is not None),
            default=None,
        ),
        "science_boundary": (
            "This screens one hour with hypothetical spark exposure and no HCN assay. "
            "A diagnostic formaldehyde mode may explicitly extrapolate below measured pH. "
            "It does not validate the full 925-minute experiment, the mechanism, "
            "deposited energy, or product amount."
        ),
    }


def main() -> None:
    reports = [
        run_probe(duration_s=3600, cycles=cycles, substeps=substeps)
        for cycles, substeps in CONFIGURATIONS
    ]
    print(json.dumps(summarize_convergence(reports), indent=2))


if __name__ == "__main__":
    main()
