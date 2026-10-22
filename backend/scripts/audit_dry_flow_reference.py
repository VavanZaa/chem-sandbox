"""Assay-free GF40 flow-reference sensitivity for the dry packet.

The Brooks GF40 manual lists 0, 15, 20 C, and 70 F reference options. The
paper identifies GF40 controllers but not their installed option codes. The
existing 25 C project assumption is retained as a separate legacy scenario;
none of these values is asserted to be the experiment's actual setting.
"""
from __future__ import annotations

import json
from typing import Any

from scripts.run_dry_flow_hcn_probe import run_dry_flow_probe


GF40_LISTED_REFERENCE_TEMPERATURES_K = (273.15, 288.15, 293.15, 294.2611111111111)
PROJECT_LEGACY_ASSUMPTION_K = 298.15


def audit_dry_flow_reference() -> dict[str, Any]:
    rows = []
    for reference_temperature_K in (
        *GF40_LISTED_REFERENCE_TEMPERATURES_K,
        PROJECT_LEGACY_ASSUMPTION_K,
    ):
        prediction = run_dry_flow_probe(
            source_forcing="evolving_eedf",
            substeps=16,
            flow_reference_temperature_K=reference_temperature_K,
        )
        conditions = prediction["conditions"]
        flow = prediction["flow_temperature_basis_diagnostic"]
        outlet = prediction["outlet_packet_diagnostics"]
        rows.append({
            "reference_temperature_K": reference_temperature_K,
            "listed_GF40_option": (
                reference_temperature_K in GF40_LISTED_REFERENCE_TEMPERATURES_K
            ),
            "inferred_feed_molar_flow_mol_s": flow["inferred_feed_molar_flow_mol_s"],
            "ideal_gas_packet_residence_s": conditions["simulated_hot_packet_residence_s"],
            "declared_energy_J_per_packet": conditions["power_times_residence_J_per_gas_packet"],
            "ch4_conversion_fraction": outlet["methane_conversion_fraction"],
            "c2h6_carbon_selectivity_fraction": outlet[
                "product_carbon_selectivity_fraction_by_species"
            ]["C2H6"],
            "hcn_amount_below_1e-12_mol_reporting_floor": (
                prediction["outlet_hcn_mol"] < 1e-12
            ),
            "material_ledger_closed": prediction["material_ledger"]["closed"],
            "collision_energy_upper_bound_closed": prediction["energy_ledger"][
                "closed_upper_bound"
            ],
            "n2_product_state_source_excluded": (
                prediction["heavy_chemistry"]["excluded_fragment_channels_by_id"].get(
                    "N2_dissociation_outcome_000", {}
                ).get("reason") == "unresolved_electronic_product_state"
            ),
        })
    return {
        "status": "flow_reference_sensitivity_not_actual_gf40_setting_identification",
        "GF40_manual_url": (
            "https://cdn.brooksinstrument.com/-/media/brooks/documentation/"
            "products/mass-flow-controllers/elastomer-sealed/gf40/"
            "installation-manual-gf40.pdf"
        ),
        "GF40_manual_location": "2023 installation manual Table 5-1, reference temperature code XIII",
        "paper_reports_controller_model": "Brooks GF40",
        "paper_reports_installed_reference_code": False,
        "actual_experiment_reference_temperature_identified": False,
        "uniform_gas_scenario_temperature_K": 555.15,
        "rows": rows,
        "measured_hcn_assay_loaded": False,
        "interpretation": (
            "The GF40 manual lists possible reference options, not the code "
            "installed in the published experiment. These scenarios alter "
            "molar throughput and packet residence while retaining 200 sccm "
            "nominal V/Q; no option is selected by product matching."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_flow_reference(), indent=2))
