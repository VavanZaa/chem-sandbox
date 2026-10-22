"""Assay-free dry-flow field sensitivity; no HCN observation enters the solver."""
from __future__ import annotations

import json
from typing import Any

from scripts.run_dry_flow_hcn_probe import run_dry_flow_probe


FIELDS_TD = (189.0, 200.0, 300.0, 435.0)


def audit_dry_field_sweep() -> dict[str, Any]:
    rows = []
    for field in FIELDS_TD:
        prediction = run_dry_flow_probe(
            source_forcing="evolving_eedf",
            substeps=16,
            reduced_electric_field_Td=field,
        )
        outlet = prediction["outlet_packet_diagnostics"]
        yields = outlet["product_carbon_yield_fraction_by_species"]
        rows.append({
            "reduced_electric_field_Td": field,
            "methane_conversion_fraction": outlet["methane_conversion_fraction"],
            "c2h6_carbon_yield_fraction": yields["C2H6"],
            "hcn_carbon_yield_fraction_below_reporting_floor": yields["HCN"],
            "hcn_amount_below_1e-12_mol_reporting_floor": prediction["outlet_hcn_mol"] < 1e-12,
            "unresolved_n2_total_outcome_opportunity_mol": prediction["plasma"][
                "unresolved_n2_total_outcome_opportunity_mol"
            ],
            "conditional_sub14_526eV_ground_n_opportunity_mol": prediction["plasma"][
                "conditional_sub14_526eV_ground_n_opportunity_mol"
            ],
            "n2_product_state_source_excluded": (
                prediction["heavy_chemistry"]["excluded_fragment_channels_by_id"].get(
                    "N2_dissociation_outcome_000", {}
                ).get("reason") == "unresolved_electronic_product_state"
            ),
            "material_ledger_closed": prediction["material_ledger"]["closed"],
            "collision_energy_upper_bound_closed": prediction["energy_ledger"]["closed_upper_bound"],
        })
    return {
        "status": "dry_field_sensitivity_not_experimental_validation",
        "field_provenance": {
            "189_and_435_Td": "paper average-field estimates, not resolved local streamer fields",
            "200_Td": "project hypothetical scenario",
            "300_Td": "paper model-selected setting fitted to product trends; not independent",
        },
        "frozen_scenario": "1% CH4, 99% N2, 555.15 K uniform gas, 298.15 K flow reference, full-packet exposure",
        "rows": rows,
        "measured_hcn_assay_loaded": False,
        "interpretation": (
            "This is a sensitivity of an uncalibrated whole-packet scenario. "
            "The unresolved N2 collision opportunity is not a propagated N-atom source, "
            "and sub-floor HCN amounts are not validated yields."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_field_sweep(), indent=2))
