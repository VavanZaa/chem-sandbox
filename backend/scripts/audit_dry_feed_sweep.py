"""Freeze an assay-free dry-flow feed sweep for external trend evaluation.

The chosen feed points span the paper's Figure 1A range but no measured
product amount, selectivity, or fitted streamer exposure enters this solver.
"""
from __future__ import annotations

import json
import math
from typing import Any

from scripts.run_dry_flow_hcn_probe import run_dry_flow_probe


METHANE_FEED_FRACTIONS = (0.001, 0.01, 0.07)
FIGURE_1A_METHANE_FEED_FRACTIONS = (0.001, 0.0025, 0.005, 0.01, 0.02, 0.035, 0.07)


def audit_dry_feed_sweep(
    *, reduced_electric_field_Td: float = 200.0,
    trace_methane_neutral_pathways: bool = False,
    methane_feed_fractions: tuple[float, ...] = METHANE_FEED_FRACTIONS,
    substeps: int = 16,
) -> dict[str, Any]:
    if (
        not methane_feed_fractions
        or any(not math.isfinite(value) or not 0 < value < 1 for value in methane_feed_fractions)
        or any(a >= b for a, b in zip(methane_feed_fractions, methane_feed_fractions[1:]))
    ):
        raise ValueError("Methane feed fractions must be strictly increasing and between zero and one")
    if not isinstance(substeps, int) or not 1 <= substeps <= 128:
        raise ValueError("Dry-flow substeps must be an integer in [1, 128]")
    rows = []
    for methane_fraction in methane_feed_fractions:
        prediction = run_dry_flow_probe(
            methane_fraction=methane_fraction,
            source_forcing="evolving_eedf",
            substeps=substeps,
            reduced_electric_field_Td=reduced_electric_field_Td,
            trace_methane_neutral_pathways=trace_methane_neutral_pathways,
        )
        outlet = prediction["outlet_packet_diagnostics"]
        n2_opportunities = prediction["plasma"][
            "n2_excitation_opportunity_mol_by_product_label"
        ]
        n2_a3_opportunities_mol = sum(
            amount for label, amount in n2_opportunities.items()
            if label.startswith("N2(A3,")
        )
        n2_aprime_opportunities_mol = n2_opportunities.get("N2(a'1)", 0.0)
        row = {
            "feed_CH4_mole_fraction": methane_fraction,
            "inlet_CH4_mol": prediction["inlet_total_mol"] * methane_fraction,
            "ch4_conversion_fraction": outlet["methane_conversion_fraction"],
            "ch4_electron_fragment_source_fraction_of_feed": outlet[
                "methane_electron_fragment_source_fraction_of_feed"
            ],
            "ch4_net_neutral_chemistry_consumption_fraction_of_feed": outlet[
                "methane_net_neutral_chemistry_consumption_fraction_of_feed"
            ],
            "effective_electron_density_m3_range": prediction["plasma"][
                "effective_electron_density_m3_range"
            ],
            "electron_density_was_power_capped": prediction["plasma"][
                "density_was_power_capped"
            ],
            "declared_discharge_J": prediction["energy_ledger"]["declared_discharge_J"],
            "resolved_fragment_J": prediction["energy_ledger"]["resolved_fragment_J"],
            "uncoupled_n2_a3_excitation_opportunities_per_feed_CH4": (
                n2_a3_opportunities_mol / (prediction["inlet_total_mol"] * methane_fraction)
            ),
            "uncoupled_n2_aprime_excitation_opportunities_per_feed_CH4": (
                n2_aprime_opportunities_mol / (prediction["inlet_total_mol"] * methane_fraction)
            ),
            "c2h6_carbon_yield_fraction": outlet[
                "product_carbon_yield_fraction_by_species"
            ]["C2H6"],
            "c2h6_carbon_selectivity_fraction": outlet[
                "product_carbon_selectivity_fraction_by_species"
            ]["C2H6"],
            "h2_hydrogen_selectivity_fraction": outlet[
                "product_hydrogen_selectivity_fraction_by_species"
            ]["H2"],
            "hcn_carbon_yield_fraction_below_reporting_floor": outlet[
                "product_carbon_yield_fraction_by_species"
            ]["HCN"],
            "hcn_amount_below_1e-12_mol_reporting_floor": (
                prediction["outlet_hcn_mol"] < 1e-12
            ),
            "n2_product_state_source_excluded": (
                prediction["heavy_chemistry"]["excluded_fragment_channels_by_id"].get(
                    "N2_dissociation_outcome_000", {}
                ).get("reason") == "unresolved_electronic_product_state"
            ),
            "material_ledger_closed": prediction["material_ledger"]["closed"],
            "collision_energy_upper_bound_closed": prediction["energy_ledger"][
                "closed_upper_bound"
            ],
        }
        if trace_methane_neutral_pathways:
            row["methane_neutral_pathways"] = prediction["heavy_chemistry"][
                "methane_neutral_pathways"
            ]
        rows.append(row)
    return {
        "status": "dry_feed_fraction_sensitivity_not_experimental_validation",
        "scenario": (
            "555.15 K uniform gas, 298.15 K assumed flow reference, "
            f"{reduced_electric_field_Td:g} Td, full-packet power-capped exposure"
        ),
        "reduced_electric_field_Td_scenario": reduced_electric_field_Td,
        "substeps": substeps,
        "trace_methane_neutral_pathways": trace_methane_neutral_pathways,
        "field_provenance": (
            "189 Td is the paper's approximate average-field estimate across "
            "feed fractions, not a local streamer field; 200 Td is a project "
            "hypothetical reference. Neither was selected by product matching."
        ),
        "rows": rows,
        "measured_product_assay_loaded": False,
        "n2_excitation_opportunities_are_product_yields": False,
        "selectivity_basis_matches_paper_exactly": False,
        "selectivity_basis_caveat": (
            "Modeled non-CH4 species include radicals and omit deposited coke; "
            "the paper uses measured stable products and coke. Carbon and "
            "hydrogen selectivity denominators follow converted methane atom "
            "counts, but compare only trends, not selectivity magnitudes."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_feed_sweep(), indent=2))
