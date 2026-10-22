"""Assay-free dry CH4/N2 residence sensitivity for held-out Figure 1D comparison."""
from __future__ import annotations

import json
from typing import Any

from scripts.run_dry_flow_hcn_probe import run_dry_flow_probe


NOMINAL_RESIDENCE_TIMES_S = (1.88, 0.94, 0.63, 0.47)


def audit_dry_residence_sweep(*, reduced_electric_field_Td: float = 184.0) -> dict[str, Any]:
    rows = []
    for nominal_residence_s in NOMINAL_RESIDENCE_TIMES_S:
        prediction = run_dry_flow_probe(
            methane_fraction=0.01,
            nominal_residence_s=nominal_residence_s,
            source_forcing="evolving_eedf",
            substeps=16,
            reduced_electric_field_Td=reduced_electric_field_Td,
        )
        outlet = prediction["outlet_packet_diagnostics"]
        rows.append({
            "nominal_residence_s": nominal_residence_s,
            "nominal_inferred_flow_sccm": prediction["flow_temperature_basis_diagnostic"][
                "nominal_volume_flow_implied_by_reported_residence_sccm"
            ],
            "simulated_hot_packet_residence_s": prediction["conditions"][
                "simulated_hot_packet_residence_s"
            ],
            "declared_discharge_J_per_packet": prediction["energy_ledger"][
                "declared_discharge_J"
            ],
            "ch4_conversion_fraction": outlet["methane_conversion_fraction"],
            "c2h6_carbon_selectivity_fraction": outlet[
                "product_carbon_selectivity_fraction_by_species"
            ]["C2H6"],
            "hcn_mol_below_1e-12_reporting_floor": prediction["outlet_hcn_mol"] < 1e-12,
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
        "status": "dry_residence_sensitivity_not_experimental_validation",
        "scenario": (
            "1% CH4, 28.1 W, 555.15 K uniform gas, 298.15 K assumed flow "
            f"reference, {reduced_electric_field_Td:g} Td, full-packet power-capped exposure"
        ),
        "reduced_electric_field_Td_scenario": reduced_electric_field_Td,
        "field_provenance": (
            "Default 184 Td is the paper supplement Figure S7B average field "
            "estimated across its residence sweep, not a local streamer measurement."
        ),
        "rows": rows,
        "measured_product_assay_loaded": False,
        "residence_convention_matches_paper_exactly": False,
        "microdischarge_exposure_independently_constrained": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_residence_sweep(), indent=2))
