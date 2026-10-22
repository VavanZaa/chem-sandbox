"""Assay-free uniform-temperature sensitivity for the dry-flow packet.

Only 555.15 K is the steady downstream-probe reading in the paper's separate
28.6 W thermometry run. None of these uniform values is a measured reactor
temperature profile or a calibrated operating condition.
"""
from __future__ import annotations

import json
from typing import Any

from scripts.run_dry_flow_hcn_probe import run_dry_flow_probe


TEMPERATURES_K = (450.0, 555.15, 700.0)


def audit_dry_temperature_sweep() -> dict[str, Any]:
    rows = []
    for temperature in TEMPERATURES_K:
        prediction = run_dry_flow_probe(
            source_forcing="evolving_eedf",
            substeps=16,
            uniform_gas_temperature_K=temperature,
        )
        conditions = prediction["conditions"]
        outlet = prediction["outlet_packet_diagnostics"]
        yields = outlet["product_carbon_yield_fraction_by_species"]
        rows.append({
            "uniform_gas_temperature_K": temperature,
            "equals_downstream_probe_value": conditions[
                "uniform_temperature_equals_downstream_probe_value"
            ],
            "ideal_gas_packet_residence_s": conditions["simulated_hot_packet_residence_s"],
            "ch4_conversion_fraction": outlet["methane_conversion_fraction"],
            "c2h6_carbon_yield_fraction": yields["C2H6"],
            "hcn_carbon_yield_fraction_below_reporting_floor": yields["HCN"],
            "hcn_amount_below_1e-12_mol_reporting_floor": prediction["outlet_hcn_mol"] < 1e-12,
            "n2_product_state_source_excluded": (
                prediction["heavy_chemistry"]["excluded_fragment_channels_by_id"].get(
                    "N2_dissociation_outcome_000", {}
                ).get("reason") == "unresolved_electronic_product_state"
            ),
            "material_ledger_closed": prediction["material_ledger"]["closed"],
            "collision_energy_upper_bound_closed": prediction["energy_ledger"]["closed_upper_bound"],
        })
    return {
        "status": "dry_uniform_temperature_sensitivity_not_profile_reconstruction",
        "source_temperature_basis": (
            "paper measured approximately 282 C 1 cm above the reactor after "
            "10 min steady operation in a separate 28.6 W thermometry run"
        ),
        "nominal_residence_s": 0.47,
        "assumed_flow_reference_temperature_K": 298.15,
        "rows": rows,
        "measured_hcn_assay_loaded": False,
        "interpretation": (
            "All rows assume a uniform reactor temperature. Only 555.15 K "
            "equals a downstream-probe reading; 450 and 700 K are sensitivity "
            "cases. The actual axial and radial reactor temperature profiles "
            "and associated residence-time distribution remain unknown."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_temperature_sweep(), indent=2))
