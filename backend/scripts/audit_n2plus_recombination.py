"""Compare private N2+ DR support with a feed-only self-consistent EEDF.

No HCN assay or expected product is read. Output contains only summarized
coefficients and coverage, never the contributor cross-section table.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from app.prebiotic.boltzmann import IST_CO2_O2_DATASET_ID, solve_self_consistent_eedf
from app.prebiotic.n2plus_recombination import audit_n2plus_recombination
from app.prebiotic.schema import BoltzmannEEDFRequest


def main() -> None:
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--energy-cells", type=int, default=800)
    options = argument_parser.parse_args()
    feed_path = (
        Path(__file__).resolve().parents[1]
        / "data" / "benchmarks" / "barth_2024_feed_only_220725001_v1.json"
    )
    feed = json.loads(feed_path.read_text(encoding="utf-8"))
    fractions = dict(feed["measured_initial_gas_mole_fractions"])
    fractions["N2"] = 1.0 - sum(fractions.values())
    assumptions = feed["diagnostic_only_assumptions"]
    request = BoltzmannEEDFRequest(
        collision_set=IST_CO2_O2_DATASET_ID,
        mole_fractions=fractions,
        gas_temperature_K=assumptions["gas_and_water_temperature_K"],
        reduced_electric_field_Td=assumptions["reduced_electric_field_Td"],
        energy_cells=options.energy_cells,
        acknowledge_two_term_limit=True,
    )
    solved = solve_self_consistent_eedf(
        request, include_same_set_cantera_reference=False
    )
    audit = audit_n2plus_recombination(solved["eedf"])
    print(json.dumps({
        "feed_id": feed["experiment_id"],
        "feed_boundary": feed["diagnostic_boundary"],
        "field_Td_is_assumed_not_measured": True,
        "energy_cells": options.energy_cells,
        "eedf_mean_energy_eV": solved["eedf"]["mean_energy_eV"],
        "eedf_converged": solved["solver"]["converged"],
        "recombination_audit": audit,
    }, indent=2))


if __name__ == "__main__":
    main()
