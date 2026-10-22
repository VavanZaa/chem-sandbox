"""Feed-only, read-only N2 product-state evidence-window sensitivity audit.

Uses the Barth initial gas composition but never reads an HCN assay. The
Tabata cross section is a *total* neutral-dissociation outcome, not a
state-resolved source. Energy-window fractions must not become atom-state
branching probabilities or HCN yields.
"""

from __future__ import annotations

import json
from typing import Any

from app.prebiotic.boltzmann import (
    IST_CO2_O2_DATASET_ID,
    solve_self_consistent_eedf,
)
from app.prebiotic.schema import BoltzmannEEDFRequest
from scripts.probe_barth_feed_pathway import FEED_CARD


def audit_barth_n2_source_energy_support() -> dict[str, Any]:
    feed = json.loads(FEED_CARD.read_text(encoding="utf-8"))
    measured = {
        species: float(value)
        for species, value in feed["measured_initial_gas_mole_fractions"].items()
    }
    mixture = {**measured, "N2": 1.0 - sum(measured.values())}
    if mixture["N2"] <= 0:
        raise ValueError("Barth initial feed leaves no N2 balance")
    temperature_K = float(feed["diagnostic_only_assumptions"]["gas_and_water_temperature_K"])
    declared_field_Td = float(feed["diagnostic_only_assumptions"]["reduced_electric_field_Td"])

    def calculate(field_Td: float, energy_max_eV: float, energy_cells: int) -> dict[str, Any]:
        request = BoltzmannEEDFRequest(
            collision_set=IST_CO2_O2_DATASET_ID,
            mole_fractions=mixture,
            gas_temperature_K=temperature_K,
            reduced_electric_field_Td=field_Td,
            energy_max_eV=energy_max_eV,
            energy_cells=energy_cells,
            acknowledge_two_term_limit=True,
        )
        result = solve_self_consistent_eedf(request)
        channel = next(
            item for item in result["rate_coefficients"]
            if item["target"] == "N2" and item.get("product_state_unresolved")
        )
        return {
            "reduced_electric_field_Td": field_Td,
            "energy_max_eV": energy_max_eV,
            "energy_cells": energy_cells,
            "eedf_converged": result["solver"]["converged"],
            "total_n2_neutral_dissociation_rate_coefficient_m3_s": channel[
                "rate_coefficient_m3_s"
            ],
            "product_state_energy_support": channel["product_state_energy_support"],
            "ground_atom_energetic_window": channel["ground_atom_energetic_window"],
            "tabata_analytic_fit_diagnostic": channel["analytic_fit_diagnostic"],
            "product_state_branching_resolved": False,
        }

    field_sweep = [calculate(field, 60.0, 400) for field in (50.0, 100.0, 150.0, 200.0)]
    grid_sweep = [
        calculate(declared_field_Td, maximum, cells)
        for maximum, cells in ((60.0, 400), (100.0, 700), (150.0, 1000))
    ]
    return {
        "status": "total_n2_source_energy_support_audit_not_product_branching",
        "feed_card": feed["id"],
        "feed_card_contains_hcn_assay": False,
        "collision_set": IST_CO2_O2_DATASET_ID,
        "temperature_K": temperature_K,
        "field_sweep": field_sweep,
        "grid_sweep_at_declared_hypothetical_field": grid_sweep,
        "source_state_branch_inferred": False,
        "conditional_ground_atom_energetic_window_used_as_source": False,
        "hcn_yield_predicted": False,
        "interpretation": (
            "Fractions partition the EEDF-folded total N2 event rate by electron "
            "energy relative to Cosby 1993's 18.5-148.5 eV fragment-state "
            "measurement and the nominal 14.526 eV first no-ground-N product "
            "threshold. The electric field is hypothetical, and neither "
            "these fractions nor the total outcome rate are a measured state "
            "branch, injected atom source, or HCN yield."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_barth_n2_source_energy_support(), indent=2))
