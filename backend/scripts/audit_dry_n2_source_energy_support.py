"""Assay-free N2 fragment-state evidence window for the dry CH4/N2 feed.

The collision channel is a total dissociation outcome. Its EEDF-folded
energy support is not a product-state branch and is never injected as N atoms.
"""
from __future__ import annotations

import json
from typing import Any

from app.prebiotic.boltzmann import COMMUNITY_CH4_DATASET_ID, solve_self_consistent_eedf
from app.prebiotic.schema import BoltzmannEEDFRequest
from scripts.run_dry_flow_hcn_probe import DOWNSTREAM_PROBE_TEMPERATURE_K


def audit_dry_n2_source_energy_support() -> dict[str, Any]:
    mixture = {"CH4": 0.01, "N2": 0.99}

    def calculate(field_Td: float, energy_max_eV: float, energy_cells: int) -> dict[str, Any]:
        result = solve_self_consistent_eedf(BoltzmannEEDFRequest(
            collision_set=COMMUNITY_CH4_DATASET_ID,
            mole_fractions=mixture,
            gas_temperature_K=DOWNSTREAM_PROBE_TEMPERATURE_K,
            reduced_electric_field_Td=field_Td,
            energy_max_eV=energy_max_eV,
            energy_cells=energy_cells,
            acknowledge_two_term_limit=True,
        ))
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
            "product_state_branching_resolved": False,
        }

    return {
        "status": "dry_total_n2_energy_support_not_product_branching",
        "feed_mole_fractions": mixture,
        "uniform_gas_temperature_K_scenario": DOWNSTREAM_PROBE_TEMPERATURE_K,
        "collision_set": COMMUNITY_CH4_DATASET_ID,
        "field_sweep": [
            calculate(field, 150.0, 1000)
            for field in (189.0, 200.0, 300.0, 435.0)
        ],
        "field_provenance": {
            "189_and_435_Td": "paper-reported average experimental estimates; not resolved streamer fields",
            "200_Td": "project hypothetical scenario",
            "300_Td": "paper model-selected condition; not independent validation input",
        },
        "grid_sweep_at_200_Td": [
            calculate(200.0, maximum, cells)
            for maximum, cells in ((60.0, 400), (100.0, 700), (150.0, 1000))
        ],
        "source_state_branch_inferred": False,
        "atom_source_injected": False,
        "hcn_assay_read": False,
        "hcn_yield_predicted": False,
        "interpretation": (
            "These fractions partition the EEDF-folded total N2 outcome across "
            "electron-energy windows, including Cosby's 18.5-148.5 eV "
            "fragment-state measurement. They are not N(4S)/N(2D) yields, "
            "a measured streamer exposure, or an HCN prediction."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_n2_source_energy_support(), indent=2))
