"""Assay-free electron-drift cost of total N2 dissociation opportunities.

This is a transport/source plausibility screen, not a nitrogen-atom yield:
the imported total neutral-dissociation cross section has unresolved atomic
electronic products, and local streamer exposure is not measured.
"""
from __future__ import annotations

import json
from typing import Any

from app.prebiotic.boltzmann import COMMUNITY_CH4_DATASET_ID, solve_self_consistent_eedf
from app.prebiotic.schema import BoltzmannEEDFRequest
from scripts.run_dry_flow_hcn_probe import DOWNSTREAM_PROBE_TEMPERATURE_K


def _intrinsic_dissociation_cost(
    *, methane_fraction: float, temperature_K: float, field_Td: float,
    energy_max_eV: float, energy_cells: int,
) -> dict[str, Any]:
    if not 0 <= methane_fraction < 1:
        raise ValueError("The N2 target fraction must be positive")
    n2_fraction = 1.0 - methane_fraction
    fractions = {"N2": n2_fraction}
    if methane_fraction:
        fractions["CH4"] = methane_fraction
    result = solve_self_consistent_eedf(BoltzmannEEDFRequest(
        collision_set=COMMUNITY_CH4_DATASET_ID,
        mole_fractions=fractions,
        gas_temperature_K=temperature_K,
        reduced_electric_field_Td=field_Td,
        energy_max_eV=energy_max_eV,
        energy_cells=energy_cells,
        acknowledge_two_term_limit=True,
    ))
    channel = next(
        item for item in result["rate_coefficients"]
        if item["target"] == "N2" and item.get("product_state_unresolved")
    )
    k_m3_s = float(channel["rate_coefficient_m3_s"])
    reduced_mobility = float(result["transport"]["reduced_electron_mobility_m_inv_V_inv_s_inv"])
    en_V_m2 = field_Td * 1e-21
    # Electron power density is e*n_e*mu*E^2 and event-rate density is
    # n_e*k*N_N2. Since mu*N_total is the reported reduced mobility,
    # the ratio in eV/event is (mu*N_total)*(E/N_total)^2/(k*x_N2).
    cost_eV_per_total_outcome = reduced_mobility * en_V_m2**2 / (k_m3_s * n2_fraction)
    return {
        "methane_mole_fraction": methane_fraction,
        "nitrogen_mole_fraction": n2_fraction,
        "gas_temperature_K": temperature_K,
        "reduced_electric_field_Td": field_Td,
        "energy_max_eV": energy_max_eV,
        "energy_cells": energy_cells,
        "eedf_converged": result["solver"]["converged"],
        "n2_total_neutral_outcome_rate_coefficient_m3_s": k_m3_s,
        "reduced_electron_mobility_m_inv_V_inv_s_inv": reduced_mobility,
        "electron_drift_eV_per_total_n2_outcome": cost_eV_per_total_outcome,
        "conditional_eV_per_atom_if_two_ground_atoms_per_outcome": cost_eV_per_total_outcome / 2,
        "conditional_eV_per_atom_if_one_ground_atom_per_outcome": cost_eV_per_total_outcome,
        "product_state_branch_known": False,
        "actual_ground_atom_energy_cost_determined": False,
    }


def audit_dry_nitrogen_atom_energy_cost() -> dict[str, Any]:
    pure_n2_grid = [
        _intrinsic_dissociation_cost(
            methane_fraction=0.0, temperature_K=300.0, field_Td=140.0,
            energy_max_eV=maximum, energy_cells=cells,
        )
        for maximum, cells in ((60.0, 400), (100.0, 700), (150.0, 1000))
    ]
    dry_scenario = _intrinsic_dissociation_cost(
        methane_fraction=0.01, temperature_K=DOWNSTREAM_PROBE_TEMPERATURE_K,
        field_Td=189.0, energy_max_eV=150.0, energy_cells=1000,
    )
    return {
        "status": "intrinsic_n2_total_outcome_energy_cost_not_atom_yield",
        "collision_set": COMMUNITY_CH4_DATASET_ID,
        "pure_n2_reference_condition_grid_check": pure_n2_grid,
        "dry_one_percent_methane_average_field_scenario": dry_scenario,
        "independent_comparison": {
            "doi": "10.1134/S1063780X13050085",
            "reported_ground_N4S_energy_cost_eV_per_atom": 337.0,
            "reported_uncertainty_eV_per_atom": 35.0,
            "conditions": "pure N2, 300 K, 1 atm, planar 2 mm Townsend DBD, approximately 140 Td",
            "measurement_basis": (
                "Popov inferred atom-production energy cost from an earlier "
                "TALIF ground-N density measurement, electrical power and "
                "modeled bulk/wall atom losses; it is not direct fragment branching."
            ),
            "like_for_like_with_dry_ch4_n2_reactor": False,
            "shares_cosby_total_cross_section_ancestry": True,
        },
        "density_or_residence_assumed": False,
        "nitrogen_atom_source_injected": False,
        "hcn_assay_loaded": False,
        "hcn_yield_predicted": False,
        "interpretation": (
            "Electron drift power divided by the EEDF-folded total neutral N2 "
            "outcome rate gives an intrinsic eV/outcome independent of electron "
            "density or residence. Dividing by two assumes two ground-state "
            "atoms per event and is an optimistic lower cost; dividing by one "
            "is only a one-ground-atom conditional case, not an upper bound. "
            "Actual ground-atom count, losses, wall chemistry and apparatus "
            "power partition remain unresolved, so the pure-N2 reference "
            "cannot validate the CH4/N2 HCN source or select a branch."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_nitrogen_atom_energy_cost(), indent=2))
