"""Reproduce the paper's independent N2 optical field-ratio calculation.

This is a read-only spectroscopy diagnostic. No HCN, CH4 conversion, or
product selectivity is used to choose an electric field or alter the solver.
"""
from __future__ import annotations

import json
from typing import Any

from scipy.constants import Boltzmann

from app.prebiotic.boltzmann import (
    COMMUNITY_CH4_DATASET_ID,
    PHELPS_N2_SHA256,
    solve_self_consistent_eedf,
)
from app.prebiotic.schema import BoltzmannEEDFRequest


PAPER_SUPPLEMENT_EQUATION = "S24"
PAPER_OPTICAL_RATIO = 3.47e-3  # Supplementary Figure S12A, separate 28.6 W optical run.
GAS_TEMPERATURE_K = 550.0  # The paper's spectroscopic rate calculation.
GAS_PRESSURE_PA = 101_325.0
FEED_MOLE_FRACTIONS = {"CH4": 0.01, "N2": 0.99}
RADIATIVE_RATE_N2_PLUS_B_S_INV = 1.14e7
RADIATIVE_RATE_N2_C_S_INV = 1.31e7
N2_QUENCH_RATE_N2_PLUS_B_CM3_S = 8.84e-10
N2_QUENCH_RATE_N2_C_CM3_S = 1.0e-11
FIELD_GRID_TD = (100.0, 189.0, 300.0, 435.0, 500.0)
SUPPLEMENT_PDF_SHA256 = "43520da71dc844752958cb586cd251e9b7f8e49b707ff53a1c91e33f3dddc084"


def audit_dry_optical_field() -> dict[str, Any]:
    nitrogen_number_density_cm3 = (
        GAS_PRESSURE_PA / (Boltzmann * GAS_TEMPERATURE_K) / 1e6
        * FEED_MOLE_FRACTIONS["N2"]
    )
    n2_plus_b_quench_factor = RADIATIVE_RATE_N2_PLUS_B_S_INV / (
        RADIATIVE_RATE_N2_PLUS_B_S_INV
        + N2_QUENCH_RATE_N2_PLUS_B_CM3_S * nitrogen_number_density_cm3
    )
    n2_c_quench_factor = RADIATIVE_RATE_N2_C_S_INV / (
        RADIATIVE_RATE_N2_C_S_INV
        + N2_QUENCH_RATE_N2_C_CM3_S * nitrogen_number_density_cm3
    )
    def calculate(
        field_td: float,
        energy_max_eV: float,
        energy_cells: int,
        mixture: dict[str, float] | None = None,
        gas_temperature_K: float = GAS_TEMPERATURE_K,
    ) -> dict[str, Any]:
        mixture = FEED_MOLE_FRACTIONS if mixture is None else mixture
        local_n2_density_cm3 = (
            GAS_PRESSURE_PA / (Boltzmann * gas_temperature_K) / 1e6
            * mixture["N2"]
        )
        local_n2_plus_b_factor = RADIATIVE_RATE_N2_PLUS_B_S_INV / (
            RADIATIVE_RATE_N2_PLUS_B_S_INV
            + N2_QUENCH_RATE_N2_PLUS_B_CM3_S * local_n2_density_cm3
        )
        local_n2_c_factor = RADIATIVE_RATE_N2_C_S_INV / (
            RADIATIVE_RATE_N2_C_S_INV
            + N2_QUENCH_RATE_N2_C_CM3_S * local_n2_density_cm3
        )
        eedf = solve_self_consistent_eedf(BoltzmannEEDFRequest(
            collision_set=COMMUNITY_CH4_DATASET_ID,
            mole_fractions=mixture,
            gas_temperature_K=gas_temperature_K,
            reduced_electric_field_Td=field_td,
            energy_max_eV=energy_max_eV,
            energy_cells=energy_cells,
            acknowledge_two_term_limit=True,
        ), include_same_set_cantera_reference=(
            mixture == FEED_MOLE_FRACTIONS
            and gas_temperature_K == GAS_TEMPERATURE_K
            and field_td == 435.0 and energy_max_eV == 150.0 and energy_cells == 1000
        ))
        rates = {
            rate["product_label"]: rate["rate_coefficient_m3_s"]
            for rate in eedf["rate_coefficients"]
            if rate["target"] == "N2" and rate["product_label"] in {
                "N2^+(B2SIGMA)", "N2(C3)"
            }
        }
        if set(rates) != {"N2^+(B2SIGMA)", "N2(C3)"}:
            raise ValueError("N2+(B) and N2(C) collision channels are both required")
        predicted_ratio = (
            rates["N2^+(B2SIGMA)"] / rates["N2(C3)"]
            * local_n2_plus_b_factor / local_n2_c_factor
        )
        comparison = eedf["same_collision_set_reference"]
        cantera_optical_ratio = None
        cantera_rate_factors = None
        if comparison.get("available_for_this_mixture"):
            optical_channels = {
                row["product_label"]: row
                for row in comparison["rate_comparison"]["selected_n2_optical_channels"]
            }
            if set(optical_channels) != {"N2^+(B2SIGMA)", "N2(C3)"}:
                raise ValueError("Cantera lacks one of the two N2 optical channels")
            cantera_optical_ratio = (
                optical_channels["N2^+(B2SIGMA)"]["cantera_eedf_rate_coefficient_m3_s"]
                / optical_channels["N2(C3)"]["cantera_eedf_rate_coefficient_m3_s"]
                * local_n2_plus_b_factor / local_n2_c_factor
            )
            cantera_rate_factors = {
                label: channel["factor_difference"]
                for label, channel in optical_channels.items()
            }
        return {
            "reduced_electric_field_Td": field_td,
            "energy_max_eV": energy_max_eV,
            "energy_cells": energy_cells,
            "gas_temperature_K": gas_temperature_K,
            "eedf_mole_fractions": mixture,
            "n2_plus_b_direct_rate_m3_s": rates["N2^+(B2SIGMA)"],
            "n2_c_direct_rate_m3_s": rates["N2(C3)"],
            "predicted_spectroscopic_intensity_ratio": predicted_ratio,
            "ratio_to_paper_measured_intensity_ratio": predicted_ratio / PAPER_OPTICAL_RATIO,
            "same_cross_section_cantera_predicted_intensity_ratio": cantera_optical_ratio,
            "same_cross_section_cantera_relative_mean_energy_difference": (
                comparison.get("relative_mean_energy_difference")
            ),
            "same_cross_section_cantera_rate_factor_differences": cantera_rate_factors,
            "same_cross_section_cantera_reference_available": comparison.get(
                "available_for_this_mixture", False
            ),
            "eedf_converged": eedf["solver"]["converged"],
            "eedf_numerical_warnings": eedf["numerical_warnings"],
        }
    rows = [calculate(field_td, 150.0, 1000) for field_td in FIELD_GRID_TD]
    grid_at_435 = [
        calculate(435.0, maximum, cells)
        for maximum, cells in ((60.0, 400), (100.0, 700), (150.0, 1000))
    ]
    return {
        "status": "independent_optical_field_cross_section_audit_not_field_calibration",
        "paper_supplement_equation": PAPER_SUPPLEMENT_EQUATION,
        "supplement_pdf_sha256": SUPPLEMENT_PDF_SHA256,
        "paper_measured_intensity_ratio": PAPER_OPTICAL_RATIO,
        "paper_optical_run_not_identical_to_figure1_nominal_power": True,
        "mixture": FEED_MOLE_FRACTIONS,
        "gas_temperature_K": GAS_TEMPERATURE_K,
        "nitrogen_number_density_cm3": nitrogen_number_density_cm3,
        "n2_plus_b_radiative_quench_factor": n2_plus_b_quench_factor,
        "n2_c_radiative_quench_factor": n2_c_quench_factor,
        "collision_set": COMMUNITY_CH4_DATASET_ID,
        "installed_cantera_phelps_sha256": PHELPS_N2_SHA256,
        "installed_cantera_phelps_export_date": "2025-08-01",
        "paper_phelps_lxcat_retrieval_date": "2023-03-01",
        "field_grid": rows,
        "energy_grid_sensitivity_at_435_Td": grid_at_435,
        "pure_n2_eedf_sensitivity_at_435_Td": calculate(
            435.0, 150.0, 1000, {"N2": 1.0}
        ),
        "room_temperature_sensitivity_at_435_Td": calculate(
            435.0, 150.0, 1000, gas_temperature_K=298.15
        ),
        "hcn_or_product_assay_used": False,
        "electric_field_selected_by_this_audit": False,
        "photodetector_response_and_full_emission_transport_validated": False,
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_optical_field(), indent=2))
