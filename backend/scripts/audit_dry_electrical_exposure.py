"""Assay-free consistency audit of the published dry-DBD electrical estimates.

The paper's supplement estimates E/N from emission and electron density from
current/area/drift speed. Neither is a resolved streamer measurement. This
audit checks its printed arithmetic and compares a uniform drift-power
requirement with measured average power. No occupancy is fed to chemistry.
"""
from __future__ import annotations

import json
import math
from typing import Any

from scipy.constants import R, elementary_charge

from app.prebiotic.boltzmann import COMMUNITY_CH4_DATASET_ID
from app.prebiotic.coupled_plasma import solve_mixture_plasma_snapshot
from app.prebiotic.schema import BoltzmannDischargeProgram
from scripts.run_dry_flow_hcn_probe import (
    AMBIENT_FLOW_REFERENCE_K,
    AVERAGE_DISCHARGE_POWER_W,
    PRESSURE_PA,
    REACTOR_VOLUME_L,
    RESIDENCE_S,
    _gri_inventory,
)


PAPER_ELECTRICAL_FIELD_TD = 435.0
PAPER_ELECTRON_DENSITY_M3 = 7.3e15  # 7.3e9 cm^-3, supplementary text.
PAPER_HALF_CYCLE_AVERAGE_CURRENT_A = 2.66e-3
PAPER_ANNULUS_AREA_CM2 = 0.314
PAPER_REPORTED_CURRENT_DENSITY_A_M2 = 407.0
PAPER_DRIFT_VELOCITY_M_S = 3.46e5  # 3.46e7 cm/s.
PAPER_EEDF_GAS_TEMPERATURE_K = 550.0
PAPER_S27_PRODUCT_FITTED_PLASMA_VOLUME_MM3 = 4.2
PAPER_S27_STATED_FILAMENT_DIAMETER_MM = 0.1
PAPER_S27_FILAMENT_LENGTH_MM = 2.0
PAPER_S27_REPORTED_MICRODISCHARGES_PER_CYCLE = 67


def audit_s27_microdischarge_count() -> dict[str, Any]:
    """Check the supplement's geometry only; never use it as solver exposure."""
    volume = PAPER_S27_PRODUCT_FITTED_PLASMA_VOLUME_MM3
    stated_diameter = PAPER_S27_STATED_FILAMENT_DIAMETER_MM
    length = PAPER_S27_FILAMENT_LENGTH_MM
    printed_expression = math.pi * stated_diameter**2 * length / volume
    count_if_diameter = volume / (math.pi * (stated_diameter / 2) ** 2 * length)
    count_if_radius = volume / (math.pi * stated_diameter**2 * length)
    return {
        "supplement_equation": "S27",
        "plasma_volume_mm3_from_product_fitted_fraction": volume,
        "stated_filament_diameter_mm": stated_diameter,
        "stated_filament_length_mm": length,
        "reported_microdischarges_per_cycle_approx": PAPER_S27_REPORTED_MICRODISCHARGES_PER_CYCLE,
        "as_typeset_pi_d_squared_l_over_volume": printed_expression,
        "volume_over_cylinder_if_stated_size_is_diameter": count_if_diameter,
        "volume_over_cylinder_if_stated_size_is_radius": count_if_radius,
        "radius_convention_recovers_reported_count": math.isclose(
            count_if_radius, PAPER_S27_REPORTED_MICRODISCHARGES_PER_CYCLE, rel_tol=0.01
        ),
        "reactor_specific_filament_diameter_measured": False,
        "plasma_volume_independent_of_products": False,
        "allowed_as_exposure_input": False,
        "interpretation": (
            "As typeset, S27 is the inverse of a volume-ratio count; a corrected "
            "ratio with the stated diameter gives about 267, whereas treating "
            "the stated size as a radius gives about 67. This arithmetic "
            "ambiguity cannot establish local filament count or exposure, "
            "especially because the plasma volume came from product fitting."
        ),
    }


def audit_dry_electrical_exposure() -> dict[str, Any]:
    area_m2 = math.pi * (0.007**2 - 0.003**2) / 4
    geometric_current_density_A_m2 = PAPER_HALF_CYCLE_AVERAGE_CURRENT_A / area_m2
    geometric_electron_density_m3 = (
        geometric_current_density_A_m2
        / (elementary_charge * PAPER_DRIFT_VELOCITY_M_S)
    )
    published_j_implied_density_m3 = (
        PAPER_REPORTED_CURRENT_DENSITY_A_M2
        / (elementary_charge * PAPER_DRIFT_VELOCITY_M_S)
    )

    volume_m3 = REACTOR_VOLUME_L * 1e-3
    inventory_mol = PRESSURE_PA * volume_m3 / (R * PAPER_EEDF_GAS_TEMPERATURE_K)
    packet_residence_s = (
        RESIDENCE_S * AMBIENT_FLOW_REFERENCE_K / PAPER_EEDF_GAS_TEMPERATURE_K
    )
    prepared, inlet = _gri_inventory(inventory_mol, 0.01)
    controls = BoltzmannDischargeProgram.model_validate({
        "collision_set": COMMUNITY_CH4_DATASET_ID,
        "reduced_electric_field_Td": PAPER_ELECTRICAL_FIELD_TD,
        "electron_number_density": {
            "value": PAPER_ELECTRON_DENSITY_M3,
            "unit": "1/m^3",
            "source": "Kamarinopoulou et al. 2024 supplement current-based estimate",
            "method": "printed current density divided by charge and BOLSIG+ drift speed",
            "basis": "literature",
        },
        "acknowledge_two_term_limit": True,
    })
    snapshot = solve_mixture_plasma_snapshot(
        controls=controls,
        identities=prepared["identities"],
        gas_amounts_mol=inlet,
        gas_volume_L=REACTOR_VOLUME_L,
        gas_temperature_K=PAPER_EEDF_GAS_TEMPERATURE_K,
        discharge_duration_s=packet_residence_s,
        duty_fraction=1.0,
        deposited_energy_J=AVERAGE_DISCHARGE_POWER_W * packet_residence_s,
    )
    drift = snapshot["drift_power_diagnostic"]
    model_drift_velocity_m_s = (
        drift["electron_mobility_m2_V_s"] * drift["electric_field_V_m"]
    )
    uniform_drift_power_at_published_density_W = drift["estimated_electron_drift_power_W"]
    uniform_drift_power_at_geometric_density_W = (
        uniform_drift_power_at_published_density_W
        * geometric_electron_density_m3 / PAPER_ELECTRON_DENSITY_M3
    )
    return {
        "status": "electrical_arithmetic_and_conditional_power_bound_audited",
        "source": {
            "doi": "10.1126/sciadv.adl4246",
            "location": "supplement, Reduced electric field and electron density measurements, eqs. S23-S26",
            "supplement_url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10980260/supplementaryFiles",
        },
        "paper_printed_values": {
            "emission_based_reduced_field_Td": PAPER_ELECTRICAL_FIELD_TD,
            "half_cycle_average_current_A": PAPER_HALF_CYCLE_AVERAGE_CURRENT_A,
            "annulus_area_cm2": PAPER_ANNULUS_AREA_CM2,
            "current_density_A_m2": PAPER_REPORTED_CURRENT_DENSITY_A_M2,
            "drift_velocity_m_s": PAPER_DRIFT_VELOCITY_M_S,
            "electron_density_m3": PAPER_ELECTRON_DENSITY_M3,
            "average_power_W_scenario": AVERAGE_DISCHARGE_POWER_W,
        },
        "supplement_s27_microdischarge_count_audit": audit_s27_microdischarge_count(),
        "arithmetic_audit": {
            "area_from_printed_diameters_cm2": area_m2 * 1e4,
            "current_density_from_printed_current_and_area_A_m2": geometric_current_density_A_m2,
            "printed_to_recomputed_current_density_ratio": (
                PAPER_REPORTED_CURRENT_DENSITY_A_M2 / geometric_current_density_A_m2
            ),
            "electron_density_from_printed_current_and_area_m3": geometric_electron_density_m3,
            "electron_density_from_printed_current_density_m3": published_j_implied_density_m3,
            "printed_current_area_density_arithmetic_consistent": math.isclose(
                geometric_current_density_A_m2,
                PAPER_REPORTED_CURRENT_DENSITY_A_M2,
                rel_tol=0.05,
            ),
            "model_drift_velocity_m_s": model_drift_velocity_m_s,
            "model_to_paper_drift_velocity_ratio": (
                model_drift_velocity_m_s / PAPER_DRIFT_VELOCITY_M_S
            ),
        },
        "conditional_uniform_drift_power": {
            "at_published_electron_density_W": uniform_drift_power_at_published_density_W,
            "at_current_area_recomputed_electron_density_W": uniform_drift_power_at_geometric_density_W,
            "power_capped_uniform_effective_electron_density_m3": snapshot[
                "electron_number_density_m3"
            ],
            "power_compatible_space_time_fraction_upper_if_published_density_local": (
                AVERAGE_DISCHARGE_POWER_W / uniform_drift_power_at_published_density_W
            ),
            "power_compatible_space_time_fraction_upper_if_recomputed_density_local": (
                AVERAGE_DISCHARGE_POWER_W / uniform_drift_power_at_geometric_density_W
            ),
            "all_measured_power_assigned_to_electron_drift_for_upper_bound": True,
            "independent_streamer_occupancy_measured": False,
            "fraction_injected_into_chemistry": False,
            "equivalence_limit": (
                "A uniform power-capped electron density has the same linear "
                "electron-target source rate as a local density times a "
                "space-time filling fraction, but not the same nonlinear "
                "radical chemistry, gradients, or repeated streamer exposure."
            ),
        },
        "hcn_assay_loaded": False,
        "interpretation": (
            "The printed current, annulus area, current density and inferred "
            "electron density do not all agree arithmetically. If either "
            "density were a local value at the paper's field, power limits "
            "its joint space-time filling fraction under a uniform drift "
            "approximation. But the published density itself assumes a "
            "homogeneous area, AC current need not equal electron conduction "
            "current, and local streamer field/density are unknown. These "
            "fractions are conditional bounds, not measured exposure."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_electrical_exposure(), indent=2))
