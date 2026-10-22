"""Assay-free molecular-diffusion timescale check for the dry DBD annulus.

This is a kinetic-theory transport estimate at a uniform scenario state, not
a measurement of filament motion, turbulence, convection, or gas residence.
"""
from __future__ import annotations

import json
import math
from typing import Any

import cantera as ct

from scripts.run_dry_flow_hcn_probe import (
    AC_FREQUENCY_HZ,
    AMBIENT_FLOW_REFERENCE_K,
    DOWNSTREAM_PROBE_TEMPERATURE_K,
    PRESSURE_PA,
    RESIDENCE_S,
)


QUARTZ_INNER_DIAMETER_M = 0.007
GROUND_ROD_DIAMETER_M = 0.003


def audit_dry_diffusive_mixing() -> dict[str, Any]:
    gas = ct.Solution("gri30.yaml")
    gas.TPX = DOWNSTREAM_PROBE_TEMPERATURE_K, PRESSURE_PA, {"CH4": 0.01, "N2": 0.99}
    methane_diffusion_m2_s = float(gas.mix_diff_coeffs[gas.species_index("CH4")])
    annular_radial_gap_m = (QUARTZ_INNER_DIAMETER_M - GROUND_ROD_DIAMETER_M) / 2
    ac_period_s = 1.0 / AC_FREQUENCY_HZ
    packet_residence_s = (
        RESIDENCE_S * AMBIENT_FLOW_REFERENCE_K / DOWNSTREAM_PROBE_TEMPERATURE_K
    )
    gap_diffusion_time_s = annular_radial_gap_m**2 / (2 * methane_diffusion_m2_s)
    steps = (16, 32, 64, 128)
    return {
        "status": "conditional_molecular_diffusion_timescale_audit",
        "source": {
            "transport_mechanism": "Cantera packaged GRI-Mech 3.0 mixture-averaged transport",
            "cantera_version": ct.__version__,
            "geometry": "Kamarinopoulou et al. 2024 reactor methods: 7 mm quartz ID, 3 mm grounding rod OD",
            "ac_frequency": "Kamarinopoulou et al. 2024 Figure 1 conditions: 22 kHz",
        },
        "scenario": {
            "temperature_K": DOWNSTREAM_PROBE_TEMPERATURE_K,
            "temperature_is_uniform_downstream_probe_proxy": True,
            "pressure_Pa": PRESSURE_PA,
            "feed_CH4_mole_fraction": 0.01,
            "assumed_flow_reference_K": AMBIENT_FLOW_REFERENCE_K,
        },
        "methane_mix_diffusion_coefficient_m2_s": methane_diffusion_m2_s,
        "annular_radial_gap_m": annular_radial_gap_m,
        "ac_period_s": ac_period_s,
        "hot_packet_residence_s": packet_residence_s,
        "one_dimensional_rms_diffusion_per_ac_period_m": math.sqrt(
            2 * methane_diffusion_m2_s * ac_period_s
        ),
        "one_dimensional_rms_diffusion_per_hot_packet_m": math.sqrt(
            2 * methane_diffusion_m2_s * packet_residence_s
        ),
        "one_dimensional_rms_gap_crossing_time_s": gap_diffusion_time_s,
        "ac_periods_per_rms_gap_crossing_time": gap_diffusion_time_s / ac_period_s,
        "rms_diffusion_by_complete_remix_encounter_count_m": {
            str(count): math.sqrt(2 * methane_diffusion_m2_s * packet_residence_s / count)
            for count in steps
        },
        "one_d_rms_diffusion_less_than_gap_per_ac_period": (
            math.sqrt(2 * methane_diffusion_m2_s * ac_period_s)
            < annular_radial_gap_m
        ),
        "one_d_rms_diffusion_less_than_gap_at_16_encounters": (
            math.sqrt(2 * methane_diffusion_m2_s * packet_residence_s / 16)
            < annular_radial_gap_m
        ),
        "measured_spatial_mixing_or_streamer_motion_available": False,
        "hcn_assay_loaded": False,
        "interpretation": (
            "One-dimensional RMS molecular diffusion is a transport scale, not "
            "a full-mixing criterion or measured gas path. In particular, "
            "convection, thermal gradients, filament motion, and electrohydrodynamic "
            "mixing could change encounters; this audit cannot infer their rates."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_diffusive_mixing(), indent=2))
