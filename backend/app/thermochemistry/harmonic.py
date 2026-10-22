from __future__ import annotations

import math
from typing import Any

from .calculator import IMAGINARY_THRESHOLD_CM_INV, R_KJ_MOL_K
from .schema import HarmonicCorrectionRequest


# N_A h c in kJ mol-1 cm, so multiplying by a wavenumber in cm-1
# gives the molar quantum energy.
WAVENUMBER_TO_KJ_MOL = 0.0119626565638697


def harmonic_vibrational_corrections(request: HarmonicCorrectionRequest) -> dict[str, Any]:
    temperature_K = request.temperature.as_unit("K")
    if temperature_K <= 0:
        raise ValueError("Temperature must be greater than 0 K")
    if request.frequencies.basis == "hypothetical":
        raise ValueError("Harmonic corrections require computed, experimental, or literature frequencies")
    imaginary = [value for value in request.frequencies.values_cm_inv if value < IMAGINARY_THRESHOLD_CM_INV]
    required_imaginary = 1 if request.stationary_point_kind == "transition_state" else 0
    if len(imaginary) != required_imaginary:
        raise ValueError(
            f"{request.stationary_point_kind} requires {required_imaginary} significant imaginary mode(s); found {len(imaginary)}"
        )
    positive = [value for value in request.frequencies.values_cm_inv if value > 0]
    if not positive:
        raise ValueError("At least one positive vibrational mode is required")
    treated = []
    replacements = []
    for value in positive:
        if request.low_frequency_treatment == "raise_to_cutoff" and value < request.low_frequency_cutoff_cm_inv:
            treated.append(request.low_frequency_cutoff_cm_inv)
            replacements.append({"original_cm_inv": value, "used_cm_inv": request.low_frequency_cutoff_cm_inv})
        else:
            treated.append(value)
    quanta = [value * WAVENUMBER_TO_KJ_MOL for value in treated]
    zpe = 0.5 * sum(quanta)
    thermal_energy = 0.0
    entropy_kj_mol_K = 0.0
    for quantum in quanta:
        x = quantum / (R_KJ_MOL_K * temperature_K)
        occupation = 0.0 if x > 700 else 1.0 / math.expm1(x)
        thermal_energy += quantum * occupation
        entropy_kj_mol_K += R_KJ_MOL_K * (
            x * occupation - math.log1p(-math.exp(-x))
        )
    return {
        "temperature_K": temperature_K,
        "stationary_point_kind": request.stationary_point_kind,
        "zero_point_energy_kj_mol": zpe,
        "vibrational_thermal_energy_excluding_zpe_kj_mol": thermal_energy,
        "vibrational_entropy_J_mol_K": entropy_kj_mol_K * 1000.0,
        "significant_imaginary_modes_cm_inv": imaginary,
        "positive_modes_used_cm_inv": treated,
        "low_frequency_replacements": replacements,
        "source": request.frequencies.source,
        "method": request.frequencies.method,
        "limitations": [
            "These are vibrational harmonic-oscillator terms only.",
            "Translational, rotational, electronic, solvation, and standard-state terms are not included.",
            "Low-frequency cutoff treatment is an approximation and is recorded explicitly.",
        ],
    }
