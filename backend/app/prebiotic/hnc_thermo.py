"""Narrow, evaluated ground-state HNC thermochemistry for mechanism audits.

The constructor adds nothing by itself; the separate opt-in gas diagnostic
can use this species for HNC reaction branches and inventory accounting. It
does not provide a gas-water transfer pathway or a validated HNC yield. The
NIST CCCBDB values are 298.15 K reference properties; constant heat capacity
is only a local approximation around room temperature.
"""

from __future__ import annotations

import math

import cantera as ct


NIST_HNC_SOURCE = "https://cccbdb.nist.gov/exp2x.asp?casno=6914074"
NIST_HNC_REFERENCE_TEMPERATURE_K = 298.15
NIST_HNC_REFERENCE_PRESSURE_PA = 100_000.0
GRI_REFERENCE_PRESSURE_PA = ct.one_atm
NIST_HNC_FORMATION_ENTHALPY_KJ_MOL = 194.33
NIST_HNC_FORMATION_ENTHALPY_UNCERTAINTY_KJ_MOL = 9.00
NIST_HNC_STANDARD_ENTROPY_J_MOL_K = 205.21
NIST_HNC_STANDARD_CP_J_MOL_K = 40.03
# This is an intentionally narrow approximation window, not a range reported
# or validated by the reference compilation.
HNC_CONSTANT_CP_DIAGNOSTIC_TEMPERATURE_RANGE_K = (290.0, 310.0)


def _nist_hnc_entropy_at_gri_reference_pressure_j_mol_k() -> float:
    # NIST CCCBDB uses a 1-bar standard pressure, while GRI species in
    # Cantera use 1 atm. For an ideal gas, S(p2) = S(p1) - R ln(p2/p1).
    return NIST_HNC_STANDARD_ENTROPY_J_MOL_K - (
        ct.gas_constant / 1000.0
        * math.log(GRI_REFERENCE_PRESSURE_PA / NIST_HNC_REFERENCE_PRESSURE_PA)
    )


def nist_hnc_room_temperature_species(temperature_K: float) -> ct.Species:
    """Return an HNC species only for a tightly bounded 298 K diagnostic.

    The NIST enthalpy, entropy and heat capacity define the reference point;
    Cantera's constant-Cp expression supplies a local interpolation. Callers
    must not use it for a hot spark zone or broader atmospheric mechanism.
    """

    low, high = HNC_CONSTANT_CP_DIAGNOSTIC_TEMPERATURE_RANGE_K
    if not math.isfinite(temperature_K) or not low <= temperature_K <= high:
        raise ValueError(
            "NIST HNC constant-Cp diagnostic is restricted to 290–310 K"
        )
    species = ct.Species("HNC", {"C": 1, "H": 1, "N": 1})
    species.thermo = ct.ConstantCp(
        low,
        high,
        GRI_REFERENCE_PRESSURE_PA,
        [
            NIST_HNC_REFERENCE_TEMPERATURE_K,
            NIST_HNC_FORMATION_ENTHALPY_KJ_MOL * 1e6,
            _nist_hnc_entropy_at_gri_reference_pressure_j_mol_k() * 1e3,
            NIST_HNC_STANDARD_CP_J_MOL_K * 1e3,
        ],
    )
    species.update_user_data({
        "reference": NIST_HNC_SOURCE,
        "source_compilation": "NIST CCCBDB release 22, HNC, Gurvich et al. 1989",
        "nist_standard_pressure_pa": NIST_HNC_REFERENCE_PRESSURE_PA,
        "cantera_reference_pressure_pa": GRI_REFERENCE_PRESSURE_PA,
        "entropy_reference_pressure_conversion": "ideal gas, 1 bar to 1 atm",
        "approximation": "constant heat capacity, restricted to 290–310 K",
        "not_a_reaction_rate_or_yield": True,
    })
    return species


def room_temperature_h2cn_h_branch_enthalpies(temperature_K: float) -> dict[str, float]:
    """Audit both H2CN + H exit-channel heats without adding either reaction."""

    hnc = nist_hnc_room_temperature_species(temperature_K)
    gri = ct.Solution("gri30.yaml")
    phase = ct.Solution(
        thermo="ideal-gas", kinetics="gas",
        species=[gri.species(name) for name in ("H", "H2", "H2CN", "HCN")] + [hnc],
        reactions=[],
    )
    phase.TP = temperature_K, ct.one_atm
    enthalpies = {
        name: float(phase.partial_molar_enthalpies[phase.species_index(name)] / 1e6)
        for name in phase.species_names
    }
    return {
        "temperature_K": temperature_K,
        "h2cn_plus_h_to_hcn_plus_h2_delta_h_kj_mol": (
            enthalpies["HCN"] + enthalpies["H2"]
            - enthalpies["H2CN"] - enthalpies["H"]
        ),
        "h2cn_plus_h_to_hnc_plus_h2_delta_h_kj_mol": (
            enthalpies["HNC"] + enthalpies["H2"]
            - enthalpies["H2CN"] - enthalpies["H"]
        ),
        "hnc_minus_hcn_isomer_gap_kj_mol": enthalpies["HNC"] - enthalpies["HCN"],
        "scope": "GRI gas thermochemistry plus NIST 298 K HNC constant-Cp diagnostic; no rates or yield",
    }
