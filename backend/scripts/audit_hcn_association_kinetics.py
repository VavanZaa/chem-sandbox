"""Read-only, assay-free applicability check for H + HCN association.

Tsang and Herron (1991), J. Phys. Chem. Ref. Data 20, 609-663,
DOI 10.1063/1.555890, pp. 654-655, explicitly call their H + HCN
association expressions estimates rather than measurements. This audit does
not replace GRI kinetics or infer a wet-flask HCN yield.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import cantera as ct
from scipy.constants import Avogadro, Boltzmann


FEED_CARD = (
    Path(__file__).resolve().parents[1] / "data" / "benchmarks"
    / "barth_2024_feed_only_220725001_v1.json"
)
GRI_ASSOCIATION_INDEX = 236
GRI_ASSOCIATION_EQUATION = "H + HCN (+M) <=> H2CN (+M)"
SOURCE_DOI = "10.1063/1.555890"


def audit_hcn_association_initial_feed(feed_path: Path = FEED_CARD) -> dict[str, Any]:
    """Compare feed-state GRI with the NIST estimate's distinct rate limits.

    The evaluated k0 terms cover N2 and CO2 only. Their mixture contribution
    is *not* a complete low-pressure rate because H2, H2O, O2, and possible
    other colliders lack terms in this intake. No falloff rate is inferred.
    """

    feed = json.loads(feed_path.read_text(encoding="utf-8"))
    assumptions = feed["diagnostic_only_assumptions"]
    temperature_K = float(assumptions["gas_and_water_temperature_K"])
    pressure_Pa = float(assumptions["initial_pressure_Pa"])
    if not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("Feed temperature must be finite and positive")
    if not math.isfinite(pressure_Pa) or pressure_Pa <= 0:
        raise ValueError("Feed pressure must be finite and positive")
    fractions = {
        name: float(value)
        for name, value in feed["measured_initial_gas_mole_fractions"].items()
    }
    if any(not math.isfinite(value) or value < 0 for value in fractions.values()):
        raise ValueError("Feed mole fractions must be finite and nonnegative")
    nitrogen_fraction = 1.0 - sum(fractions.values())
    if not 0 < nitrogen_fraction < 1:
        raise ValueError("Reported feed fractions leave no positive N2 balance")
    fractions["N2"] = nitrogen_fraction

    phase = ct.Solution("gri30.yaml")
    if phase.reaction_equations()[GRI_ASSOCIATION_INDEX] != GRI_ASSOCIATION_EQUATION:
        raise ValueError("The installed GRI H + HCN reaction index has changed")
    if any(name not in phase.species_names for name in fractions):
        raise ValueError("GRI lacks a feed species needed for the comparison")
    phase.TPX = temperature_K, pressure_Pa, fractions
    gri_effective_rate = (
        float(phase.forward_rate_constants[GRI_ASSOCIATION_INDEX])
        * 1000.0 / Avogadro
    )
    # Tsang & Herron, pp. 654-655.  k_inf has cm3 molecule-1 s-1 units;
    # k0 has cm6 molecule-2 s-1 units.  The source specifies 300-2500 K.
    tsang_high_pressure_limit = 5.5e-11 * math.exp(-2438.0 / temperature_K)
    low_pressure_coefficients = {
        "N2": 4.4e-24 * temperature_K ** -2.73 * math.exp(-3855.0 / temperature_K),
        "CO2": 4.3e-24 * temperature_K ** -2.63 * math.exp(-3791.0 / temperature_K),
    }
    number_density_molecule_cm3 = pressure_Pa / (Boltzmann * temperature_K) / 1e6
    reviewed_terms = {
        name: coefficient * number_density_molecule_cm3 * fractions[name]
        for name, coefficient in low_pressure_coefficients.items()
    }
    source_range_K = [300.0, 2500.0]
    phase.TPX = source_range_K[0], pressure_Pa, fractions
    gri_at_source_minimum = (
        float(phase.forward_rate_constants[GRI_ASSOCIATION_INDEX])
        * 1000.0 / Avogadro
    )
    tsang_high_pressure_at_source_minimum = (
        5.5e-11 * math.exp(-2438.0 / source_range_K[0])
    )
    return {
        "status": "read_only_unvalidated_kinetic_applicability_audit",
        "source": "Tsang and Herron 1991, J. Phys. Chem. Ref. Data 20, 609-663, pp. 654-655",
        "source_doi": SOURCE_DOI,
        "source_url": "https://srd.nist.gov/jpcrdreprint/1.555890.pdf",
        "source_rate_type": "evaluated estimate; authors state no direct rate measurements",
        "source_temperature_range_K": source_range_K,
        "declared_feed_temperature_K": temperature_K,
        "declared_feed_pressure_Pa": pressure_Pa,
        "source_temperature_in_range": source_range_K[0] <= temperature_K <= source_range_K[1],
        "source_temperature_extrapolation_below_minimum_K": max(
            source_range_K[0] - temperature_K, 0.0
        ),
        "initial_N2_balance_mole_fraction": nitrogen_fraction,
        "gri_effective_rate_cm3_molecule_s": gri_effective_rate,
        "tsang_estimated_high_pressure_limit_cm3_molecule_s": tsang_high_pressure_limit,
        "gri_over_tsang_estimated_high_pressure_limit": (
            gri_effective_rate / tsang_high_pressure_limit
        ),
        "source_minimum_temperature_reference": {
            "temperature_K": source_range_K[0],
            "pressure_Pa": pressure_Pa,
            "same_hypothetical_initial_composition": True,
            "gri_effective_rate_cm3_molecule_s": gri_at_source_minimum,
            "tsang_estimated_high_pressure_limit_cm3_molecule_s": (
                tsang_high_pressure_at_source_minimum
            ),
            "gri_over_tsang_estimated_high_pressure_limit": (
                gri_at_source_minimum / tsang_high_pressure_at_source_minimum
            ),
        },
        "tsang_estimated_k0_cm6_molecule2_s_by_reviewed_collider": low_pressure_coefficients,
        "tsang_reviewed_collider_k0_times_density_cm3_molecule_s": reviewed_terms,
        "tsang_reviewed_collider_k0_times_density_sum_cm3_molecule_s": sum(
            reviewed_terms.values()
        ),
        "reviewed_collider_terms_are_complete_mixture_rate": False,
        "missing_collider_terms": sorted(set(fractions) - set(low_pressure_coefficients)),
        "tsang_high_pressure_limit_is_a_measured_bound": False,
        "tsang_stated_low_temperature_uncertainty": (
            "factor 100 at 500 K; not a quantified confidence interval at 298.15 K"
        ),
        "gri_mechanism_modified": False,
        "experimental_hcn_assay_loaded": False,
        "interpretation": (
            "This compares an installed GRI effective rate with an estimated "
            "high-pressure limit and incomplete N2/CO2 low-pressure terms. "
            "It is not a source-supported full wet-feed falloff law, a calibrated "
            "uncertainty interval, or an HCN yield."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_hcn_association_initial_feed(), indent=2, sort_keys=True))
