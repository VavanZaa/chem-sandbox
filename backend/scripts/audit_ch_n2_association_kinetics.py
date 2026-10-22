"""Assay-free applicability audit for GRI's CH + N2 association entry.

Brownsword et al., J. Chem. Soc., Faraday Trans. 92, 1087 (1996),
DOI 10.1039/FT9969201087, measured pressure-dependent CH + N2 loss over
202-584 K in argon and report a fitted Ar low-pressure association limit.
The measurement is predominantly association, not an HCNN product assay or
a measured wet-gas effective rate, product split, or HCN yield.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import cantera as ct
from scipy.constants import Avogadro

from scripts.probe_barth_feed_pathway import FEED_CARD


GRI_ASSOCIATION_INDEX = 240
GRI_ASSOCIATION_EQUATION = "CH + N2 (+M) <=> HCNN (+M)"
SOURCE_DOI = "10.1039/FT9969201087"
HIGH_PRESSURE_SOURCE_DOI = "10.1063/1.472383"


def audit_ch_n2_association(feed_path: Path = FEED_CARD) -> dict[str, Any]:
    feed = json.loads(feed_path.read_text(encoding="utf-8"))
    assumptions = feed["diagnostic_only_assumptions"]
    temperature_K = float(assumptions["gas_and_water_temperature_K"])
    pressure_Pa = float(assumptions["initial_pressure_Pa"])
    if not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("Feed temperature must be finite and positive")
    if not math.isfinite(pressure_Pa) or pressure_Pa <= 0:
        raise ValueError("Feed pressure must be finite and positive")
    fractions = {
        species: float(fraction)
        for species, fraction in feed["measured_initial_gas_mole_fractions"].items()
    }
    if any(not math.isfinite(value) or value < 0 for value in fractions.values()):
        raise ValueError("Feed gas fractions must be finite and nonnegative")
    fractions["N2"] = 1.0 - sum(fractions.values())
    if fractions["N2"] <= 0:
        raise ValueError("Feed fractions leave no positive N2 balance")

    mechanism_file = next(
        (
            Path(directory) / "gri30.yaml"
            for directory in ct.get_data_directories()
            if (Path(directory) / "gri30.yaml").is_file()
        ),
        None,
    )
    if mechanism_file is None:
        raise ValueError("Cantera's packaged GRI mechanism is unavailable")
    phase = ct.Solution(str(mechanism_file))
    reaction = phase.reaction(GRI_ASSOCIATION_INDEX)
    if reaction.equation != GRI_ASSOCIATION_EQUATION:
        raise ValueError("The installed GRI CH + N2 reaction index has changed")
    if reaction.third_body.efficiencies.get("AR") != 1.0:
        raise ValueError("The GRI argon third-body efficiency has changed")

    # Cantera low_rate uses m6 kmol-2 s-1; one kmol contains 1000*N_A
    # molecules and one m6 is 1e12 cm6, hence the factor 1e6/N_A^2.
    gri_low_pressure_Ar_cm6_molecule2_s = (
        float(reaction.rate.low_rate(temperature_K)) * 1e6 / Avogadro**2
    )
    # Brownsword et al., abstract: k0,Ar(T) = 1.7e-31 (T/298)^-2.3
    # cm6 molecule-2 s-1 for the association-dominated low-P branch.
    measured_fitted_Ar_low_pressure_cm6_molecule2_s = (
        1.7e-31 * (temperature_K / 298.0) ** -2.3
    )
    # Fulle and Hippler, J. Chem. Phys. 105, 5423 (1996), extracted this
    # asymptotic limit from 1-150 bar helium experiments at 200-715 K.
    # It is an experimentally fitted high-P limit, not a direct wet 1-bar rate.
    gri_high_pressure_cm3_molecule_s = (
        float(reaction.rate.high_rate(temperature_K)) * 1000.0 / Avogadro
    )
    fulle_high_pressure_cm3_molecule_s = (
        4.1e-11 * (temperature_K / 300.0) ** -0.15
    )
    fulle_high_pressure_reported_coefficient_half_width_cm3_molecule_s = (
        0.8e-11 * (temperature_K / 300.0) ** -0.15
    )
    phase.TPX = temperature_K, pressure_Pa, fractions
    effective_third_body_kmol_m3 = sum(
        float(phase.concentrations[index]) * reaction.third_body.efficiencies.get(
            species, reaction.third_body.default_efficiency,
        )
        for index, species in enumerate(phase.species_names)
    )
    gri_wet_feed_low_pressure_projection_cm3_molecule_s = (
        float(reaction.rate.low_rate(temperature_K))
        * effective_third_body_kmol_m3 * 1000.0 / Avogadro
    )
    gri_wet_feed_effective_cm3_molecule_s = (
        float(phase.forward_rate_constants[GRI_ASSOCIATION_INDEX])
        * 1000.0 / Avogadro
    )
    return {
        "status": "read_only_argon_low_and_helium_high_pressure_association_comparison",
        "source": "Brownsword et al., J. Chem. Soc., Faraday Trans. 92, 1087-1094 (1996)",
        "source_doi": SOURCE_DOI,
        "source_url": "https://doi.org/10.1039/FT9969201087",
        "source_temperature_range_K": [202.0, 584.0],
        "source_argon_pressure_range_bar": [0.00533, 0.53],
        "source_argon_pressure_range_reference": (
            "https://kinetics.nist.gov/kinetics/Detail?id=1996BRO%2FHER1087-1094%3A2"
        ),
        "high_pressure_source": "Fulle and Hippler, J. Chem. Phys. 105, 5423-5430 (1996)",
        "high_pressure_source_doi": HIGH_PRESSURE_SOURCE_DOI,
        "high_pressure_source_helium_range_bar": [1.0, 150.0],
        "high_pressure_source_temperature_range_K": [200.0, 715.0],
        "high_pressure_source_falloff_curve_temperatures_K": [200.0, 250.0, 300.0, 400.0, 500.0],
        "declared_feed_temperature_in_high_pressure_source_range": 200.0 <= temperature_K <= 715.0,
        "declared_feed_temperature_K": temperature_K,
        "source_temperature_in_range": 202.0 <= temperature_K <= 584.0,
        "declared_feed_pressure_Pa": pressure_Pa,
        "declared_feed_pressure_above_source_argon_measurement_range": (
            pressure_Pa > 0.53e5
        ),
        "cantera_version": ct.__version__,
        "gri_file_sha256": hashlib.sha256(mechanism_file.read_bytes()).hexdigest(),
        "gri_reaction_index": GRI_ASSOCIATION_INDEX,
        "gri_reaction_equation": GRI_ASSOCIATION_EQUATION,
        "gri_argon_low_pressure_limit_cm6_molecule2_s": (
            gri_low_pressure_Ar_cm6_molecule2_s
        ),
        "brownsword_fitted_argon_low_pressure_limit_cm6_molecule2_s": (
            measured_fitted_Ar_low_pressure_cm6_molecule2_s
        ),
        "gri_over_brownsword_argon_low_pressure_limit": (
            gri_low_pressure_Ar_cm6_molecule2_s
            / measured_fitted_Ar_low_pressure_cm6_molecule2_s
        ),
        "gri_high_pressure_limit_cm3_molecule_s": gri_high_pressure_cm3_molecule_s,
        "fulle_fitted_high_pressure_limit_cm3_molecule_s": fulle_high_pressure_cm3_molecule_s,
        "fulle_fitted_high_pressure_reported_coefficient_half_width_cm3_molecule_s": (
            fulle_high_pressure_reported_coefficient_half_width_cm3_molecule_s
        ),
        "gri_over_fulle_high_pressure_limit": (
            gri_high_pressure_cm3_molecule_s / fulle_high_pressure_cm3_molecule_s
        ),
        "gri_wet_feed_effective_rate_cm3_molecule_s": (
            gri_wet_feed_effective_cm3_molecule_s
        ),
        "gri_wet_feed_effective_over_gri_high_pressure_limit": (
            gri_wet_feed_effective_cm3_molecule_s / gri_high_pressure_cm3_molecule_s
        ),
        "gri_wet_feed_effective_third_body_kmol_m3": effective_third_body_kmol_m3,
        "gri_wet_feed_low_pressure_projection_cm3_molecule_s": (
            gri_wet_feed_low_pressure_projection_cm3_molecule_s
        ),
        "gri_wet_feed_falloff_rate_over_low_pressure_projection": (
            gri_wet_feed_effective_cm3_molecule_s
            / gri_wet_feed_low_pressure_projection_cm3_molecule_s
        ),
        "argon_low_pressure_comparison_validates_wet_feed_effective_rate": False,
        "brownsword_global_fit_high_pressure_values_are_independent_validation": False,
        "fulle_high_pressure_limit_validates_wet_feed_effective_rate": False,
        "association_rate_validates_hcnn_o_product_branch": False,
        "single_hcnn_adduct_represents_all_ch_n2_isomers": False,
        "competing_adduct_theory_doi": "10.1021/acs.jpca.6c01935",
        "association_rate_validates_hcn_yield": False,
        "gri_mechanism_modified": False,
        "experimental_hcn_assay_loaded": False,
        "limitation": (
            "The primary measurement constrains Ar-supported CH+N2 association "
            "in its low-pressure limit. The authors themselves describe the "
            "high-pressure global-fit values as fitting parameters, not true "
            "rate constants. An independent He-bath high-pressure experiment "
            "constrains the other asymptotic limit. Neither comparison "
            "determines the 1-bar wet N2/CO2/H2/H2O/O2 falloff, the "
            "HCNN/HNCN adduct assignment and fate, the CH+N2 NCN branch, "
            "or the HCN amount."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_ch_n2_association(), indent=2, sort_keys=True))
