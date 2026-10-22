"""Assay-free initial-feed CH + H2O competition and kinetic-domain audit.

Blitz et al., J. Phys. Chem. A 103, 5699 (1999), DOI 10.1021/jp9903716,
measured CH + H2O disappearance near room temperature. Bergeat et al.,
Chem. Phys. Lett. 480, 21 (2009), DOI 10.1016/j.cplett.2009.07.080,
measured H-atom production, not the CH2O versus CO product split.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import cantera as ct
from scipy.constants import Avogadro

from scripts.audit_ch_n2_association_kinetics import audit_ch_n2_association
from scripts.probe_barth_feed_pathway import FEED_CARD


GRI_CH_WATER_INDEX = 126
GRI_CH_WATER_EQUATION = "CH + H2O <=> CH2O + H"


def audit_ch_water_competition(feed_path: Path = FEED_CARD) -> dict[str, Any]:
    feed = json.loads(feed_path.read_text(encoding="utf-8"))
    assumptions = feed["diagnostic_only_assumptions"]
    temperature_K = float(assumptions["gas_and_water_temperature_K"])
    pressure_Pa = float(assumptions["initial_pressure_Pa"])
    if not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("Feed temperature must be finite and positive")
    if not math.isfinite(pressure_Pa) or pressure_Pa <= 0:
        raise ValueError("Feed pressure must be finite and positive")
    fractions = {
        species: float(value)
        for species, value in feed["measured_initial_gas_mole_fractions"].items()
    }
    if any(not math.isfinite(value) or value < 0 for value in fractions.values()):
        raise ValueError("Feed gas fractions must be finite and nonnegative")
    fractions["N2"] = 1.0 - sum(fractions.values())
    if fractions["N2"] <= 0 or fractions.get("H2O", 0.0) <= 0:
        raise ValueError("Feed must contain positive N2 balance and water vapour")

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
    reaction = phase.reaction(GRI_CH_WATER_INDEX)
    if reaction.equation != GRI_CH_WATER_EQUATION:
        raise ValueError("The installed GRI CH + H2O reaction index has changed")
    phase.TPX = temperature_K, pressure_Pa, fractions
    gri_ch_water_cm3_molecule_s = (
        float(phase.forward_rate_constants[GRI_CH_WATER_INDEX]) * 1000.0 / Avogadro
    )
    # Blitz et al. (1999) measured k = 1.56e-11 (T/298)^-1.42,
    # cm3 molecule-1 s-1; the coefficient's stated width is 0.25e-11.
    blitz_ch_water_cm3_molecule_s = 1.56e-11 * (temperature_K / 298.0) ** -1.42
    blitz_coefficient_half_width_cm3_molecule_s = (
        0.25e-11 * (temperature_K / 298.0) ** -1.42
    )
    n2_audit = audit_ch_n2_association(feed_path)
    gri_ch_n2_wet_cm3_molecule_s = float(
        n2_audit["gri_wet_feed_effective_rate_cm3_molecule_s"]
    )
    # At the same number density, the CH loss hazards are proportional to
    # bimolecular rate times target mole fraction. These are INITIAL-state
    # rates only, not integrated branch fractions or product yields.
    gri_initial_water_over_n2_ch_loss_hazard = (
        gri_ch_water_cm3_molecule_s * fractions["H2O"]
        / (gri_ch_n2_wet_cm3_molecule_s * fractions["N2"])
    )
    return {
        "status": "read_only_initial_feed_ch_water_competition",
        "source": "Blitz et al., J. Phys. Chem. A 103, 5699-5704 (1999)",
        "source_doi": "10.1021/jp9903716",
        "source_room_temperature_pressure_range_Torr": [20.0, 200.0],
        "independent_product_study_doi": "10.1016/j.cplett.2009.07.080",
        "independent_product_study_H_atom_branch_at_330K": 1.0,
        "independent_product_study_resolves_CH2O_vs_CO_branch": False,
        "declared_feed_temperature_K": temperature_K,
        "cantera_version": ct.__version__,
        "gri_file_sha256": hashlib.sha256(mechanism_file.read_bytes()).hexdigest(),
        "declared_feed_pressure_Pa": pressure_Pa,
        "declared_feed_pressure_above_measured_room_temperature_range": (
            pressure_Pa > 200.0 * 101325.0 / 760.0
        ),
        "feed_mole_fraction_H2O": fractions["H2O"],
        "feed_mole_fraction_N2_balance": fractions["N2"],
        "gri_reaction_index": GRI_CH_WATER_INDEX,
        "gri_reaction_equation": GRI_CH_WATER_EQUATION,
        "gri_ch_water_rate_cm3_molecule_s": gri_ch_water_cm3_molecule_s,
        "blitz_measured_fit_ch_water_rate_cm3_molecule_s": blitz_ch_water_cm3_molecule_s,
        "blitz_reported_coefficient_half_width_cm3_molecule_s": (
            blitz_coefficient_half_width_cm3_molecule_s
        ),
        "gri_over_blitz_ch_water_rate": (
            gri_ch_water_cm3_molecule_s / blitz_ch_water_cm3_molecule_s
        ),
        "gri_ch_n2_wet_feed_rate_cm3_molecule_s": gri_ch_n2_wet_cm3_molecule_s,
        "gri_initial_water_over_n2_ch_loss_hazard": (
            gri_initial_water_over_n2_ch_loss_hazard
        ),
        "initial_hazard_ratio_is_trajectory_branching_or_yield": False,
        "gri_CH2O_product_assignment_independently_verified": False,
        "gri_mechanism_modified": False,
        "experimental_hcn_assay_loaded": False,
        "limitation": (
            "The 1999 experiment constrains CH loss by water, not the "
            "CH2O/CO product split. Its room-temperature pressure test "
            "stopped at 200 Torr, below the declared 1-bar feed. The 2009 "
            "H-atom observation does not "
            "resolve that split. The ratio uses only the declared initial "
            "feed and installed GRI kinetics; neither wet-flask time-varying "
            "branching nor an HCN amount is determined."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_ch_water_competition(), indent=2, sort_keys=True))
