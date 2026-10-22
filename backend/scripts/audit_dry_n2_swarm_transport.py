"""Assay-free pure-N2 electron drift check for the dry-flow EEDF machinery.

Experimental points are a small transcription of Dutton (1975), table 1.17,
page 621, which identifies the individual original measurements. They are
never passed to the EEDF solver or used to choose the dry reactor's field.
Historical overlap with the imported Phelps cross-section fit is unknown.
"""
from __future__ import annotations

import json
from typing import Any

from app.prebiotic.boltzmann import COMMUNITY_CH4_DATASET_ID, solve_self_consistent_eedf
from app.prebiotic.schema import BoltzmannEEDFRequest


# Dutton table units: E/N in 10^-17 V cm2 (one Td), W in 10^7 cm/s
# (one such unit is 10^5 m/s). Source codes retain the table's attribution.
_HELD_OUT_N2_POINTS = (
    (164.2, 1.58, "d"),
    (208.0, 1.96, "f"),
    (212.2, 2.24, "b"),
    (242.5, 2.53, "b"),
    (243.0, 2.23, "f"),
    (330.0, 3.30, "d"),
    (333.0, 3.10, "g"),
    (446.0, 4.50, "g"),
    (450.0, 4.50, "d"),
)

_STUDIES = {
    "b": "H. A. Blevin et al., Australian Journal of Physics 20, 741 (1967), DOI 10.1071/PH670741",
    "d": "K. H. Wagner, Zeitschrift fuer Physik 178, 64 (1964)",
    "f": "L. Frommhold, Zeitschrift fuer Physik 160, 554 (1960)",
    "g": "K. H. Wagner and H. Raether, Zeitschrift fuer Physik 170, 540 (1962)",
}


def _model_pure_n2_drift_velocity(
    field_Td: float, *, energy_max_eV: float = 150.0, energy_cells: int = 1000,
) -> dict[str, Any]:
    result = solve_self_consistent_eedf(BoltzmannEEDFRequest(
        collision_set=COMMUNITY_CH4_DATASET_ID,
        mole_fractions={"N2": 1.0},
        gas_temperature_K=295.0,
        reduced_electric_field_Td=field_Td,
        energy_max_eV=energy_max_eV,
        energy_cells=energy_cells,
        acknowledge_two_term_limit=True,
    ))
    reduced_mobility = float(result["transport"]["reduced_electron_mobility_m_inv_V_inv_s_inv"])
    return {
        "field_Td": field_Td,
        "energy_max_eV": energy_max_eV,
        "energy_cells": energy_cells,
        "eedf_converged": result["solver"]["converged"],
        # (mu*N) * (E/N) = mu*E, independent of the chosen gas number density.
        "model_two_term_drift_m_s": reduced_mobility * field_Td * 1e-21,
    }


def audit_dry_n2_swarm_transport() -> dict[str, Any]:
    rows = []
    for field_Td, measured_units, study in _HELD_OUT_N2_POINTS:
        prediction = _model_pure_n2_drift_velocity(field_Td)
        measured_m_s = measured_units * 1e5
        rows.append({
            **prediction,
            "measured_drift_m_s": measured_m_s,
            "dutton_table_study_code": study,
            "model_minus_measured_fraction": (
                prediction["model_two_term_drift_m_s"] / measured_m_s - 1
            ),
        })
    high_field_grid = [
        _model_pure_n2_drift_velocity(450.0, energy_max_eV=maximum, energy_cells=cells)
        for maximum, cells in ((60.0, 400), (100.0, 700), (150.0, 1000))
    ]
    return {
        "status": "external_pure_n2_transport_residuals_not_discharge_calibration",
        "reference": {
            "title": "A Survey of Electron Swarm Data",
            "author": "J. Dutton",
            "year": 1975,
            "url": "https://srd.nist.gov/JPCRD/jpcrd69.pdf",
            "location": "printed page 621, table 1.17",
            "selected_point_count": len(_HELD_OUT_N2_POINTS),
            "original_experimental_studies_by_table_code": _STUDIES,
            "temperature_basis": (
                "Dutton labels the table's N2 data at ambient temperatures; "
                "295 K is a declared model representative, not every study's measured T."
            ),
            "phelps_swarm_fit_training_overlap_ruled_out": False,
        },
        "collision_set": COMMUNITY_CH4_DATASET_ID,
        "pure_n2_295K_comparison": rows,
        "pure_n2_450Td_grid_check": high_field_grid,
        "measured_drift_used_as_solver_input": False,
        "field_selected_from_hcn_or_products": False,
        "hcn_assay_loaded": False,
        "revised_dry_reactor_exposure_from_survey": False,
        "interpretation": (
            "The imported N2 collision set and two-term mobility underpredict "
            "several high-field historical electron-swarm drift values. These "
            "measurements are pure N2 near ambient temperature, whereas the "
            "HCN target is hot 1% CH4/N2. Their measured drift-speed "
            "definitions and ionization effects may not match the BOLOS "
            "mobility-derived speed. Some historical "
            "swarm data may also have informed the Phelps set, so this is an "
            "external residual check rather than independent validation of "
            "local streamer field, dissociation branching, or HCN yield."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_n2_swarm_transport(), indent=2))
