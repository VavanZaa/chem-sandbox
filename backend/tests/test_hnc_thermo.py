"""Do not silently mix the anomalous published HNC record into GRI."""

import math

import cantera as ct
import pytest

from app.prebiotic.heavy_gas import KAMARINOPOULOU_MECHANISM
from app.prebiotic.hnc_thermo import (
    NIST_HNC_FORMATION_ENTHALPY_KJ_MOL,
    NIST_HNC_FORMATION_ENTHALPY_UNCERTAINTY_KJ_MOL,
    NIST_HNC_REFERENCE_PRESSURE_PA,
    NIST_HNC_STANDARD_CP_J_MOL_K,
    NIST_HNC_STANDARD_ENTROPY_J_MOL_K,
    nist_hnc_room_temperature_species,
    room_temperature_h2cn_h_branch_enthalpies,
)


def test_nist_hnc_reference_state_and_gri_hcn_isomer_gap() -> None:
    gri = ct.Solution("gri30.yaml")
    assert "HNC" not in gri.species_names
    hnc = nist_hnc_room_temperature_species(298.15)
    phase = ct.Solution(
        thermo="ideal-gas", kinetics="gas",
        species=[gri.species("HCN"), hnc], reactions=[],
    )
    phase.TPX = 298.15, ct.one_atm, {"HNC": 1}
    assert phase.enthalpy_mole / 1e6 == pytest.approx(NIST_HNC_FORMATION_ENTHALPY_KJ_MOL)
    assert phase.entropy_mole / 1e3 < NIST_HNC_STANDARD_ENTROPY_J_MOL_K
    assert phase.cp_mole / 1e3 == pytest.approx(NIST_HNC_STANDARD_CP_J_MOL_K)
    phase.TPX = 298.15, NIST_HNC_REFERENCE_PRESSURE_PA, {"HNC": 1}
    assert phase.entropy_mole / 1e3 == pytest.approx(NIST_HNC_STANDARD_ENTROPY_J_MOL_K)
    phase.TPX = 298.15, ct.one_atm, {"HCN": 1}
    isomer_gap_kj_mol = NIST_HNC_FORMATION_ENTHALPY_KJ_MOL - phase.enthalpy_mole / 1e6
    assert isomer_gap_kj_mol == pytest.approx(63.51477513986397)
    assert isomer_gap_kj_mol > NIST_HNC_FORMATION_ENTHALPY_UNCERTAINTY_KJ_MOL
    # Independent ATcT v1.122 formation enthalpies: HNC 192.18, HCN 129.277
    # kJ/mol at 298.15 K. This is a cross-check, never an apparatus rate/yield.
    atct_isomer_gap_kj_mol = 192.18 - 129.277
    assert abs(isomer_gap_kj_mol - atct_isomer_gap_kj_mol) < (
        NIST_HNC_FORMATION_ENTHALPY_UNCERTAINTY_KJ_MOL
    )
    assert hnc.input_data["reference"].startswith("https://cccbdb.nist.gov/")


def test_nist_hnc_constant_cp_cannot_be_used_outside_narrow_window() -> None:
    for temperature in (0.0, 289.99, 310.01, 1000.0, math.nan, math.inf):
        with pytest.raises(ValueError, match="restricted to 290–310 K"):
            nist_hnc_room_temperature_species(temperature)
    assert nist_hnc_room_temperature_species(290.0).name == "HNC"
    assert nist_hnc_room_temperature_species(310.0).name == "HNC"


def test_published_hnc_record_disagrees_with_evaluated_ground_state() -> None:
    published = ct.Solution(str(KAMARINOPOULOU_MECHANISM))
    published.TPX = 298.15, ct.one_atm, {"HNC": 1}
    published_hnc_enthalpy_kj_mol = published.enthalpy_mole / 1e6
    published_hnc_cp_j_mol_k = published.cp_mole / 1e3
    published_hnc_entropy_j_mol_k = published.entropy_mole / 1e3
    assert published_hnc_enthalpy_kj_mol == pytest.approx(638.1481659734416)
    assert published_hnc_enthalpy_kj_mol - NIST_HNC_FORMATION_ENTHALPY_KJ_MOL > 400
    assert abs(published_hnc_cp_j_mol_k - NIST_HNC_STANDARD_CP_J_MOL_K) > 5
    assert abs(published_hnc_entropy_j_mol_k - NIST_HNC_STANDARD_ENTROPY_J_MOL_K) > 20


def test_h2cn_h_exit_channels_are_exothermic_and_isomer_consistent() -> None:
    result = room_temperature_h2cn_h_branch_enthalpies(298.15)
    hcn = result["h2cn_plus_h_to_hcn_plus_h2_delta_h_kj_mol"]
    hnc = result["h2cn_plus_h_to_hnc_plus_h2_delta_h_kj_mol"]
    assert hcn == pytest.approx(-334.51497283279406)
    assert hnc == pytest.approx(-271.0001976929301)
    assert hcn < hnc < 0
    assert hnc - hcn == pytest.approx(result["hnc_minus_hcn_isomer_gap_kj_mol"])
    with pytest.raises(ValueError, match="restricted to 290–310 K"):
        room_temperature_h2cn_h_branch_enthalpies(400.0)
