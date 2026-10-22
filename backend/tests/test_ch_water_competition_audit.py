"""The CH + water rate check cannot be mistaken for HCN yield validation."""

import pytest

from scripts.audit_ch_water_competition import audit_ch_water_competition


def test_ch_water_initial_feed_competition_is_assay_free_and_nonpropagating() -> None:
    report = audit_ch_water_competition()
    assert report["gri_reaction_equation"] == "CH + H2O <=> CH2O + H"
    assert len(report["gri_file_sha256"]) == 64
    assert report["cantera_version"]
    assert report["declared_feed_pressure_above_measured_room_temperature_range"] is True
    assert report["gri_ch_water_rate_cm3_molecule_s"] == pytest.approx(
        3.390814614272909e-11, rel=1e-7,
    )
    assert report["blitz_measured_fit_ch_water_rate_cm3_molecule_s"] == pytest.approx(
        1.5588856451762595e-11, rel=1e-7,
    )
    assert report["gri_over_blitz_ch_water_rate"] == pytest.approx(
        2.1751528887094973, rel=1e-7,
    )
    assert report["gri_initial_water_over_n2_ch_loss_hazard"] == pytest.approx(
        0.2426068452815605, rel=1e-7,
    )
    assert report["independent_product_study_resolves_CH2O_vs_CO_branch"] is False
    assert report["initial_hazard_ratio_is_trajectory_branching_or_yield"] is False
    assert report["gri_CH2O_product_assignment_independently_verified"] is False
    assert report["experimental_hcn_assay_loaded"] is False
    assert report["gri_mechanism_modified"] is False
