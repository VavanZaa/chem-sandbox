"""Measured Ar association is not silently generalized to a wet HCN yield."""

import pytest

from scripts.audit_ch_n2_association_kinetics import audit_ch_n2_association


def test_ch_n2_argon_low_pressure_audit_is_assay_free_and_bounded() -> None:
    report = audit_ch_n2_association()
    assert report["source_temperature_in_range"] is True
    assert report["declared_feed_pressure_above_source_argon_measurement_range"] is True
    assert report["gri_reaction_equation"] == "CH + N2 (+M) <=> HCNN (+M)"
    assert report["gri_argon_low_pressure_limit_cm6_molecule2_s"] == pytest.approx(
        1.5588480074404923e-31, rel=1e-7,
    )
    assert 0.8 < report["gri_over_brownsword_argon_low_pressure_limit"] < 1.2
    assert report["gri_wet_feed_effective_rate_cm3_molecule_s"] > 0
    assert report["gri_wet_feed_falloff_rate_over_low_pressure_projection"] == pytest.approx(
        0.5615992983371366, rel=1e-7,
    )
    assert report["high_pressure_source_doi"] == "10.1063/1.472383"
    assert report["declared_feed_temperature_in_high_pressure_source_range"] is True
    assert report["gri_high_pressure_limit_cm3_molecule_s"] == pytest.approx(
        1.2099739309968349e-11, rel=1e-7,
    )
    assert report["fulle_fitted_high_pressure_limit_cm3_molecule_s"] == pytest.approx(
        4.1038060072938136e-11, rel=1e-7,
    )
    assert report["gri_over_fulle_high_pressure_limit"] == pytest.approx(
        0.29484189282980555, rel=1e-7,
    )
    assert report["gri_wet_feed_effective_over_gri_high_pressure_limit"] == pytest.approx(
        0.2015216181584579, rel=1e-7,
    )
    assert report["argon_low_pressure_comparison_validates_wet_feed_effective_rate"] is False
    assert report["brownsword_global_fit_high_pressure_values_are_independent_validation"] is False
    assert report["fulle_high_pressure_limit_validates_wet_feed_effective_rate"] is False
    assert report["association_rate_validates_hcnn_o_product_branch"] is False
    assert report["single_hcnn_adduct_represents_all_ch_n2_isomers"] is False
    assert report["competing_adduct_theory_doi"] == "10.1021/acs.jpca.6c01935"
    assert report["association_rate_validates_hcn_yield"] is False
    assert report["gri_mechanism_modified"] is False
    assert report["experimental_hcn_assay_loaded"] is False
    assert len(report["gri_file_sha256"]) == 64
