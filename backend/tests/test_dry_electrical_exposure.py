"""Electrical audit uses primary-source numbers, never product assays."""

import pytest

from scripts.audit_dry_electrical_exposure import (
    audit_dry_electrical_exposure,
    audit_s27_microdischarge_count,
)


def test_s27_count_is_ambiguous_and_not_independent_exposure() -> None:
    report = audit_s27_microdischarge_count()
    assert report["as_typeset_pi_d_squared_l_over_volume"] == pytest.approx(0.01496, rel=1e-3)
    assert report["volume_over_cylinder_if_stated_size_is_diameter"] == pytest.approx(267.38, rel=1e-3)
    assert report["volume_over_cylinder_if_stated_size_is_radius"] == pytest.approx(66.845, rel=1e-3)
    assert report["radius_convention_recovers_reported_count"] is True
    assert report["plasma_volume_independent_of_products"] is False
    assert report["reactor_specific_filament_diameter_measured"] is False
    assert report["allowed_as_exposure_input"] is False


def test_dry_electrical_numbers_show_conditional_not_measured_exposure() -> None:
    report = audit_dry_electrical_exposure()
    arithmetic = report["arithmetic_audit"]
    power = report["conditional_uniform_drift_power"]
    assert report["hcn_assay_loaded"] is False
    assert report["supplement_s27_microdischarge_count_audit"] == audit_s27_microdischarge_count()
    assert arithmetic["area_from_printed_diameters_cm2"] == pytest.approx(0.314, rel=1e-3)
    assert arithmetic["printed_current_area_density_arithmetic_consistent"] is False
    assert arithmetic["printed_to_recomputed_current_density_ratio"] > 4
    assert arithmetic["electron_density_from_printed_current_density_m3"] == pytest.approx(
        report["paper_printed_values"]["electron_density_m3"], rel=0.02
    )
    assert arithmetic["model_to_paper_drift_velocity_ratio"] == pytest.approx(1, rel=0.1)
    assert power["at_published_electron_density_W"] > report["paper_printed_values"]["average_power_W_scenario"]
    assert 0 < power["power_compatible_space_time_fraction_upper_if_published_density_local"] < (
        power["power_compatible_space_time_fraction_upper_if_recomputed_density_local"] < 1
    )
    assert power["power_capped_uniform_effective_electron_density_m3"] == pytest.approx(
        report["paper_printed_values"]["electron_density_m3"]
        * power["power_compatible_space_time_fraction_upper_if_published_density_local"]
    )
    assert power["independent_streamer_occupancy_measured"] is False
    assert power["fraction_injected_into_chemistry"] is False
