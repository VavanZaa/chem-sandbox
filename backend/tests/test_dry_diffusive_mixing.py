"""Molecular diffusion is a bounded scale, not measured streamer mixing."""

import math

import pytest

from scripts.audit_dry_diffusive_mixing import audit_dry_diffusive_mixing


def test_dry_diffusion_scale_does_not_claim_complete_mixing() -> None:
    audit = audit_dry_diffusive_mixing()
    diffusion = audit["methane_mix_diffusion_coefficient_m2_s"]
    period = audit["ac_period_s"]
    gap = audit["annular_radial_gap_m"]
    assert audit["hcn_assay_loaded"] is False
    assert audit["scenario"]["temperature_is_uniform_downstream_probe_proxy"] is True
    assert audit["measured_spatial_mixing_or_streamer_motion_available"] is False
    assert diffusion == pytest.approx(6.607561922335203e-5, rel=1e-5)
    assert gap == pytest.approx(0.002)
    assert period == pytest.approx(1 / 22_000)
    assert audit["one_dimensional_rms_diffusion_per_ac_period_m"] == pytest.approx(
        math.sqrt(2 * diffusion * period)
    )
    assert audit["one_dimensional_rms_gap_crossing_time_s"] == pytest.approx(
        gap**2 / (2 * diffusion)
    )
    assert audit["one_d_rms_diffusion_less_than_gap_per_ac_period"] is True
    assert audit["one_d_rms_diffusion_less_than_gap_at_16_encounters"] is True
    assert audit["rms_diffusion_by_complete_remix_encounter_count_m"]["16"] == pytest.approx(
        0.0014438989185841956, rel=1e-5
    )
