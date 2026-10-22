"""External electron-swarm residuals stay separate from HCN calibration."""

import pytest

from scripts.audit_dry_n2_swarm_transport import audit_dry_n2_swarm_transport


@pytest.fixture(scope="module")
def audit() -> dict:
    return audit_dry_n2_swarm_transport()


def test_nitrogen_survey_residual_is_explicit_not_a_calibration(audit: dict) -> None:
    rows = audit["pure_n2_295K_comparison"]
    assert len(rows) == audit["reference"]["selected_point_count"] == 9
    assert all(row["eedf_converged"] for row in rows)
    assert all(row["measured_drift_m_s"] > 0 for row in rows)
    assert all(
        row["model_minus_measured_fraction"]
        == pytest.approx(row["model_two_term_drift_m_s"] / row["measured_drift_m_s"] - 1)
        for row in rows
    )
    at_450 = next(row for row in rows if row["field_Td"] == 450.0)
    assert at_450["measured_drift_m_s"] == pytest.approx(450000.0)
    assert -0.28 < at_450["model_minus_measured_fraction"] < -0.20
    assert audit["reference"]["phelps_swarm_fit_training_overlap_ruled_out"] is False
    assert audit["measured_drift_used_as_solver_input"] is False
    assert audit["field_selected_from_hcn_or_products"] is False
    assert audit["hcn_assay_loaded"] is False
    assert audit["revised_dry_reactor_exposure_from_survey"] is False


def test_high_field_residual_is_not_eedf_grid_artifact(audit: dict) -> None:
    grid = audit["pure_n2_450Td_grid_check"]
    assert all(row["eedf_converged"] for row in grid)
    velocities = [row["model_two_term_drift_m_s"] for row in grid]
    assert max(velocities) / min(velocities) < 1.005
    assert grid[-1]["model_two_term_drift_m_s"] == pytest.approx(
        next(row for row in audit["pure_n2_295K_comparison"] if row["field_Td"] == 450.0)[
            "model_two_term_drift_m_s"
        ]
    )
