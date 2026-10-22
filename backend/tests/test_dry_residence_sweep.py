"""The Figure 1D flow trend is an external check, not a source setting."""

import pytest

from evaluation.reference_only.dry_hcn_2024 import evaluate_frozen_dry_residence_sweep
from scripts.audit_dry_residence_sweep import (
    NOMINAL_RESIDENCE_TIMES_S,
    audit_dry_residence_sweep,
)
from scripts.run_dry_flow_hcn_probe import run_dry_flow_probe


def test_residence_sweep_closes_ledgers_but_misses_plot_magnitude() -> None:
    sweep = audit_dry_residence_sweep()
    assert sweep["measured_product_assay_loaded"] is False
    assert sweep["reduced_electric_field_Td_scenario"] == 184.0
    assert sweep["residence_convention_matches_paper_exactly"] is False
    assert sweep["microdischarge_exposure_independently_constrained"] is False
    rows = sweep["rows"]
    assert [row["nominal_residence_s"] for row in rows] == list(
        NOMINAL_RESIDENCE_TIMES_S
    )
    assert [row["nominal_inferred_flow_sccm"] for row in rows] == pytest.approx(
        [50.10638298, 100.21276596, 149.52380952, 200.42553191], rel=1e-7
    )
    assert [row["ch4_conversion_fraction"] for row in rows] == pytest.approx(
        [0.5185616479, 0.3065714771, 0.2163728601, 0.1652167156], rel=1e-5
    )
    assert all(row["hcn_mol_below_1e-12_reporting_floor"] for row in rows)
    assert all(row["material_ledger_closed"] for row in rows)
    assert all(row["collision_energy_upper_bound_closed"] for row in rows)
    assert all(row["n2_product_state_source_excluded"] for row in rows)
    report = evaluate_frozen_dry_residence_sweep(sweep)
    assert report["nominal_residence_points_match_reference"] is True
    assert report["conversion_directional_trend_matches"] is True
    assert report["long_residence_coarse_plot_bound_matches"] is False
    assert report["short_residence_coarse_plot_bound_matches"] is False
    assert report["coarse_plot_bounds_are_statistical_uncertainty"] is False
    assert report["like_for_like_quantitative_comparison_ready"] is False
    assert report["reference_used_by_solver"] is False
    assert report["numeric_experimental_conversions_digitized"] is False
    assert len(report["frozen_prediction_sha256"]) == 64


def test_residence_control_rejects_out_of_scope_values() -> None:
    with pytest.raises(ValueError, match="Nominal residence time"):
        run_dry_flow_probe(nominal_residence_s=0.0)
