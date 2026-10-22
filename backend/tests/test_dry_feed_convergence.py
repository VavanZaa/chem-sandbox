"""Feed-trend discrepancy must not be blamed on chemistry substep size."""

import pytest

from scripts.audit_dry_feed_convergence import audit_dry_feed_convergence
from scripts.audit_dry_feed_sweep import audit_dry_feed_sweep


def test_feed_trend_and_selectivities_stable_as_substeps_refine() -> None:
    audit = audit_dry_feed_convergence()
    assert audit["substep_counts"] == (16, 32, 64)
    assert audit["measured_product_assay_loaded"] is False
    assert audit["streamer_pulses_resolved"] is False
    assert audit["all_material_ledgers_closed"] is True
    assert audit["all_collision_energy_upper_bounds_closed"] is True
    assert audit["all_n2_product_state_sources_excluded"] is True
    assert all(delta > 0 for delta in audit["ch4_conversion_feed_trend_by_substeps"].values())
    for row in audit["rows"]:
        assert row["last_pair_ch4_conversion_absolute_change"] < 1e-5
        assert row["last_pair_ethane_selectivity_absolute_change"] < 2e-6
        assert row["last_pair_hydrogen_selectivity_absolute_change"] < 1e-6
    assert audit["ch4_conversion_feed_trend_by_substeps"][64] == pytest.approx(
        0.0038801791, rel=1e-4
    )


@pytest.mark.parametrize("substeps", [0, 129, 1.5, "16"])
def test_feed_sweep_rejects_invalid_substeps(substeps: object) -> None:
    with pytest.raises(ValueError, match="substeps"):
        audit_dry_feed_sweep(substeps=substeps)


@pytest.mark.parametrize("counts", [(16,), (32, 16), (16, 16), (0, 32)])
def test_convergence_audit_rejects_invalid_grid(counts: tuple[int, ...]) -> None:
    with pytest.raises(ValueError, match="Substep counts"):
        audit_dry_feed_convergence(substep_counts=counts)
