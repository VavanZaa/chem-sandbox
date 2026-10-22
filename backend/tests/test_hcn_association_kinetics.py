"""Source-form audit of the sensitive H + HCN gas-kinetic step."""

import json
import math

import pytest

from scripts.audit_hcn_association_kinetics import (
    FEED_CARD,
    audit_hcn_association_initial_feed,
)


def test_tsang_estimate_is_read_only_and_keeps_missing_colliders_visible() -> None:
    audit = audit_hcn_association_initial_feed()
    assert audit["source_doi"] == "10.1063/1.555890"
    assert audit["source_rate_type"].startswith("evaluated estimate")
    assert audit["source_temperature_in_range"] is False
    assert audit["source_temperature_extrapolation_below_minimum_K"] == pytest.approx(1.85)
    assert audit["initial_N2_balance_mole_fraction"] == pytest.approx(0.917114)
    assert audit["gri_effective_rate_cm3_molecule_s"] == pytest.approx(1.6314699675e-12)
    assert audit["tsang_estimated_high_pressure_limit_cm3_molecule_s"] == pytest.approx(
        1.5455986383e-14
    )
    assert audit["gri_over_tsang_estimated_high_pressure_limit"] == pytest.approx(
        105.5558621
    )
    reference = audit["source_minimum_temperature_reference"]
    assert reference["temperature_K"] == 300.0
    assert reference["same_hypothetical_initial_composition"] is True
    assert reference["gri_effective_rate_cm3_molecule_s"] == pytest.approx(
        1.6197313816e-12
    )
    assert reference["gri_over_tsang_estimated_high_pressure_limit"] == pytest.approx(
        99.6430005
    )
    assert audit["tsang_reviewed_collider_k0_times_density_sum_cm3_molecule_s"] == pytest.approx(
        4.7336015339e-17
    )
    assert audit["missing_collider_terms"] == ["H2", "H2O", "O2"]
    assert audit["reviewed_collider_terms_are_complete_mixture_rate"] is False
    assert audit["tsang_high_pressure_limit_is_a_measured_bound"] is False
    assert audit["gri_mechanism_modified"] is False
    assert audit["experimental_hcn_assay_loaded"] is False


def test_tsang_formula_is_source_domain_checked_at_300K(tmp_path) -> None:
    feed = json.loads(FEED_CARD.read_text(encoding="utf-8"))
    feed["diagnostic_only_assumptions"]["gas_and_water_temperature_K"] = 300.0
    path = tmp_path / "feed.json"
    path.write_text(json.dumps(feed), encoding="utf-8")
    audit = audit_hcn_association_initial_feed(path)
    assert audit["source_temperature_in_range"] is True
    assert audit["tsang_estimated_high_pressure_limit_cm3_molecule_s"] == pytest.approx(
        5.5e-11 * math.exp(-2438.0 / 300.0)
    )
    assert sum(audit["tsang_reviewed_collider_k0_times_density_cm3_molecule_s"].values()) == pytest.approx(
        audit["tsang_reviewed_collider_k0_times_density_sum_cm3_molecule_s"]
    )


def test_tsang_audit_rejects_invalid_feed_before_rate_comparison(tmp_path) -> None:
    feed = json.loads(FEED_CARD.read_text(encoding="utf-8"))
    feed["measured_initial_gas_mole_fractions"]["H2"] = 2.0
    path = tmp_path / "invalid_feed.json"
    path.write_text(json.dumps(feed), encoding="utf-8")
    with pytest.raises(ValueError, match="no positive N2 balance"):
        audit_hcn_association_initial_feed(path)
