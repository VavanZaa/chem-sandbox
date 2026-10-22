"""Guard the private N2+ record against extrapolation or false source terms."""

import hashlib

import pytest

from app.prebiotic.n2plus_recombination import (
    PRIVATE_RECORD,
    REVIEWED_SHA256,
    audit_n2plus_recombination,
    load_private_n2plus_recombination,
)


def test_private_n2plus_record_is_hash_pinned_and_low_energy_only() -> None:
    assert hashlib.sha256(PRIVATE_RECORD.read_bytes()).hexdigest() == REVIEWED_SHA256
    energy, sigma_low, sigma_high = load_private_n2plus_recombination()
    assert len(energy) == 31899
    assert energy[0] == pytest.approx(1e-5)
    assert energy[-1] == pytest.approx(0.39999)
    assert sigma_low[0] > sigma_low[-1] > 0
    assert (sigma_high >= sigma_low).all()
    assert (sigma_high > sigma_low).sum() > 0


def test_n2plus_audit_reports_only_supported_window_without_product_source() -> None:
    eedf = {
        "energy_levels_eV": [0.1, 0.2, 0.4, 1.0, 2.0],
        "probability_density_eV_inv": [0.2, 0.3, 0.4, 0.2, 0.1],
    }
    result = audit_n2plus_recombination(eedf)
    assert eedf["probability_density_eV_inv"] == [0.2, 0.3, 0.4, 0.2, 0.1]
    assert result["common_energy_window_eV"] == pytest.approx([0.1, 0.39999])
    assert result["unresolved_low_energy_window_eV"] == pytest.approx([1e-5, 0.1])
    assert result["record_ends_below_eedf_maximum"] is True
    assert 0 < result["eedf_probability_in_common_window"] < 1
    low, high = result["partial_rate_coefficient_envelope_m3_s"]
    assert 0 < low <= high
    assert result["rounded_duplicate_energy_rows"] == 8100
    assert result["total_rate_coefficient_m3_s"] is None
    assert result["neutral_product_states_resolved"] is False
    assert result["allowed_as_chemical_source"] is False


def test_n2plus_audit_rejects_invalid_eedf() -> None:
    with pytest.raises(ValueError, match="increasing tabulated EEDF"):
        audit_n2plus_recombination({
            "energy_levels_eV": [0.1, 0.1],
            "probability_density_eV_inv": [1.0, 1.0],
        })
