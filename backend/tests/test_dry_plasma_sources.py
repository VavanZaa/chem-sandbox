import pytest

from scripts.audit_dry_plasma_sources import (
    PAVAN_MECHANISM_SHA256,
    _cantera_pavan_ch3_rate_coefficient,
    _post_solve_ch3_rate_coefficient,
    _published_ch3_curve,
    _published_fixed_electron_rates,
    audit_dry_plasma_sources,
)


def test_published_electron_reference_is_hash_checked_and_parsed() -> None:
    rates = _published_fixed_electron_rates()
    assert rates == {"CH4=>CH3+H": 416.0, "N2=>2N": 13.5}


def test_published_ch3_curve_stays_a_post_solve_diagnostic() -> None:
    curve = _published_ch3_curve()
    assert curve[0][0] == 7.5
    assert curve[1][0] == 0.0
    assert curve[0][-1] == 100.0
    assert curve[1][-1] == pytest.approx(1.2743e-20)
    assert len(curve[0]) == 35
    below_threshold = {
        "eedf": {
            "energy_levels_eV": [0.5, 2.5, 5.0],
            "distribution_eV_minus_3_over_2": [1.0, 1.0, 1.0],
            "normalization": 1.0,
        }
    }
    assert _post_solve_ch3_rate_coefficient(below_threshold, curve) == 0.0
    below_threshold["eedf"]["energy_levels_eV"] = [0.5, 10.0, 101.0]
    with pytest.raises(ValueError, match="valid domain"):
        _post_solve_ch3_rate_coefficient(below_threshold, curve)


def test_packaged_cantera_methane_channel_has_physical_rate_units() -> None:
    coefficient, digest, cantera_version, mean_energy = _cantera_pavan_ch3_rate_coefficient(0.02)
    assert 1e-15 < coefficient < 1e-13  # m^3/molecule/s, not m^3/kmol/s
    assert digest == PAVAN_MECHANISM_SHA256
    assert cantera_version.startswith("3.")
    assert 1 < mean_energy < 20


def test_dry_source_audit_exposes_missing_methane_fragments_without_fitting() -> None:
    audit = audit_dry_plasma_sources()
    assert audit["published_product_yields_used"] is False
    assert audit["reference_conditions"]["MKM_exposure_time_used_for_comparison"] is False
    assert audit["independent_CH3_curve"]["published_full_collision_set_used"] is False
    assert audit["independent_CH3_curve"]["Morgan_EEDF_recomputed_with_CH3_curve"] is False
    assert audit["cantera_pavan_shadow"]["eedf_solved_separately"] is True
    assert len(audit["cases"]) == 4
    for case in audit["cases"]:
        assert case["eedf_converged"] is True
        assert case["N2_to_N_plus_N_source_frequency_s_inv"] > 0
        assert case["CH4_unresolved_dissociation_channel_count"] > 0
        assert case["CH4_to_CH3_plus_H_explicit_channel_count"] == 0
        assert case["published_CH3_curve_post_solve_diagnostic_frequency_s_inv"] > 0
        assert case["cantera_pavan_CH3_plus_H_frequency_s_inv"] > 0
        assert case["cantera_pavan_mean_electron_energy_eV"] > 0
        assert case["unresolved_CH4_sum_is_not_a_CH3_plus_H_rate"] is True
