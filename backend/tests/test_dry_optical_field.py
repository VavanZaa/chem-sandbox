"""Independent optical ratio cannot silently calibrate the dry field."""

import pytest

from app.prebiotic.boltzmann import _load_phelps_n2_processes
from scripts.audit_dry_optical_field import audit_dry_optical_field


def test_installed_optical_channel_spot_values_match_phelps_public_note() -> None:
    # Phelps' revised ELECTRON.TXT note (2008), pages 16 and 18; four
    # independently published spot values, not a replacement collision set.
    # The note uses 1e-16 cm2 = 1e-20 m2 units for its tabulated values.
    expected = {
        "N2(C3)": {20.0: 0.2500e-20, 30.0: 0.1110e-20},
        "N2^+(B2SIGMA)": {30.0: 0.0720e-20, 50.0: 0.1470e-20},
    }
    processes = {
        process.get("product"): dict(process["data"])
        for process in _load_phelps_n2_processes()
        if process.get("product") in expected
    }
    assert set(processes) == set(expected)
    for product, samples in expected.items():
        for energy_eV, cross_section_m2 in samples.items():
            assert processes[product][energy_eV] == pytest.approx(cross_section_m2)


def test_phelps_optical_ratio_misses_paper_diagnostic_without_hcn_input() -> None:
    audit = audit_dry_optical_field()
    assert audit["paper_supplement_equation"] == "S24"
    assert audit["installed_cantera_phelps_export_date"] == "2025-08-01"
    assert audit["paper_phelps_lxcat_retrieval_date"] == "2023-03-01"
    assert len(audit["installed_cantera_phelps_sha256"]) == 64
    assert audit["hcn_or_product_assay_used"] is False
    assert audit["electric_field_selected_by_this_audit"] is False
    assert audit["photodetector_response_and_full_emission_transport_validated"] is False
    assert audit["paper_optical_run_not_identical_to_figure1_nominal_power"] is True
    assert audit["paper_measured_intensity_ratio"] == pytest.approx(3.47e-3)
    assert audit["n2_plus_b_radiative_quench_factor"] < audit["n2_c_radiative_quench_factor"]
    fields = audit["field_grid"]
    assert [row["reduced_electric_field_Td"] for row in fields] == [
        100.0, 189.0, 300.0, 435.0, 500.0
    ]
    assert all(row["eedf_converged"] for row in fields)
    assert all(
        not any("clipped a negative reconstructed elastic" in warning for warning in row["eedf_numerical_warnings"])
        for row in fields
    )
    assert all(
        a["predicted_spectroscopic_intensity_ratio"]
        < b["predicted_spectroscopic_intensity_ratio"]
        for a, b in zip(fields, fields[1:])
    )
    assert fields[3]["predicted_spectroscopic_intensity_ratio"] == pytest.approx(
        5.9691e-4, rel=1e-3
    )
    assert fields[3]["same_cross_section_cantera_reference_available"] is True
    assert fields[3]["same_cross_section_cantera_predicted_intensity_ratio"] == pytest.approx(
        4.64226e-4, rel=1e-3
    )
    assert abs(fields[3]["same_cross_section_cantera_relative_mean_energy_difference"]) < 0.05
    assert fields[3]["same_cross_section_cantera_predicted_intensity_ratio"] < audit[
        "paper_measured_intensity_ratio"
    ] / 7
    assert fields[-1]["ratio_to_paper_measured_intensity_ratio"] < 0.26
    grid = audit["energy_grid_sensitivity_at_435_Td"]
    assert all(row["eedf_converged"] for row in grid)
    ratios = [row["predicted_spectroscopic_intensity_ratio"] for row in grid]
    assert max(ratios) / min(ratios) < 1.01
    pure_n2 = audit["pure_n2_eedf_sensitivity_at_435_Td"]
    assert pure_n2["eedf_mole_fractions"] == {"N2": 1.0}
    assert pure_n2["eedf_converged"] is True
    assert pure_n2["predicted_spectroscopic_intensity_ratio"] == pytest.approx(
        6.12021e-4, rel=1e-3
    )
    assert pure_n2["ratio_to_paper_measured_intensity_ratio"] < 0.18
    room_temperature = audit["room_temperature_sensitivity_at_435_Td"]
    assert room_temperature["gas_temperature_K"] == pytest.approx(298.15)
    assert room_temperature["eedf_converged"] is True
    assert room_temperature["predicted_spectroscopic_intensity_ratio"] == pytest.approx(
        5.72503e-4, rel=1e-3
    )
    assert room_temperature["ratio_to_paper_measured_intensity_ratio"] < 0.17
