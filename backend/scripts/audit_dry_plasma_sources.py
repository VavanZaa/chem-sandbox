"""Compare blind EEDF collision sources with a dry CH4/N2 plasma reference.

This is a source-coverage diagnostic, not a reactor or product-yield benchmark.
The fixed electron rates in Kamarinopoulou et al. are read only after checking
the preserved source-file hash. No published product yield enters this module.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import cantera as ct
import numpy as np
from scipy.constants import Avogadro, e as elementary_charge, m_e as electron_mass

from app.prebiotic.boltzmann import solve_self_consistent_eedf
from app.prebiotic.schema import BoltzmannEEDFRequest


REFERENCE_DIR = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "reference"
    / "kamarinopoulou_2024"
)
METHANE_FEED_FRACTIONS = (0.001, 0.005, 0.02, 0.07)
ELECTRON_DENSITY_M3 = 1e17  # Published MKM control: 10^11 cm^-3.
PAVAN_MECHANISM_SHA256 = "3789dc895d9ad0bd5437bb84386b9c274ca553e9ca52f98417f3a35d8f04bbcc"
CH3_CURVE_PATH = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "reference"
    / "bouwman_2021"
    / "ch3_neutral_dissociation_table_b2.json"
)


def _published_ch3_curve() -> tuple[np.ndarray, np.ndarray, str]:
    """Read the paper's Table B2, never as another EEDF loss process."""
    raw = CH3_CURVE_PATH.read_bytes()
    record = json.loads(raw)
    energy = np.asarray(record["tabulated_energy_eV"], dtype=float)
    sigma = np.asarray(record["tabulated_cross_section_m2"], dtype=float)
    if (
        energy.shape != sigma.shape
        or energy.size < 2
        or not np.all(np.isfinite(energy))
        or not np.all(np.isfinite(sigma))
        or not np.all(np.diff(energy) > 0)
        or np.any(sigma < 0)
        or record["target"] != "CH4"
        or record["identified_fragment"] != "CH3"
        or energy[0] <= record["reported_threshold_eV"]
    ):
        raise ValueError("The published CH3 diagnostic curve is invalid")
    energy = np.insert(energy, 0, float(record["reported_threshold_eV"]))
    sigma = np.insert(sigma, 0, 0.0)
    return energy, sigma, hashlib.sha256(raw).hexdigest()


def _post_solve_ch3_rate_coefficient(solved: dict, curve: tuple[np.ndarray, np.ndarray, str]) -> float:
    """Fold the published outcome curve through the *unchanged* solved EEDF."""
    tabulated_energy, tabulated_sigma, _ = curve
    eedf = solved["eedf"]
    energy = np.asarray(eedf["energy_levels_eV"], dtype=float)
    distribution = np.asarray(eedf["distribution_eV_minus_3_over_2"], dtype=float)
    if energy.size != distribution.size or energy.size < 2 or energy[-1] > tabulated_energy[-1]:
        raise ValueError("The EEDF grid exceeds the published CH3 table's valid domain")
    normalization = float(eedf["normalization"])
    if normalization <= 0 or not np.isfinite(normalization):
        raise ValueError("The EEDF normalization must be positive and finite")
    probability = np.sqrt(energy) * distribution / normalization
    electron_speed = np.sqrt(2 * elementary_charge * energy / electron_mass)
    sigma = np.interp(energy, tabulated_energy, tabulated_sigma, left=0.0)
    return float(np.trapezoid(sigma * electron_speed * probability, energy))


def _cantera_pavan_ch3_rate_coefficient(methane_fraction: float) -> tuple[float, str, str, float]:
    """Independent packaged EEDF + channel result, not a wet-reactor source."""
    mechanism = next(
        (
            Path(directory) / "example_data" / "methane-plasma-pavan-2023.yaml"
            for directory in ct.get_data_directories()
            if (Path(directory) / "example_data" / "methane-plasma-pavan-2023.yaml").is_file()
        ),
        None,
    )
    if mechanism is None:
        raise ValueError("Cantera's packaged Pavan plasma mechanism is unavailable")
    digest = hashlib.sha256(mechanism.read_bytes()).hexdigest()
    if digest != PAVAN_MECHANISM_SHA256:
        raise ValueError("The installed Cantera/Pavan collision mechanism failed its hash check")
    phase = ct.Solution(str(mechanism), "gri30-plasma", transport_model=None)
    phase.TPX = 550.0, ct.one_atm, {
        "CH4": methane_fraction,
        "N2": 1.0 - methane_fraction,
        "e": 1e-12,
    }
    phase.reduced_electric_field = 300.0e-21
    phase.update_electron_energy_distribution()
    expected = "CH4 + e => CH3 + H + e"
    matches = [
        index for index, reaction in enumerate(phase.reactions())
        if reaction.equation == expected and reaction.reaction_type == "electron-collision-plasma"
    ]
    if len(matches) != 1:
        raise ValueError("Expected one Cantera/Pavan electron-impact CH3 + H channel")
    # Cantera reports bimolecular electron-collision k in m^3/kmol/s.
    # Divide by 1000 * Avogadro for m^3/molecule/s.
    coefficient = float(phase.forward_rate_constants[matches[0]] / (1000 * Avogadro))
    return coefficient, digest, ct.__version__, float(phase.mean_electron_energy)


def _published_fixed_electron_rates() -> dict[str, float]:
    manifest = json.loads((REFERENCE_DIR / "import_manifest.json").read_text(encoding="utf-8"))
    path = REFERENCE_DIR / "kamarinopoulou_2024_fixed_electron.inp"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    expected = manifest["generated_files"][path.name]
    if digest != expected:
        raise ValueError("The published fixed-electron reference file failed its hash check")
    wanted = {"CH4=>CH3+H", "N2=>2N"}
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] in wanted:
            result[parts[0]] = float(parts[1])
    if set(result) != wanted:
        raise ValueError("The published fixed-electron reference reactions are missing")
    return result


def audit_dry_plasma_sources() -> dict:
    """Expose comparable source frequencies and every non-comparability gate."""
    published = _published_fixed_electron_rates()
    ch3_curve = _published_ch3_curve()
    cases = []
    cantera_mechanism_sha256 = ""
    cantera_version = ""
    for methane_fraction in METHANE_FEED_FRACTIONS:
        request = BoltzmannEEDFRequest(
            collision_set="hybrid_morgan_phelps_tabata_2026-09-23",
            mole_fractions={"CH4": methane_fraction, "N2": 1 - methane_fraction},
            gas_temperature_K=550.0,
            reduced_electric_field_Td=300.0,
            energy_max_eV=60.0,
            energy_cells=400,
            acknowledge_two_term_limit=True,
        )
        solved = solve_self_consistent_eedf(request)
        rates = solved["rate_coefficients"]
        unresolved_methane = [
            rate
            for rate in rates
            if rate["target"] == "CH4" and "dissoc" in rate["product_label"].lower()
        ]
        explicit_neutral_methane = [
            rate
            for rate in rates
            if rate["target"] == "CH4"
            and rate["kind"] == "dissociation"
            and rate["product_label"].replace(" ", "") == "CH3+H"
        ]
        nitrogen = [
            rate
            for rate in rates
            if rate["target"] == "N2"
            and rate["product_label"].replace(" ", "") == "N+N"
        ]
        if len(nitrogen) != 1:
            raise ValueError("Expected exactly one reviewed N2 neutral-dissociation outcome")
        n2_frequency_s_inv = nitrogen[0]["rate_coefficient_m3_s"] * ELECTRON_DENSITY_M3
        ch3_rate_coefficient = _post_solve_ch3_rate_coefficient(solved, ch3_curve)
        ch3_frequency_s_inv = ch3_rate_coefficient * ELECTRON_DENSITY_M3
        cantera_coefficient, cantera_mechanism_sha256, cantera_version, cantera_mean_energy = (
            _cantera_pavan_ch3_rate_coefficient(methane_fraction)
        )
        cases.append(
            {
                "initial_CH4_mole_fraction": methane_fraction,
                "eedf_converged": solved["solver"]["converged"],
                "mean_electron_energy_eV": solved["eedf"]["mean_energy_eV"],
                "cantera_pavan_mean_electron_energy_eV": cantera_mean_energy,
                "N2_to_N_plus_N_source_frequency_s_inv": n2_frequency_s_inv,
                "N2_to_N_plus_N_ratio_to_published_fixed_rate": (
                    n2_frequency_s_inv / published["N2=>2N"]
                ),
                "CH4_unresolved_dissociation_channel_count": len(unresolved_methane),
                "CH4_unresolved_dissociation_frequency_sum_s_inv": (
                    sum(rate["rate_coefficient_m3_s"] for rate in unresolved_methane)
                    * ELECTRON_DENSITY_M3
                ),
                "CH4_to_CH3_plus_H_explicit_channel_count": len(explicit_neutral_methane),
                "published_CH3_curve_post_solve_diagnostic_frequency_s_inv": ch3_frequency_s_inv,
                "cantera_pavan_CH3_plus_H_frequency_s_inv": (
                    cantera_coefficient * ELECTRON_DENSITY_M3
                ),
                "published_CH3_curve_ratio_to_published_fixed_rate": (
                    ch3_frequency_s_inv / published["CH4=>CH3+H"]
                ),
                "published_CH4_to_CH3_plus_H_frequency_s_inv": published["CH4=>CH3+H"],
                "unresolved_CH4_sum_is_not_a_CH3_plus_H_rate": True,
            }
        )
    return {
        "status": "dry_plasma_source_coverage_diagnostic_not_product_validation",
        "source": "Kamarinopoulou et al., Science Advances 10 (2024), adl4246",
        "doi": "10.1126/sciadv.adl4246",
        "reference_conditions": {
            "reduced_electric_field_Td": 300.0,
            "gas_temperature_K": 550.0,
            "electron_density_m3": ELECTRON_DENSITY_M3,
            "MKM_exposure_time_used_for_comparison": False,
        },
        "published_product_yields_used": False,
        "measured_discharge_power_used_to_fit_sources": False,
        "source_file_hash_checked": True,
        "independent_CH3_curve": {
            "source": "Bouwman et al. (2021), Table B2; https://ir.cwi.nl/pub/30824/",
            "transcription_sha256": ch3_curve[2],
            "published_full_collision_set_used": False,
            "Morgan_EEDF_recomputed_with_CH3_curve": False,
            "result_scope": "post-solve fragment-rate sensitivity, not production chemistry or yield",
        },
        "cantera_pavan_shadow": {
            "mechanism": "Cantera example_data/methane-plasma-pavan-2023.yaml",
            "mechanism_sha256": cantera_mechanism_sha256,
            "cantera_version": cantera_version,
            "eedf_solved_separately": True,
            "result_scope": "packaged dry plasma source sensitivity, not independent experimental evidence or a wet yield",
        },
        "cases": cases,
        "limitations": [
            "Different complete electron-collision sets and two-term solvers need not agree at identical E/N.",
            "The published N2 fixed electron rate is a cross-set diagnostic, not an experimental rate measurement.",
            "Morgan CH4 dissociation labels do not resolve CH3 + H; their sum cannot be compared to that channel.",
            "The CH3 diagnostic folds a published outcome curve through a Morgan-based EEDF; it is not a consistent whole-set replacement or a validated rate.",
            "Cantera/Pavan uses its own collision/EEDF data but may share upstream inputs with published compilations; agreement is not independent validation.",
            "The published plasma exposure time was inferred by fitting reactor data and is not used here.",
            "A dry DBD source comparison is not a wet recirculating Miller-Urey yield prediction.",
        ],
    }


def main() -> None:
    print(json.dumps(audit_dry_plasma_sources(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
