from __future__ import annotations

from collections import Counter
from functools import lru_cache
import hashlib
import io
from importlib.metadata import PackageNotFoundError, version
import json
import logging
from pathlib import Path
import re
from typing import Any
import warnings

import cantera as ct
import numpy as np
import scipy.integrate as scipy_integrate
import yaml
from scipy.constants import Boltzmann, elementary_charge
from scipy.constants import atomic_mass, electron_mass

from .collision_intake import (
    normalize_lxcat_bolsig_text,
    normalize_lxcat_processes,
    split_lxcat_products,
)
from .schema import BoltzmannEEDFRequest


# BOLOS 0.2 predates SciPy's removal of the ``simps`` alias. The numerical
# implementation is unchanged; expose the new function under the legacy name
# until BOLOS publishes a compatibility release.
if not hasattr(scipy_integrate, "simps"):
    scipy_integrate.simps = scipy_integrate.simpson  # type: ignore[attr-defined]


DATASET_ID = "lxcat_morgan_2026-09-22"
DATASET_SHA256 = "86e256c3ec690cb3ff759529f099959df777f4d0a331f769925e5265c9e8efd6"
HYBRID_DATASET_ID = "hybrid_morgan_phelps_2026-09-23"
PHELPS_N2_SHA256 = "be36740e0db01d996668b2c3ff6cbf12a0f296155086a16b7335d26df43f3ce2"
HYBRID_DATASET_SHA256 = "c099abe09c174bd0f789bd840dca4861799b085d49b38445201eb24f0188d84e"
REACTIVE_N2_DATASET_ID = "hybrid_morgan_phelps_tabata_2026-09-23"
TABATA_N2_SHA256 = "cba384a025d4003deecb7fe0ae83ad8e6b5836c0608b6084f5861972695513d4"
REACTIVE_N2_DATASET_SHA256 = "ef8213da28c9ea0db3b8969a4d95ec5dddd4fc1d4042cb03c7f93df6518e9ae1"
COMMUNITY_CH4_DATASET_ID = "hybrid_morgan_community_phelps_tabata_2026-09-23"
COMMUNITY_CH4_SHA256 = "e890473214d4373deeb6735df460e2ef9558d4f3008fbc312029fa42f606ed03"
COMMUNITY_CH4_DATASET_SHA256 = "2d70dfc503c18ea735618b7e3b1eb97984bc52442c2feaed6de0565546248cc0"
IST_H2_COMMUNITY_CH4_DATASET_ID = "hybrid_ist_h2_community_ch4_phelps_tabata_2026-09-23"
IST_H2_SHA256 = "2721e8cb5b0683b614c9094c279c104f38e987fcb29455fb168cb39e0e511468"
IST_H2_COMMUNITY_CH4_DATASET_SHA256 = "b18769ffdd0b4d827ab7c376d7332dcc3863bd301df0a5cdc13048cd5fac83e5"
IST_CO2_O2_DATASET_ID = "hybrid_ist_h2_co2_o2_phelps_tabata_2026-09-26"
IST_CO2_O2_SHA256 = "521e279b77a7f030935bfd4ff3676d1b4d2cfe9654fd6479591dc4c8b2b884c7"
DATASET_TARGETS = {
    DATASET_ID: ("H2", "H2O", "CH4", "NH3"),
    HYBRID_DATASET_ID: ("H2", "H2O", "CH4", "NH3", "N2"),
    REACTIVE_N2_DATASET_ID: ("H2", "H2O", "CH4", "NH3", "N2"),
    COMMUNITY_CH4_DATASET_ID: ("H2", "H2O", "CH4", "NH3", "N2"),
    IST_H2_COMMUNITY_CH4_DATASET_ID: ("H2", "H2O", "CH4", "NH3", "N2"),
    IST_CO2_O2_DATASET_ID: ("H2", "H2O", "CH4", "NH3", "N2", "CO2", "O2"),
}
DATASET_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "collision"
    / "private"
    / "Morgan_complete_H2_H2O_CH4_NH3_2026-09-22"
    / "Cross section.txt"
)
TABATA_N2_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "collision"
    / "qst_tabata_n2_neutral_dissociation_2006.json"
)
COMMUNITY_CH4_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "collision"
    / "private"
    / "Community_Song_Bouwman_CH4_2026-09-23"
    / "Cross section.txt"
)
IST_H2_PATH = (
    Path(__file__).resolve().parents[2]
    / "data" / "collision" / "private" / "IST_Lisbon_H2_2026-09-23"
    / "Cross section.txt"
)
IST_CO2_O2_PATH = (
    Path(__file__).resolve().parents[2]
    / "data" / "collision" / "private" / "IST_Lisbon_CO2_O2_2026-09-26"
    / "Cross section.txt"
)
SOURCE_REFERENCE = (
    "Morgan database, www.lxcat.net, retrieved on September 22, 2026"
)
PHELPS_SOURCE_REFERENCE = (
    "A. V. Phelps electron cross-section compilation distributed with Cantera 3.2; "
    "original contributor data available through LXCat"
)
TABATA_SOURCE_REFERENCE = (
    "T. Tabata, T. Shirai, M. Sataka, and H. Kubo, Atomic Data and Nuclear Data "
    "Tables 92, 375-406 (2006), DOI 10.1016/j.adt.2006.02.002; official QST "
    "total N + N production table recommended by Cosby; electronic product states unresolved"
)
COMMUNITY_CH4_SOURCE_REFERENCE = (
    "Community database, www.lxcat.net, retrieved on September 23, 2026; "
    "D. Bouwman et al., Plasma Sources Science and Technology 30, 075012 (2021), "
    "DOI 10.1088/1361-6595/ac0b2b; M.-Y. Song et al., Journal of Physical "
    "and Chemical Reference Data 44, 023101 (2015), DOI 10.1063/1.4918630"
)
IST_H2_SOURCE_REFERENCE = (
    "IST-Lisbon database, www.lxcat.net, retrieved September 23, 2026; "
    "S. Baghel, L. Marques, and L. L. Alves, Plasma Sources Science and "
    "Technology 35, 085023 (2026), DOI 10.1088/1361-6595/ae918a; "
    "contributor validation states 1-200 Td without additional rotational data"
)
IST_CO2_O2_SOURCE_REFERENCE = (
    "IST-Lisbon database, www.lxcat.net, retrieved September 26, 2026; "
    "M. Grofulovic, L. L. Alves and V. Guerra, J. Phys. D 49, 395207 (2016); "
    "Y. Liu et al., Plasma Sources Sci. Technol. 34, 035003 (2025); "
    "G. Gousset et al., J. Phys. D 24, 290 (1991); "
    "L. L. Alves et al., Eur. Phys. J. D 70, 124 (2016)"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_morgan_processes() -> list[dict[str, Any]]:
    if not DATASET_PATH.is_file():
        raise ValueError(
            "The private LXCat/Morgan collision set is not installed. "
            "Download it through LXCat's complete-set workflow before solving the EEDF."
        )
    actual = _sha256(DATASET_PATH)
    if actual != DATASET_SHA256:
        raise ValueError(
            "The installed collision set does not match its reviewed SHA-256 manifest; "
            f"expected {DATASET_SHA256}, received {actual}."
        )
    try:
        from bolos import parser
    except ImportError as error:
        raise ValueError("BOLOS 0.2 is required for the BOLSIG+ collision-set parser") from error
    with DATASET_PATH.open("r", encoding="utf-8") as stream:
        return parser.parse(stream)


def _load_community_ch4_processes() -> list[dict[str, Any]]:
    """Load the complete, private Song-Bouwman methane set as one replacement."""

    if not COMMUNITY_CH4_PATH.is_file():
        raise ValueError(
            "The private LXCat Community/Song-Bouwman CH4 set is not installed. "
            "Obtain it through LXCat's complete-set workflow and retain it locally."
        )
    actual = _sha256(COMMUNITY_CH4_PATH)
    if actual != COMMUNITY_CH4_SHA256:
        raise ValueError(
            "The private LXCat Community CH4 file does not match its reviewed SHA-256; "
            f"expected {COMMUNITY_CH4_SHA256}, received {actual}."
        )
    from bolos import parser

    with COMMUNITY_CH4_PATH.open("r", encoding="utf-8") as stream:
        parsed = parser.parse(stream)
    processes, changes = normalize_lxcat_processes(parsed, "CH4")
    selected = [process for process in processes if process.get("target") == "CH4"]
    counts = Counter(process["kind"] for process in selected)
    expected = Counter(
        {"ATTACHMENT": 2, "ELASTIC": 1, "EXCITATION": 6, "IONIZATION": 7}
    )
    if counts != expected or len(changes) != 15:
        raise ValueError(
            "The private Community CH4 set does not match the reviewed 16-process "
            "structure and normalization log"
        )
    return selected


def _load_ist_h2_processes() -> list[dict[str, Any]]:
    """Load the complete private H2 set, preserving the raw contributor file."""

    if not IST_H2_PATH.is_file():
        raise ValueError(
            "The private IST-Lisbon H2 collision set is not installed; obtain "
            "the complete contributor file from LXCat and retain it locally."
        )
    actual = _sha256(IST_H2_PATH)
    if actual != IST_H2_SHA256:
        raise ValueError(
            "The private IST-Lisbon H2 file does not match its reviewed SHA-256; "
            f"expected {IST_H2_SHA256}, received {actual}."
        )
    from bolos import parser

    original = IST_H2_PATH.read_text(encoding="utf-8")
    normalized_text, text_changes = normalize_lxcat_bolsig_text(original, "H2")
    parsed = parser.parse(io.StringIO(normalized_text))
    processes, process_changes = normalize_lxcat_processes(parsed, "H2")
    selected = [process for process in processes if process.get("target") == "H2"]
    expected = Counter({"ELASTIC": 1, "EXCITATION": 19, "IONIZATION": 2})
    if (
        len(selected) != 22
        or Counter(process["kind"] for process in selected) != expected
        or len(text_changes) != 1
        or len(process_changes) != 2
    ):
        raise ValueError(
            "The private IST-Lisbon H2 set differs from its reviewed 22-process "
            "structure or three in-memory label normalizations"
        )
    return selected


def _load_ist_co2_o2_processes() -> list[dict[str, Any]]:
    """Load only the reviewed complete CO2 and O2 sets from a private export."""

    if not IST_CO2_O2_PATH.is_file():
        raise ValueError(
            "The private IST-Lisbon CO2/O2 collision file is not installed; "
            "obtain the two complete contributor sets from LXCat."
        )
    actual = _sha256(IST_CO2_O2_PATH)
    if actual != IST_CO2_O2_SHA256:
        raise ValueError(
            "The private IST-Lisbon CO2/O2 file does not match its reviewed SHA-256; "
            f"expected {IST_CO2_O2_SHA256}, received {actual}."
        )
    from bolos import parser

    with IST_CO2_O2_PATH.open("r", encoding="utf-8") as stream:
        processes = parser.parse(stream)
    processes, oxygen_changes = normalize_lxcat_processes(processes, "O2")
    processes, carbon_changes = normalize_lxcat_processes(processes, "CO2")
    by_target = {
        target: Counter(item["kind"] for item in processes if item["target"] == target)
        for target in ("O2", "CO2")
    }
    if (
        len(processes) != 31
        or set(item["target"] for item in processes) != {"O2", "CO2"}
        or by_target["O2"] != Counter({
            "ATTACHMENT": 1, "EFFECTIVE": 1, "EXCITATION": 11, "IONIZATION": 1,
        })
        or by_target["CO2"] != Counter({
            "ATTACHMENT": 1, "ELASTIC": 1, "EXCITATION": 14, "IONIZATION": 1,
        })
        or Counter(change["kind"] for change in oxygen_changes) != {"product_label": 2}
        or carbon_changes
    ):
        raise ValueError(
            "The private IST-Lisbon CO2/O2 export differs from the reviewed "
            "14-channel O2 and 17-channel CO2 complete sets."
        )
    # BOLOS' BOLSIG ATTACHMENT parser defaults the process threshold to zero,
    # although this hash-pinned CO2 record explicitly declares E = 3.9 eV in
    # PARAM. Preserve that physical onset in the in-memory process so the
    # collision-energy lower-bound report cannot silently treat it as free.
    co2_attachment = next(
        item for item in processes
        if item["target"] == "CO2" and item["kind"] == "ATTACHMENT"
    )
    if (
        co2_attachment.get("product") != "CO(X) + O(-,X)"
        or float(co2_attachment.get("threshold", 0.0)) != 0.0
        or float(co2_attachment["data"][0][0]) != 3.9
        or float(co2_attachment["data"][0][1]) != 0.0
        or "PARAM.:  E = 3.9 eV, complete set" not in IST_CO2_O2_PATH.read_text(encoding="utf-8")
    ):
        raise ValueError("CO2 attachment does not match its reviewed 3.9 eV onset")
    co2_attachment["threshold"] = 3.9
    return processes


def _phelps_path() -> Path:
    path = next(
        (
            Path(directory) / "example_data" / "air-plasma-Phelps.yaml"
            for directory in ct.get_data_directories()
            if (Path(directory) / "example_data" / "air-plasma-Phelps.yaml").is_file()
        ),
        None,
    )
    if path is None:
        raise ValueError("Cantera's packaged Phelps N2 collision set is unavailable")
    return path


def _load_phelps_n2_processes() -> list[dict[str, Any]]:
    path = _phelps_path()
    actual = _sha256(path)
    if actual != PHELPS_N2_SHA256:
        raise ValueError(
            "The installed Cantera Phelps file does not match the reviewed SHA-256; "
            f"expected {PHELPS_N2_SHA256}, received {actual}."
        )
    try:
        import yaml
    except ImportError as error:
        raise ValueError("PyYAML is required to import Cantera's N2 collision table") from error
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    mass_ratio = electron_mass / (28.0134 * atomic_mass)
    processes = []
    for item in document.get("electron-collisions", []):
        if item.get("target") != "N2":
            continue
        process = {
            "target": "N2",
            "kind": str(item["kind"]).upper(),
            "data": list(
                zip(
                    (float(value) for value in item["energy-levels"]),
                    (float(value) for value in item["cross-sections"]),
                    strict=True,
                )
            ),
            "comment": PHELPS_SOURCE_REFERENCE,
        }
        if process["kind"] in {"ELASTIC", "EFFECTIVE"}:
            process["mass_ratio"] = mass_ratio
        else:
            process["threshold"] = float(item.get("threshold", item["energy-levels"][0]))
            process["product"] = item.get("product", "N2*")
        processes.append(process)
    if not processes:
        raise ValueError("Cantera's packaged Phelps file contains no N2 collision processes")
    return processes


def _load_tabata_n2_dissociation_process() -> dict[str, Any]:
    """Load the reviewed neutral N2 dissociation outcome cross section.

    The channel is folded through the converged EEDF after the Phelps energy-loss
    solve. It is deliberately marked outcome-only because the Phelps electronic
    manifold already carries the associated energy loss; loading it into BOLOS
    again would double count that loss. The published total N + N cross
    section does not supply state-selective branching. Cosby's measured
    primary pathway is N(2D) + N(4S), so N + N must not be interpreted as
    two ground-state atoms in a ground-state chemical mechanism.
    """

    if not TABATA_N2_PATH.is_file():
        raise ValueError("The reviewed QST/Tabata N2 dissociation table is unavailable")
    actual = _sha256(TABATA_N2_PATH)
    if actual != TABATA_N2_SHA256:
        raise ValueError(
            "The installed QST/Tabata N2 dissociation table does not match its reviewed "
            f"SHA-256; expected {TABATA_N2_SHA256}, received {actual}."
        )
    document = json.loads(TABATA_N2_PATH.read_text(encoding="utf-8"))
    energies = [float(value) for value in document["tabulated_energy_eV"]]
    cross_sections = [
        float(value) * 1e-4 for value in document["tabulated_cross_section_cm2"]
    ]
    if len(energies) != len(cross_sections) or len(energies) < 2:
        raise ValueError("The QST/Tabata N2 dissociation table is malformed")
    if any(right <= left for left, right in zip(energies, energies[1:])):
        raise ValueError("The QST/Tabata N2 dissociation energy grid must increase")
    if any(value < 0 for value in cross_sections):
        raise ValueError("The QST/Tabata N2 dissociation cross sections must be nonnegative")
    return {
        "target": "N2",
        "kind": "DISSOCIATION",
        "threshold": float(document["threshold_eV_from_analytic_fit"]),
        "product": "N + N",
        "data": list(zip(energies, cross_sections, strict=True)),
        "comment": TABATA_SOURCE_REFERENCE,
        "source_sha256": actual,
        "outcome_only": True,
        "product_state_unresolved": True,
        "tabata_analytic_fit_diagnostic_available": True,
    }


def _integrate_rate_window(
    energies_eV: np.ndarray, integrand: np.ndarray, lower_eV: float, upper_eV: float
) -> float:
    """Integrate an EEDF-folded collision rate over an exact energy window.

    Endpoints are interpolated so a coarse energy grid does not move the
    literature-window boundary to the nearest computational cell.
    """

    start = max(float(energies_eV[0]), lower_eV)
    stop = min(float(energies_eV[-1]), upper_eV)
    if stop <= start:
        return 0.0
    interior = (energies_eV > start) & (energies_eV < stop)
    window_energies = np.concatenate(([start], energies_eV[interior], [stop]))
    window_integrand = np.interp(window_energies, energies_eV, integrand)
    return float(np.trapezoid(window_integrand, window_energies))


def _tabata_n2_analytic_fit_m2(energies_eV: np.ndarray) -> np.ndarray:
    """QST's published F3 N+N fit; diagnostic only, not the imported table.

    Source: https://www.qst.go.jp/uploaded/attachment/8709.f . The original
    Fortran uses keV for energy and cm2 for sigma, and explicitly limits
    the fit to 12-200 eV. Below 12 eV the fit is zero even though its
    analytic threshold parameter is 9.76 eV.
    """

    energy = np.asarray(energies_eV, dtype=float)
    excess_keV = np.maximum((energy - 9.76) * 1e-3, 0.0)
    numerator = 1e-20 * 1.870 * (excess_keV / 0.01361) ** 3.030
    denominator = (
        1.0
        + (excess_keV / 0.01230) ** (3.030 + 0.0470)
        + (excess_keV / 0.03500) ** (3.030 + 1.100)
    )
    return np.where((energy >= 12.0) & (energy <= 200.0), numerator / denominator, 0.0)


def _cosby_product_state_energy_support(
    energies_eV: np.ndarray, rate_integrand: np.ndarray, total_rate: float,
) -> dict[str, Any]:
    """Partition a total N2 rate by source energy, never by product state."""

    lower_eV, upper_eV = 18.5, 148.5
    if total_rate > 0:
        fractions = {
            "below": _integrate_rate_window(
                energies_eV, rate_integrand, float(energies_eV[0]), lower_eV
            ) / total_rate,
            "within": _integrate_rate_window(
                energies_eV, rate_integrand, lower_eV, upper_eV
            ) / total_rate,
            "above": _integrate_rate_window(
                energies_eV, rate_integrand, upper_eV, float(energies_eV[-1])
            ) / total_rate,
        }
    else:
        fractions = {"below": None, "within": None, "above": None}
    return {
        "reference": "Cosby, J. Chem. Phys. 98, 9544 (1993), DOI 10.1063/1.464385",
        "measurement_energy_window_eV": [lower_eV, upper_eV],
        "computed_grid_max_energy_eV": float(energies_eV[-1]),
        "measurement_window_extends_above_computed_grid": bool(
            energies_eV[-1] < upper_eV
        ),
        "fraction_of_computed_rate_below_measured_window": fractions["below"],
        "fraction_of_computed_rate_within_measured_window": fractions["within"],
        "fraction_of_computed_rate_above_measured_window": fractions["above"],
        "state_branching_fraction_inferred": False,
    }


def _n2_ground_atom_energetic_window(
    energies_eV: np.ndarray, rate_integrand: np.ndarray, total_rate: float,
) -> dict[str, Any]:
    """Conditional ground-N necessity below the first no-ground channel.

    Ground-state N2 cannot form N(2D)+N(2D) below 14.526 eV. Lower neutral
    N+N channels each have at least one N(4S). This is an energy-access
    constraint on a *total* cross section, not a measured state branch or an
    atom source: Tabata's near-threshold values are below Cosby's measured
    electron-energy window, and the second fragment remains unresolved.
    """

    threshold_eV = 14.526
    fraction = (
        _integrate_rate_window(
            energies_eV, rate_integrand, float(energies_eV[0]), threshold_eV
        ) / total_rate
        if total_rate > 0 else None
    )
    return {
        "first_neutral_channel_without_ground_N4S_threshold_eV": threshold_eV,
        "rate_fraction_below_that_threshold": fraction,
        "conditional_minimum_ground_N4S_atoms_per_event_in_window": 1,
        "condition": (
            "Ground-state N2 target; neutral N+N total cross section is valid "
            "below 14.526 eV; incident electron energy bounds available "
            "fragment electronic plus translational energy."
        ),
        "threshold_source": (
            "N2 dissociation limits, N(4S)+N(2D) 12.143 eV, "
            "N(4S)+N(2P) 13.335 eV, N(2D)+N(2D) 14.526 eV; "
            "DOI 10.1039/D3CP04854C"
        ),
        "window_is_below_cosby_measured_electron_energies": True,
        "state_branching_measured_or_inferred": False,
        "ground_atom_or_partner_injected_into_chemistry": False,
    }


def collision_set_targets(dataset_id: str) -> tuple[str, ...]:
    try:
        return DATASET_TARGETS[dataset_id]
    except KeyError as error:
        raise ValueError(f"Unsupported collision set: {dataset_id}") from error


def _load_processes(dataset_id: str = DATASET_ID) -> list[dict[str, Any]]:
    processes = _load_morgan_processes()
    if dataset_id in {IST_H2_COMMUNITY_CH4_DATASET_ID, IST_CO2_O2_DATASET_ID}:
        processes = [item for item in processes if item["target"] != "H2"]
        processes.extend(_load_ist_h2_processes())
    if dataset_id in {COMMUNITY_CH4_DATASET_ID, IST_H2_COMMUNITY_CH4_DATASET_ID}:
        processes = [item for item in processes if item["target"] != "CH4"]
        processes.extend(_load_community_ch4_processes())
    if dataset_id in {
        HYBRID_DATASET_ID,
        REACTIVE_N2_DATASET_ID,
        COMMUNITY_CH4_DATASET_ID,
        IST_H2_COMMUNITY_CH4_DATASET_ID,
        IST_CO2_O2_DATASET_ID,
    }:
        processes.extend(_load_phelps_n2_processes())
        if dataset_id in {
            REACTIVE_N2_DATASET_ID, COMMUNITY_CH4_DATASET_ID,
            IST_H2_COMMUNITY_CH4_DATASET_ID, IST_CO2_O2_DATASET_ID,
        }:
            processes.append(_load_tabata_n2_dissociation_process())
    elif dataset_id != DATASET_ID:
        raise ValueError(f"Unsupported collision set: {dataset_id}")
    if dataset_id == IST_CO2_O2_DATASET_ID:
        processes.extend(_load_ist_co2_o2_processes())
    return processes


def _dataset_metadata(dataset_id: str) -> tuple[str, str]:
    if dataset_id == IST_CO2_O2_DATASET_ID:
        components = (
            DATASET_SHA256, IST_H2_SHA256, PHELPS_N2_SHA256,
            TABATA_N2_SHA256, IST_CO2_O2_SHA256,
        )
        composite_sha256 = hashlib.sha256("|".join(components).encode("ascii")).hexdigest()
        return (
            composite_sha256,
            f"{SOURCE_REFERENCE}; {IST_H2_SOURCE_REFERENCE}; "
            f"{PHELPS_SOURCE_REFERENCE}; {TABATA_SOURCE_REFERENCE}; "
            f"{IST_CO2_O2_SOURCE_REFERENCE}",
        )
    if dataset_id == IST_H2_COMMUNITY_CH4_DATASET_ID:
        return (
            IST_H2_COMMUNITY_CH4_DATASET_SHA256,
            f"{SOURCE_REFERENCE}; {IST_H2_SOURCE_REFERENCE}; "
            f"{COMMUNITY_CH4_SOURCE_REFERENCE}; {PHELPS_SOURCE_REFERENCE}; "
            f"{TABATA_SOURCE_REFERENCE}",
        )
    if dataset_id == COMMUNITY_CH4_DATASET_ID:
        return (
            COMMUNITY_CH4_DATASET_SHA256,
            f"{SOURCE_REFERENCE}; {COMMUNITY_CH4_SOURCE_REFERENCE}; "
            f"{PHELPS_SOURCE_REFERENCE}; {TABATA_SOURCE_REFERENCE}",
        )
    if dataset_id == REACTIVE_N2_DATASET_ID:
        return (
            REACTIVE_N2_DATASET_SHA256,
            f"{SOURCE_REFERENCE}; {PHELPS_SOURCE_REFERENCE}; {TABATA_SOURCE_REFERENCE}",
        )
    if dataset_id == HYBRID_DATASET_ID:
        return (
            HYBRID_DATASET_SHA256,
            f"{SOURCE_REFERENCE}; {PHELPS_SOURCE_REFERENCE}",
        )
    if dataset_id == DATASET_ID:
        return DATASET_SHA256, SOURCE_REFERENCE
    raise ValueError(f"Unsupported collision set: {dataset_id}")


def _has_explicit_neutral_fragments(process: dict[str, Any]) -> bool:
    """Only count declared neutral products, not charge markers or proposed branches."""

    label = str(process.get("product", "")).strip()
    comment = str(process.get("comment", "")).lower()
    fragments = split_lxcat_products(label)
    if "suggested dissociation channel" in comment or len(fragments) <= 1:
        return False
    return all(
        fragment
        and not re.search(r"\^[+-]|\([+-],|[+-]$", fragment)
        for fragment in fragments
    )


def _has_resolved_ground_fragments(process: dict[str, Any]) -> bool:
    """Distinguish declared fragment labels from ground-state injectability."""

    if process.get("product_state_unresolved") or not _has_explicit_neutral_fragments(process):
        return False
    for fragment in split_lxcat_products(str(process["product"])):
        if "*" in fragment:
            return False
        for note in re.findall(r"\(([^)]*)\)", fragment):
            if re.fullmatch(r"\d+(?:\.\d+)?\s*eV", note, flags=re.IGNORECASE):
                continue
            if fragment.startswith("H(") and note.lower() == "1s":
                continue
            if (fragment.startswith("O(") and note == "3P") or (
                fragment.startswith(("CO(", "O2(")) and note == "X"
            ) or (fragment.startswith("CO2(") and note == "X,v=000"):
                continue
            return False
    return True


def _coverage(processes: list[dict[str, Any]], targets_required: tuple[str, ...]) -> dict[str, Any]:
    counts = Counter((item["target"], item["kind"]) for item in processes)
    targets: dict[str, Any] = {}
    for target in targets_required:
        target_processes = [item for item in processes if item["target"] == target]
        neutral_candidates = [
            item for item in target_processes
            if item["kind"] in {"EXCITATION", "DISSOCIATION"}
        ]
        explicit_neutral_fragments = [
            str(item.get("product", ""))
            for item in neutral_candidates
            if _has_explicit_neutral_fragments(item)
        ]
        resolved_ground_fragments = [
            str(item.get("product", ""))
            for item in neutral_candidates
            if _has_resolved_ground_fragments(item)
        ]
        unresolved_dissociations = [
            str(item.get("product", ""))
            for item in neutral_candidates
            if "dissoc" in str(item.get("product", "")).lower()
            and str(item.get("product", "")) not in explicit_neutral_fragments
        ]
        channel_counts = {
            kind.lower(): counts[(target, kind)]
            for kind in (
                "ELASTIC",
                "EFFECTIVE",
                "MOMENTUM",
                "EXCITATION",
                "DISSOCIATION",
                "IONIZATION",
                "ATTACHMENT",
            )
        }
        momentum = (
            channel_counts["elastic"]
            + channel_counts["effective"]
            + channel_counts["momentum"]
        ) > 0
        inelastic = channel_counts["excitation"] > 0
        ionization = channel_counts["ionization"] > 0
        targets[target] = {
            "complete": momentum and inelastic and ionization,
            "completeness_scope": "electron-energy-distribution input only",
            "momentum_transfer_present": momentum,
            "excitation_present": inelastic,
            "neutral_dissociation_present": (
                channel_counts["dissociation"] > 0 or bool(explicit_neutral_fragments)
            ),
            "explicit_neutral_fragment_labels": explicit_neutral_fragments,
            "explicit_neutral_fragment_channel_count": len(explicit_neutral_fragments),
            "has_explicit_neutral_fragment_source": bool(explicit_neutral_fragments),
            "resolved_ground_state_fragment_labels": resolved_ground_fragments,
            "has_resolved_ground_state_fragment_source": bool(resolved_ground_fragments),
            "unresolved_dissociation_labels": unresolved_dissociations,
            "unresolved_dissociation_channel_count": len(unresolved_dissociations),
            "ionization_present": ionization,
            "attachment_present": channel_counts["attachment"] > 0,
            "channel_counts": channel_counts,
            "process_count": sum(channel_counts.values()),
        }
    return {
        "all_required_targets_complete": all(item["complete"] for item in targets.values()),
        "completeness_scope": "EEDF input coverage; reactive fragment and heavy-chemistry closure are separate gates",
        "any_unresolved_dissociation_labels": any(
            item["unresolved_dissociation_channel_count"] for item in targets.values()
        ),
        "all_required_targets_have_explicit_neutral_fragment_source": all(
            item["has_explicit_neutral_fragment_source"] for item in targets.values()
        ),
        "all_required_targets_have_resolved_ground_state_fragment_source": all(
            item["has_resolved_ground_state_fragment_source"] for item in targets.values()
        ),
        "targets": targets,
        "process_count": len(processes),
        "criteria": [
            "elastic/effective/momentum transfer (including BOLSIG MOMENTUM)",
            "inelastic excitation or dissociation",
            "ionization",
            "attachment retained when supplied by the contributor",
            "neutral dissociation outcome retained when independently tabulated",
        ],
    }


def collision_dataset_status(dataset_id: str = DATASET_ID) -> dict[str, Any]:
    try:
        processes = _load_processes(dataset_id)
        targets = collision_set_targets(dataset_id)
        coverage = _coverage(processes, targets)
        sha256, source = _dataset_metadata(dataset_id)
        return {
            "available": True,
            "dataset_id": dataset_id,
            "sha256": sha256,
            "source": source,
            "format": (
                "hybrid private LXCat Morgan/IST-Lisbon H2/IST-Lisbon CO2/O2, "
                "Cantera Phelps N2, and QST/Tabata N2 outcome table"
                if dataset_id == IST_CO2_O2_DATASET_ID
                else
                "hybrid private LXCat Morgan, private LXCat Community/Song-Bouwman CH4, "
                "Cantera Phelps N2, and QST/Tabata N2 outcome table"
                if dataset_id == COMMUNITY_CH4_DATASET_ID
                else "hybrid private LXCat Morgan/IST-Lisbon H2/Community CH4, Cantera Phelps N2, and QST/Tabata N2 outcome table"
                if dataset_id == IST_H2_COMMUNITY_CH4_DATASET_ID
                else
                "hybrid LXCat BOLSIG+ text, Cantera Phelps YAML, and QST/Tabata outcome table"
                if dataset_id == REACTIVE_N2_DATASET_ID
                else "hybrid LXCat BOLSIG+ text plus Cantera Phelps YAML"
                if dataset_id == HYBRID_DATASET_ID
                else "LXCat BOLSIG+ complete-set text"
            ),
            "redistribution": (
                "private Morgan and IST-Lisbon contributor copies; no LXCat tables redistributed"
                if dataset_id == IST_CO2_O2_DATASET_ID
                else
                "private Morgan and Community contributor copies; no LXCat contributor "
                "tables may be redistributed with source"
                if dataset_id == COMMUNITY_CH4_DATASET_ID
                else "private Morgan, IST-Lisbon H2 and Community CH4 contributor copies; no LXCat tables redistributed"
                if dataset_id == IST_H2_COMMUNITY_CH4_DATASET_ID
                else
                "private Morgan copy, Cantera-installed Phelps data, and a public QST recommended outcome table"
                if dataset_id == REACTIVE_N2_DATASET_ID
                else "private Morgan copy plus Cantera-installed Phelps data; no contributor data copied into the project"
                if dataset_id == HYBRID_DATASET_ID
                else "private local copy; contributor data are not bundled with source"
            ),
            "coverage": coverage,
            "detail": "Hash verified and structurally parsed by BOLOS",
            "scientific_scope": (
                "Opt-in CO2/O2 complete-set EEDF input only; the H2/O2 rotational-free "
                "overlap is restricted to 10-200 Td. BOLOS does not use the CO2(010) "
                "statistical-weight flag for superelastic kinetics, so low-field behavior "
                "is outside this implementation. Excited/ionic products are not "
                "automatically grounded or chemically propagated."
                if dataset_id == IST_CO2_O2_DATASET_ID
                else
                "Opt-in H2 whole-set replacement only. Contributor-reported swarm "
                "agreement applies to 1-200 Td without rotational channels; "
                "an independent project swarm comparison is pending. The default "
                "BOLOS grid can yield nonfinite rates at 1-10 Td and rejects them. Explicit "
                "H2 neutral fragment labels contain excited H and are excluded "
                "from ground-state chemistry pending quenching kinetics."
                if dataset_id == IST_H2_COMMUNITY_CH4_DATASET_ID
                else "EEDF input coverage does not validate product branching or wet-apparatus yield."
            ),
        }
    except (ValueError, OSError) as error:
        return {
            "available": False,
            "dataset_id": dataset_id,
            "sha256": _dataset_metadata(dataset_id)[0],
            "source": _dataset_metadata(dataset_id)[1],
            "detail": str(error),
        }


def _effective_clamp_targets(
    processes: list[dict[str, Any]], targets_required: tuple[str, ...]
) -> list[str]:
    """Predict BOLOS' effective-to-elastic clipping diagnostic.

    Some contributor sets define an effective momentum curve that is slightly
    below the sum of separately tabulated inelastic channels at a few energies.
    BOLOS clips the reconstructed elastic curve to zero there. This is exposed
    in results instead of being hidden in a log message.
    """

    clipped: list[str] = []
    for target in targets_required:
        effective = [
            item
            for item in processes
            if item["target"] == target and item["kind"] in {"EFFECTIVE", "MOMENTUM"}
        ]
        inelastic = [
            item
            for item in processes
            if item["target"] == target
            and item["kind"] in {"EXCITATION", "IONIZATION", "ATTACHMENT"}
        ]
        if len(effective) != 1:
            continue
        grid = np.asarray(effective[0]["data"], dtype=float)[:, 0]
        reconstructed = np.asarray(effective[0]["data"], dtype=float)[:, 1].copy()
        for item in inelastic:
            data = np.asarray(item["data"], dtype=float)
            reconstructed -= np.interp(grid, data[:, 0], data[:, 1], left=0.0, right=0.0)
        if float(reconstructed.min()) < 0:
            clipped.append(target)
    return clipped


def _cantera_reference(
    fractions: dict[str, float], temperature_K: float, field_Td: float
) -> dict[str, Any]:
    unsupported = sorted(set(fractions) - {"CH4", "H2O", "N2"})
    if unsupported:
        return {
            "available_for_this_mixture": False,
            "engine": "Cantera 3.2 Boltzmann-two-term",
            "unsupported_targets": unsupported,
            "reason": "The installed Pavan reference mechanism has no H2 or NH3 electron-collision channels.",
        }
    mechanism = next(
        (
            Path(directory) / "example_data" / "methane-plasma-pavan-2023.yaml"
            for directory in ct.get_data_directories()
            if (Path(directory) / "example_data" / "methane-plasma-pavan-2023.yaml").is_file()
        ),
        None,
    )
    if mechanism is None:
        return {
            "available_for_this_mixture": False,
            "engine": "Cantera 3.2 Boltzmann-two-term",
            "reason": "Cantera's packaged Pavan methane-plasma reference is unavailable.",
        }
    phase = ct.Solution(str(mechanism), "gri30-plasma", transport_model=None)
    composition = dict(fractions)
    composition["e"] = 1e-12
    phase.TPX = temperature_K, ct.one_atm, composition
    phase.reduced_electric_field = field_Td * 1e-21
    phase.update_electron_energy_distribution()
    return {
        "available_for_this_mixture": True,
        "engine": "Cantera 3.2 Boltzmann-two-term",
        "mechanism": mechanism.name,
        "mean_energy_eV": float(phase.mean_electron_energy),
        "comparison_scope": "independent trend check; different collision data and chemical scope",
    }


@lru_cache(maxsize=24)
def _same_set_cantera_yaml(
    dataset_id: str, targets: tuple[str, ...], energy_max_eV: float, energy_cells: int
) -> str:
    """Convert hash-checked collision inputs to a transient Cantera phase.

    No contributor tables are written to disk. The post-solve Tabata outcome
    is intentionally omitted here, just as it is omitted from BOLOS' EEDF
    energy-loss assembly.
    """

    selected = [
        item for item in _load_processes(dataset_id)
        if item["target"] in targets and not item.get("outcome_only")
    ]
    collisions = []
    for item in selected:
        rows = item["data"]
        collision = {
            "target": item["target"],
            "kind": str(item["kind"]).lower(),
            "energy-levels": [float(row[0]) for row in rows],
            "cross-sections": [float(row[1]) for row in rows],
        }
        if item.get("product") is not None:
            collision["product"] = str(item["product"])
        if item.get("threshold") is not None:
            collision["threshold"] = float(item["threshold"])
        collisions.append(collision)
    elements = sorted({element for target in targets for element in {
        "H2": {"H"}, "H2O": {"H", "O"}, "CH4": {"C", "H"},
        "NH3": {"N", "H"}, "N2": {"N"}, "CO2": {"C", "O"},
        "O2": {"O"},
    }[target]})
    document = {
        "phases": [{
            "name": "same-collision-set-eedf",
            "thermo": "plasma",
            "kinetics": "gas",
            "reactions": "none",
            "elements": [*elements, "E"],
            "species": [{"nasa_gas.yaml/species": ["Electron", *targets]}],
            "state": {"T": 300.0, "P": 101325.0, "X": {targets[0]: 1.0}},
            "electron-energy-distribution": {
                "type": "Boltzmann-two-term",
                "energy-levels": [
                    energy_max_eV * index / energy_cells
                    for index in range(energy_cells + 1)
                ],
            },
        }],
        "electron-collisions": collisions,
    }
    return yaml.safe_dump(document, sort_keys=False, allow_unicode=True)


def _same_set_cantera_reference(
    dataset_id: str,
    fractions: dict[str, float],
    temperature_K: float,
    field_Td: float,
    energy_max_eV: float,
    energy_cells: int,
    bolos_mean_energy_eV: float,
    bolos_rates: list[dict[str, Any]],
) -> dict[str, Any]:
    """Compare independent BOLOS/Cantera solvers on identical collision data."""

    targets = tuple(target for target in collision_set_targets(dataset_id) if target in fractions)
    try:
        source = _same_set_cantera_yaml(dataset_id, targets, energy_max_eV, energy_cells)
        phase = ct.Solution(
            yaml=source, name="same-collision-set-eedf", transport_model=None
        )
        composition = {**fractions, "Electron": 1e-12}
        phase.TPX = temperature_K, ct.one_atm, composition
        phase.reduced_electric_field = field_Td * 1e-21
        phase.update_electron_energy_distribution()
        mean_energy = float(phase.mean_electron_energy)
        energies = np.asarray(phase.electron_energy_levels, dtype=float)
        distribution = np.asarray(phase.electron_energy_distribution, dtype=float)
        probability = np.sqrt(energies) * distribution
        probability_norm = float(np.trapezoid(probability, energies))
        if not np.isfinite(probability_norm) or probability_norm <= 0:
            raise ValueError("Cantera EEDF normalization is not finite and positive")
        probability /= probability_norm
        electron_speed = np.sqrt(2.0 * elementary_charge * energies / electron_mass)
        process_lookup = {
            (item["target"], str(item["kind"]).lower(), str(item.get("product")).strip()): item
            for item in _load_processes(dataset_id)
            if item["target"] in targets and str(item["kind"]).upper() not in {"ELASTIC", "EFFECTIVE"}
        }
        rate_rows = []
        for rate in bolos_rates:
            key = (rate["target"], rate["kind"], str(rate["product_label"]).strip())
            process = process_lookup.get(key)
            if process is None:
                continue
            curve = np.asarray(process["data"], dtype=float)
            cross_section = np.interp(
                energies, curve[:, 0], curve[:, 1], left=0.0, right=0.0
            )
            cantera_rate = float(np.trapezoid(
                cross_section * electron_speed * probability, energies
            ))
            bolos_rate = float(rate["rate_coefficient_m3_s"])
            if not np.isfinite(cantera_rate) or cantera_rate < 0:
                raise ValueError("Cantera EEDF produced a nonfinite or negative collision rate")
            factor = (
                max(cantera_rate, bolos_rate) / min(cantera_rate, bolos_rate)
                if min(cantera_rate, bolos_rate) > 0 else None
            )
            rate_rows.append({
                "target": rate["target"],
                "kind": rate["kind"],
                "product_label": rate["product_label"],
                "bolos_rate_coefficient_m3_s": bolos_rate,
                "cantera_eedf_rate_coefficient_m3_s": cantera_rate,
                "factor_difference": factor,
                "source_weighted_rate_m3_s": max(cantera_rate, bolos_rate) * fractions[rate["target"]],
                "outcome_only": bool(rate.get("outcome_only")),
            })
        material_rate_floor_m3_s = 1e-20
        material = [
            row for row in rate_rows
            if max(row["bolos_rate_coefficient_m3_s"], row["cantera_eedf_rate_coefficient_m3_s"])
            >= material_rate_floor_m3_s
        ]
        within_ten = [
            row for row in material
            if row["factor_difference"] is not None and row["factor_difference"] <= 10.0
        ]
        rate_comparison = {
            "matched_process_count": len(rate_rows),
            "material_rate_floor_m3_s": material_rate_floor_m3_s,
            "material_process_count": len(material),
            "material_within_factor_ten_count": len(within_ten),
            "material_within_factor_ten_fraction": (
                len(within_ten) / len(material) if material else None
            ),
            "largest_source_weighted_channels": sorted(
                rate_rows, key=lambda row: row["source_weighted_rate_m3_s"], reverse=True
            )[:8],
            "largest_reactive_channel_by_target": [
                max(
                    (
                        row for row in rate_rows
                        if row["target"] == target
                        and (
                            row["kind"] != "excitation"
                            or "+" in str(row["product_label"])
                        )
                    ),
                    key=lambda row: row["source_weighted_rate_m3_s"],
                )
                for target in targets
                if any(
                    row["target"] == target
                    and (row["kind"] != "excitation" or "+" in str(row["product_label"]))
                    for row in rate_rows
                )
            ],
            "largest_material_disagreements": sorted(
                material,
                key=lambda row: row["factor_difference"] if row["factor_difference"] is not None else float("inf"),
                reverse=True,
            )[:8],
            "selected_n2_optical_channels": [
                row for row in rate_rows
                if row["target"] == "N2" and row["product_label"] in {
                    "N2(C3)", "N2^+(B2SIGMA)"
                }
            ],
            "comparison_scope": (
                "Collision rates use each independently solved EEDF folded through the "
                "same imported cross sections. The 1e-20 m3/s floor excludes only "
                "numerically tiny channels from the factor-ten summary; it does "
                "not affect either solver or source terms."
            ),
        }
    except (ct.CanteraError, ValueError) as error:
        return {
            "available_for_this_mixture": False,
            "engine": f"Cantera {ct.__version__} Boltzmann-two-term",
            "collision_set_dataset_id": dataset_id,
            "reason": str(error),
        }
    return {
        "available_for_this_mixture": True,
        "engine": f"Cantera {ct.__version__} Boltzmann-two-term",
        "collision_set_dataset_id": dataset_id,
        "collision_set_sha256": _dataset_metadata(dataset_id)[0],
        "selected_eedf_process_count": sum(
            1 for item in _load_processes(dataset_id)
            if item["target"] in targets and not item.get("outcome_only")
        ),
        "mean_energy_eV": mean_energy,
        "relative_mean_energy_difference": (
            bolos_mean_energy_eV - mean_energy
        ) / mean_energy,
        "rate_comparison": rate_comparison,
        "comparison_scope": (
            "independent two-term numerical implementation using the same "
            "collision curves, mixture, E/N, energy range and nominal bin "
            "resolution (the solvers discretize that range differently); not an external "
            "swarm experiment or product-yield validation"
        ),
        "private_collision_data_persisted": False,
    }


def solve_self_consistent_eedf(
    request: BoltzmannEEDFRequest, *, include_same_set_cantera_reference: bool = True
) -> dict[str, Any]:
    try:
        from bolos import grid, solver
    except ImportError as error:
        raise ValueError("BOLOS 0.2 is not installed") from error

    processes = _load_processes(request.collision_set)
    targets_required = collision_set_targets(request.collision_set)
    coverage = _coverage(processes, targets_required)
    dataset_sha256, dataset_source = _dataset_metadata(request.collision_set)
    fractions = {
        target: float(value)
        for target, value in request.mole_fractions.items()
        if value > 0
    }
    total = sum(fractions.values())
    fractions = {target: value / total for target, value in fractions.items()}
    incomplete = [
        target for target in fractions if not coverage["targets"][target]["complete"]
    ]
    if incomplete:
        raise ValueError(f"Collision coverage is incomplete for {sorted(incomplete)}")

    selected = [item for item in processes if item["target"] in fractions]
    eedf_processes = [item for item in selected if not item.get("outcome_only")]
    outcome_processes = [item for item in selected if item.get("outcome_only")]
    energy_grid = grid.LinearGrid(0.0, request.energy_max_eV, request.energy_cells)
    boltzmann = solver.BoltzmannSolver(energy_grid)
    root_logger = logging.getLogger()
    previous_level = root_logger.level
    try:
        root_logger.setLevel(logging.ERROR)
        boltzmann.load_collisions(eedf_processes)
        for target, fraction in fractions.items():
            boltzmann.set_density(target, fraction)
        boltzmann.kT = request.gas_temperature_K * solver.KB / solver.ELECTRONVOLT
        boltzmann.EN = request.reduced_electric_field_Td * solver.TOWNSEND
        boltzmann.init()
    finally:
        root_logger.setLevel(previous_level)
    initial = boltzmann.maxwell(request.initial_electron_temperature_eV)
    with warnings.catch_warnings(record=True) as convergence_warnings:
        warnings.simplefilter("always", RuntimeWarning)
        distribution, iterations, error = boltzmann.converge(
            initial,
            maxn=request.maximum_iterations,
            rtol=request.relative_tolerance,
            full=True,
        )

    energies = np.asarray(energy_grid.c, dtype=float)
    distribution = np.asarray(distribution, dtype=float)
    probability = np.sqrt(energies) * distribution
    normalization = float(scipy_integrate.simpson(probability, x=energies))
    mean_energy = float(boltzmann.mean_energy(distribution))
    if not np.isfinite(normalization) or normalization <= 0 or not np.isfinite(mean_energy) or mean_energy <= 0:
        raise ValueError("The converged EEDF did not produce finite positive normalization and mean energy")
    with warnings.catch_warnings(record=True) as mobility_warnings:
        # BOLOS evaluates nu/sqrt(E) at the zero-energy boundary before
        # explicitly zeroing that boundary's integrand in mobility().
        warnings.simplefilter("always", RuntimeWarning)
        reduced_mobility = float(boltzmann.mobility(distribution))
    if not np.isfinite(reduced_mobility) or reduced_mobility <= 0:
        raise ValueError("The converged EEDF did not produce a positive finite electron mobility")
    electron_temperature_K = 2.0 * mean_energy * elementary_charge / (3.0 * Boltzmann)

    rates = []
    rate_runtime_warnings: set[str] = set()
    for index, (target, process) in enumerate(boltzmann.iter_inelastic()):
        data = np.asarray(process.data, dtype=float)
        nonzero = np.flatnonzero(data[:, 1] > 0)
        with warnings.catch_warnings(record=True) as observed_warnings:
            warnings.simplefilter("always", RuntimeWarning)
            rate_coefficient = float(boltzmann.rate(distribution, process))
        rate_runtime_warnings.update(
            str(item.message) for item in observed_warnings
            if issubclass(item.category, RuntimeWarning)
        )
        if not np.isfinite(rate_coefficient) or rate_coefficient < 0:
            raise ValueError(
                "The converged EEDF produced a nonfinite or negative collision rate "
                f"for {target.name} at {request.reduced_electric_field_Td:g} Td "
                f"on the {request.energy_cells}-cell/{request.energy_max_eV:g} eV grid; "
                "the reported convergence flag alone is insufficient. Change the "
                "field/grid and verify numerical sensitivity before using rates."
            )
        rates.append(
            {
                "id": f"{target.name}_{process.kind.lower()}_{index:03d}",
                "target": target.name,
                "kind": process.kind.lower(),
                "product_label": process.product,
                "threshold_eV": float(process.threshold),
                "first_nonzero_energy_eV": (
                    float(data[nonzero[0], 0]) if nonzero.size else None
                ),
                "rate_coefficient_m3_s": rate_coefficient,
                "outcome_only": False,
            }
        )

    normalized_probability = probability / normalization
    electron_speed = np.sqrt(2.0 * elementary_charge * energies / electron_mass)
    for index, process in enumerate(outcome_processes):
        data = np.asarray(process["data"], dtype=float)
        cross_section = np.interp(
            energies, data[:, 0], data[:, 1], left=0.0, right=0.0
        )
        rate_integrand = cross_section * electron_speed * normalized_probability
        rate_coefficient = float(
            np.trapezoid(rate_integrand, energies)
        )
        product_state_energy_support = None
        ground_atom_energetic_window = None
        analytic_fit_diagnostic = None
        if process.get("product_state_unresolved"):
            # Cosby 1993 observed neutral N2 fragment translational states at
            # 18.5-148.5 eV. This partitions the *computed total event rate*
            # by source electron energy; it is not a branching measurement.
            product_state_energy_support = _cosby_product_state_energy_support(
                energies, rate_integrand, rate_coefficient
            )
            ground_atom_energetic_window = _n2_ground_atom_energetic_window(
                energies, rate_integrand, rate_coefficient
            )
        if process.get("tabata_analytic_fit_diagnostic_available"):
            analytic_sigma = _tabata_n2_analytic_fit_m2(energies)
            analytic_integrand = analytic_sigma * electron_speed * normalized_probability
            analytic_rate = float(np.trapezoid(analytic_integrand, energies))
            analytic_fit_diagnostic = {
                "source": "QST/Tabata official F3 fit, https://www.qst.go.jp/uploaded/attachment/8709.f",
                "rate_coefficient_m3_s": analytic_rate,
                "tabular_linear_over_analytic_fit_ratio": (
                    rate_coefficient / analytic_rate if analytic_rate > 0 else None
                ),
                "product_state_energy_support": _cosby_product_state_energy_support(
                    energies, analytic_integrand, analytic_rate
                ),
                "ground_atom_energetic_window": _n2_ground_atom_energetic_window(
                    energies, analytic_integrand, analytic_rate
                ),
                "used_for_eedf_or_chemical_source": False,
                "independent_experimental_measurement": False,
            }
        nonzero = np.flatnonzero(data[:, 1] > 0)
        rates.append(
            {
                "id": f"{process['target']}_dissociation_outcome_{index:03d}",
                "target": process["target"],
                "kind": "dissociation",
                "product_label": process["product"],
                "threshold_eV": float(process["threshold"]),
                "first_nonzero_energy_eV": (
                    float(data[nonzero[0], 0]) if nonzero.size else None
                ),
                "rate_coefficient_m3_s": rate_coefficient,
                "outcome_only": True,
                "product_state_unresolved": bool(process.get("product_state_unresolved")),
                "product_state_energy_support": product_state_energy_support,
                "ground_atom_energetic_window": ground_atom_energetic_window,
                "analytic_fit_diagnostic": analytic_fit_diagnostic,
                "source": process["comment"],
                "source_sha256": process["source_sha256"],
                "coupling": "post-solve outcome cross section; excluded from EEDF energy-loss assembly",
            }
        )

    cantera_reference = _cantera_reference(
        fractions, request.gas_temperature_K, request.reduced_electric_field_Td
    )
    if cantera_reference.get("available_for_this_mixture"):
        reference_mean = float(cantera_reference["mean_energy_eV"])
        cantera_reference["relative_mean_energy_difference"] = (
            mean_energy - reference_mean
        ) / reference_mean
    same_set_reference = (
        _same_set_cantera_reference(
            request.collision_set,
            fractions,
            request.gas_temperature_K,
            request.reduced_electric_field_Td,
            request.energy_max_eV,
            request.energy_cells,
            mean_energy,
            rates,
        )
        if include_same_set_cantera_reference
        else {
            "available_for_this_mixture": False,
            "engine": f"Cantera {ct.__version__} Boltzmann-two-term",
            "reason": "Per-cycle independent comparison omitted from the apparatus loop; run the direct EEDF endpoint to compare engines on the same collision set.",
        }
    )

    clipped = [
        target
        for target in _effective_clamp_targets(eedf_processes, targets_required)
        if target in fractions
    ]
    numerical_warnings = []
    convergence_runtime_warnings = sorted({
        str(item.message) for item in convergence_warnings
        if issubclass(item.category, RuntimeWarning)
    })
    if convergence_runtime_warnings:
        numerical_warnings.append(
            "BOLOS emitted floating-point warnings during EEDF iteration "
            f"({'; '.join(convergence_runtime_warnings)}); verify grid convergence."
        )
    mobility_runtime_warnings = sorted({
        str(item.message) for item in mobility_warnings
        if issubclass(item.category, RuntimeWarning)
        and str(item.message) != "divide by zero encountered in divide"
    })
    if mobility_runtime_warnings:
        numerical_warnings.append(
            "BOLOS emitted floating-point warnings during electron-mobility "
            f"integration ({'; '.join(mobility_runtime_warnings)}); verify grid convergence."
        )
    if clipped:
        numerical_warnings.append(
            "BOLOS clipped a negative reconstructed elastic momentum cross section to zero "
            f"for {', '.join(clipped)} where the effective curve was below summed inelastic curves."
        )
    if rate_runtime_warnings:
        numerical_warnings.append(
            "BOLOS emitted floating-point warnings while integrating collision rates "
            f"({'; '.join(sorted(rate_runtime_warnings))}); returned rates are finite, "
            "but this EEDF grid/field combination needs a numerical sensitivity check."
        )
    if (
        same_set_reference.get("available_for_this_mixture")
        and abs(float(same_set_reference["relative_mean_energy_difference"])) > 0.05
    ):
        numerical_warnings.append(
            "BOLOS and Cantera differ by more than 5% in mean electron energy "
            "on the same collision set; treat this EEDF point as solver-sensitive."
        )
    if same_set_reference.get("available_for_this_mixture"):
        rate_check = same_set_reference["rate_comparison"]
        if rate_check["material_within_factor_ten_fraction"] != 1.0:
            numerical_warnings.append(
                "At least one material collision channel differs by more than a factor "
                "of ten between BOLOS and the same-set Cantera EEDF. Mean-energy "
                "agreement alone is not enough to validate source rates."
            )

    try:
        bolos_version = version("bolos")
    except PackageNotFoundError:
        bolos_version = "unknown"
    return {
        "status": "self_consistent_eedf_solved",
        "solver": {
            "engine": "BOLOS",
            "version": bolos_version,
            "method": "steady-state two-term electron Boltzmann equation",
            "iterations": int(iterations),
            "final_relative_error": float(error),
            "converged": float(error) < request.relative_tolerance,
        },
        "conditions": {
            "mole_fractions": fractions,
            "gas_temperature_K": request.gas_temperature_K,
            "reduced_electric_field_Td": request.reduced_electric_field_Td,
            "energy_max_eV": request.energy_max_eV,
            "energy_cells": request.energy_cells,
            "note": "Pressure cancels when the field is specified as E/N; absolute rates still require number density.",
        },
        "collision_set": {
            "dataset_id": request.collision_set,
            "sha256": dataset_sha256,
            "source": dataset_source,
            "format": (
                "hybrid private LXCat Morgan/IST-Lisbon H2/IST-Lisbon CO2/O2, "
                "Cantera Phelps N2, and QST/Tabata N2 outcome table"
                if request.collision_set == IST_CO2_O2_DATASET_ID
                else
                "hybrid private LXCat Morgan, private LXCat Community/Song-Bouwman CH4, "
                "Cantera Phelps N2, and QST/Tabata N2 outcome table"
                if request.collision_set == COMMUNITY_CH4_DATASET_ID
                else "hybrid private LXCat Morgan/IST-Lisbon H2/Community CH4, Cantera Phelps N2, and QST/Tabata N2 outcome table"
                if request.collision_set == IST_H2_COMMUNITY_CH4_DATASET_ID
                else
                "hybrid LXCat BOLSIG+ text, Cantera Phelps YAML, and QST/Tabata outcome table"
                if request.collision_set == REACTIVE_N2_DATASET_ID
                else "hybrid LXCat BOLSIG+ text plus Cantera Phelps YAML"
                if request.collision_set == HYBRID_DATASET_ID
                else "LXCat BOLSIG+ complete-set text"
            ),
            "selected_process_count": len(selected),
            "eedf_process_count": len(eedf_processes),
            "post_solve_outcome_process_count": len(outcome_processes),
            "coverage": coverage,
        },
        "eedf": {
            "energy_levels_eV": energies.tolist(),
            "distribution_eV_minus_3_over_2": distribution.tolist(),
            "probability_density_eV_inv": probability.tolist(),
            "normalization": normalization,
            "mean_energy_eV": mean_energy,
            "electron_temperature_K": electron_temperature_K,
        },
        "transport": {
            "reduced_electron_mobility_m_inv_V_inv_s_inv": reduced_mobility,
            "definition": "electron mobility multiplied by gas number density",
        },
        "rate_coefficients": rates,
        "independent_reference": cantera_reference,
        "same_collision_set_reference": same_set_reference,
        "scientific_firewall": {
            "known_products_accessible_to_solver": False,
            "reaction_recipe_used": False,
            "self_consistent_eedf_solved": True,
            "complete_starting_gas_collision_coverage": coverage[
                "all_required_targets_complete"
            ],
            "heavy_species_chemistry_complete": False,
            "quantitative_yield_claim_allowed": False,
        },
        "numerical_warnings": numerical_warnings,
        "limitations": [
            "A spatially uniform steady-state two-term EEDF does not resolve streamer, sheath or electrode geometry.",
            "Collision rates are not product yields; a charge-balanced heavy-species mechanism and reactor coupling remain required.",
            "The Morgan and Phelps contributor sets retain provenance but have not yet passed an external swarm-coefficient benchmark in this project.",
            "The QST/Tabata total N + N product cross section is evaluated after EEDF convergence to avoid double-counting Phelps electronic energy loss. Its electronic product-state branching is unresolved; the outcome must not be injected as two ground-state atoms.",
            *(
                [
                    "The opt-in IST-Lisbon H2 set replaces every Morgan H2 curve. Its excited H fragment labels remain unresolved chemistry, and independent project swarm comparison is pending; contributor-reported validation without rotational channels covers 1-200 Td, while this BOLOS/grid implementation rejects nonfinite low-field rates."
                ]
                if request.collision_set == IST_H2_COMMUNITY_CH4_DATASET_ID
                else []
            ),
            *(
                [
                    "The IST-Lisbon CO2/O2 export supplies EEDF energy-loss and rate curves, "
                    "not a validated wet-reactor mechanism. Its CO2(010) statistical-weight "
                    "flag is not converted by BOLOS into a separate superelastic collision; "
                    "CO2/O2 excited and charged fragments require explicit downstream kinetics."
                ]
                if request.collision_set == IST_CO2_O2_DATASET_ID
                else []
            ),
        ],
    }
