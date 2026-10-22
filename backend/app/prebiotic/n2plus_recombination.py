"""Bounded N2+ dissociative-recombination audit, not a chemical source.

The private LXCat record stops below 0.4 eV and does not identify neutral
fragment electronic states. Folding it through a solved EEDF can diagnose only
the rate contribution in their common energy interval. It cannot supply an
N-atom source term or a total spark recombination rate.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import numpy as np
from scipy.constants import electron_mass, elementary_charge


PRIVATE_RECORD = (
    Path(__file__).resolve().parents[2]
    / "data" / "collision" / "private"
    / "eMol_LeHavre_N2plus_2026-09-27" / "Cross section.txt"
)
REVIEWED_SHA256 = "00ad61fca159d92a5adb00f7f80fb2cb9cf463e7431c9bb7b6211747f40e2967"
SOURCE = (
    "eMol-LeHavre database, www.lxcat.net, retrieved September 27, 2026; "
    "D. A. Little et al., Phys. Rev. A 90, 052705 (2014), "
    "DOI 10.1103/PhysRevA.90.052705"
)


def load_private_n2plus_recombination() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read the reviewed record; bracket conflicting rounded-energy rows.

    The LXCat export repeats 8,100 energy values, usually with different
    cross sections. Neither first-wins interpolation nor averaging is a
    defensible reconstruction of the unpublished precision. Min/max values
    at each rounded energy give a diagnostic envelope, not a confidence band.
    """

    if not PRIVATE_RECORD.is_file():
        raise ValueError("The private LXCat eMol-LeHavre N2+ record is not installed")
    actual = hashlib.sha256(PRIVATE_RECORD.read_bytes()).hexdigest()
    if actual != REVIEWED_SHA256:
        raise ValueError("The private LXCat N2+ record differs from its reviewed SHA-256")

    from bolos import parser

    with PRIVATE_RECORD.open(encoding="utf-8") as stream:
        processes = parser.parse(stream)
    if len(processes) != 1:
        raise ValueError("Expected one reviewed N2+ recombination curve")
    process = processes[0]
    data = np.asarray(process["data"], dtype=float)
    if (
        process.get("target") != "N2+"
        or process.get("kind") != "ATTACHMENT"
        or data.shape != (39999, 2)
        or not np.all(np.isfinite(data))
        or np.any(np.diff(data[:, 0]) < 0)
        or np.any(data[:, 1] < 0)
        or not np.isclose(data[0, 0], 1e-5)
        or not np.isclose(data[-1, 0], 0.39999)
    ):
        raise ValueError("The reviewed N2+ record has an unexpected process or energy grid")
    unique_energy, starts, counts = np.unique(
        data[:, 0], return_index=True, return_counts=True
    )
    if len(unique_energy) != 31899 or int(np.sum(counts - 1)) != 8100:
        raise ValueError("The reviewed N2+ energy rounding pattern changed")
    sigma_low = np.minimum.reduceat(data[:, 1], starts)
    sigma_high = np.maximum.reduceat(data[:, 1], starts)
    return unique_energy, sigma_low, sigma_high


def audit_n2plus_recombination(eedf: dict[str, Any]) -> dict[str, Any]:
    """Integrate only where a solved EEDF and the private record overlap.

    Missing sub-grid electron energies and the unmeasured higher-energy curve
    are never extrapolated. The result is a conditional *partial* coefficient,
    not a complete recombination rate or an HCN source.
    """

    record_energy, sigma_low, sigma_high = load_private_n2plus_recombination()
    eedf_energy = np.asarray(eedf["energy_levels_eV"], dtype=float)
    density = np.asarray(eedf["probability_density_eV_inv"], dtype=float).copy()
    if (
        eedf_energy.ndim != 1
        or eedf_energy.size < 2
        or density.shape != eedf_energy.shape
        or not np.all(np.isfinite(eedf_energy))
        or not np.all(np.isfinite(density))
        or not np.all(np.diff(eedf_energy) > 0)
        or eedf_energy[0] < 0
        or np.any(density < 0)
    ):
        raise ValueError("A finite, nonnegative, increasing tabulated EEDF is required")
    norm = float(np.trapezoid(density, eedf_energy))
    if norm <= 0:
        raise ValueError("The EEDF has zero probability")
    density /= norm

    lower = max(float(record_energy[0]), float(eedf_energy[0]))
    upper = min(float(record_energy[-1]), float(eedf_energy[-1]))
    if upper <= lower:
        partial_low = 0.0
        partial_high = 0.0
        probability = 0.0
    else:
        points = np.unique(np.concatenate((
            [lower, upper],
            record_energy[(record_energy > lower) & (record_energy < upper)],
            eedf_energy[(eedf_energy > lower) & (eedf_energy < upper)],
        )))
        p = np.interp(points, eedf_energy, density)
        speed = np.sqrt(2 * elementary_charge * points / electron_mass)
        partial_low = float(np.trapezoid(
            np.interp(points, record_energy, sigma_low) * speed * p, points
        ))
        partial_high = float(np.trapezoid(
            np.interp(points, record_energy, sigma_high) * speed * p, points
        ))
        probability = float(np.trapezoid(p, points))

    return {
        "record_sha256": REVIEWED_SHA256,
        "reference": SOURCE,
        "record_energy_window_eV": [float(record_energy[0]), float(record_energy[-1])],
        "rounded_duplicate_energy_rows": 8100,
        "duplicate_handling": "min/max envelope at each rounded energy; no averaging",
        "solved_eedf_energy_window_eV": [float(eedf_energy[0]), float(eedf_energy[-1])],
        "common_energy_window_eV": [lower, upper] if upper > lower else None,
        "unresolved_low_energy_window_eV": (
            [float(record_energy[0]), lower] if lower > record_energy[0] else None
        ),
        "record_ends_below_eedf_maximum": bool(record_energy[-1] < eedf_energy[-1]),
        "eedf_probability_in_common_window": probability,
        "partial_rate_coefficient_envelope_m3_s": [partial_low, partial_high],
        "total_rate_coefficient_m3_s": None,
        "neutral_product_states_resolved": False,
        "ion_vibrational_population_resolved": False,
        "allowed_as_chemical_source": False,
        "interpretation": (
            "Rounded-energy min/max envelope for the represented common interval only; "
            "not a statistical uncertainty, total rate, branching ratio, N-atom source, or yield."
        ),
    }
