from __future__ import annotations

import hashlib
import json
from typing import Any

import cantera as ct
import numpy as np
from scipy.constants import Avogadro, electron_mass, elementary_charge

from .schema import (
    ElectronCollisionDatasetRequest,
    ElectronCollisionRateRequest,
    ElectronEnergyDistribution,
)


def _eedf_values(
    eedf: ElectronEnergyDistribution, energy_eV: np.ndarray
) -> tuple[np.ndarray, float, float]:
    """Return probability density on an energy grid, its mean, and support coverage.

    The density is normalized per electron volt. Maxwellian mode is an explicit
    approximation; tabulated mode preserves the user's declared shape after
    normalization.
    """

    if eedf.mode == "maxwellian":
        mean_eV = eedf.mean_energy.as_unit("eV")
        temperature_eV = 2.0 * mean_eV / 3.0
        probability = (
            2.0
            / np.sqrt(np.pi)
            * np.sqrt(np.maximum(energy_eV, 0.0))
            / temperature_eV**1.5
            * np.exp(-energy_eV / temperature_eV)
        )
        # Coverage is integrated over the cross-section support, not silently
        # renormalized to it. Missing high-energy cross sections therefore stay visible.
        coverage = float(np.trapezoid(probability, energy_eV))
        return probability, mean_eV, coverage

    source_energy = np.asarray(eedf.energy_levels_eV, dtype=float)
    source_density = np.asarray(eedf.probability_density_eV_inv, dtype=float)
    normalization = float(np.trapezoid(source_density, source_energy))
    if normalization <= 0:
        raise ValueError("Tabulated EEDF has zero integrated probability")
    normalized = source_density / normalization
    mean_eV = float(np.trapezoid(source_energy * normalized, source_energy))
    probability = np.interp(energy_eV, source_energy, normalized, left=0.0, right=0.0)
    coverage = float(np.trapezoid(probability, energy_eV))
    return probability, mean_eV, coverage


def calculate_rate_coefficient(
    dataset: ElectronCollisionDatasetRequest,
    eedf: ElectronEnergyDistribution,
) -> dict[str, Any]:
    """Calculate the binary electron-target coefficient <sigma v>.

    Cross sections are in m², electron speed in m/s, and the EEDF is a
    probability density per eV, yielding m³/s per electron-target pair.
    """

    validate_electron_collision_dataset(dataset)
    source_energy = np.asarray(dataset.energy_levels_eV, dtype=float)
    upper = float(source_energy[-1])
    dense = np.linspace(float(source_energy[0]), upper, 8193)
    if eedf.mode == "tabulated":
        eedf_energy = np.asarray(eedf.energy_levels_eV, dtype=float)
        inside = eedf_energy[(eedf_energy >= source_energy[0]) & (eedf_energy <= upper)]
        integration_energy = np.unique(np.concatenate((dense, source_energy, inside)))
    else:
        integration_energy = np.unique(np.concatenate((dense, source_energy)))
    cross_section = np.interp(
        integration_energy,
        source_energy,
        np.asarray(dataset.cross_sections_m2, dtype=float),
        left=0.0,
        right=0.0,
    )
    probability, mean_eV, support_probability = _eedf_values(eedf, integration_energy)
    electron_speed = np.sqrt(2.0 * elementary_charge * integration_energy / electron_mass)
    coefficient_m3_s = float(
        np.trapezoid(cross_section * electron_speed * probability, integration_energy)
    )
    return {
        "rate_coefficient_m3_s": coefficient_m3_s,
        "rate_coefficient_L_mol_s": coefficient_m3_s * 1000.0 * Avogadro,
        "eedf_mode": eedf.mode,
        "eedf_mean_energy_eV": mean_eV,
        "cross_section_support_probability": support_probability,
        "integration_points": int(integration_energy.size),
        "equation": "<sigma*v> = integral sigma(E)*sqrt(2*e*E/m_e)*p(E) dE",
    }


def calculate_electron_collision_rate(request: ElectronCollisionRateRequest) -> dict[str, Any]:
    validated = validate_electron_collision_dataset(request.dataset)
    coefficient = calculate_rate_coefficient(request.dataset, request.eedf)
    electron_density_m3 = request.electron_number_density.as_unit("1/m^3")
    pseudo_first_order_s_inv = coefficient["rate_coefficient_m3_s"] * electron_density_m3
    lifetime_s = 1.0 / pseudo_first_order_s_inv if pseudo_first_order_s_inv > 0 else None
    support_complete = coefficient["cross_section_support_probability"] >= 0.999
    complete_target_set = request.dataset.coverage == "complete_target_set"
    return {
        "status": "electron_collision_rate_calculated",
        "dataset": validated,
        "eedf": {
            "mode": request.eedf.mode,
            "mean_energy_eV": coefficient["eedf_mean_energy_eV"],
            "source": request.eedf.source,
            "method": request.eedf.method,
            "basis": request.eedf.basis,
            "cross_section_support_probability": coefficient[
                "cross_section_support_probability"
            ],
        },
        "electron_number_density_m3": electron_density_m3,
        "rate_coefficient_m3_s": coefficient["rate_coefficient_m3_s"],
        "rate_coefficient_L_mol_s": coefficient["rate_coefficient_L_mol_s"],
        "pseudo_first_order_rate_s_inv": pseudo_first_order_s_inv,
        "target_lifetime_s": lifetime_s,
        "integration_points": coefficient["integration_points"],
        "equation": coefficient["equation"],
        "scientific_firewall": {
            "known_products_accessible_to_discovery": False,
            "rate_coefficient_calculated": True,
            "cross_section_support_complete": support_complete,
            "complete_target_collision_set_declared": complete_target_set,
            "quantitative_yield_claim_allowed": False,
        },
        "limitations": [
            "The coefficient is conditional on the declared EEDF and electron density.",
            "One channel does not establish complete plasma chemistry or a predictive product yield.",
        ],
    }


def validate_electron_collision_dataset(
    request: ElectronCollisionDatasetRequest,
) -> dict[str, Any]:
    """Validate one provenance-labelled collision channel with Cantera.

    Validation establishes schema and runtime compatibility only. A collision
    rate still requires an electron-energy distribution and a complete plasma
    phase/mechanism.
    """

    energy = np.asarray(request.energy_levels_eV, dtype=float)
    cross_section = np.asarray(request.cross_sections_m2, dtype=float)
    try:
        rate = ct.ElectronCollisionPlasmaRate(energy, cross_section)
    except (ct.CanteraError, TypeError, ValueError) as error:
        raise ValueError(f"Cantera rejected the electron-collision dataset: {error}") from error
    canonical = json.dumps(
        request.model_dump(mode="json", exclude_none=True),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    nonzero = np.flatnonzero(cross_section > 0)
    threshold = float(energy[nonzero[0]]) if nonzero.size else None
    return {
        "status": "validated_cross_section_dataset",
        "dataset_id": request.id,
        "sha256": hashlib.sha256(canonical).hexdigest(),
        "cantera_version": ct.__version__,
        "cantera_rate_type": type(rate).__name__,
        "points": int(energy.size),
        "energy_range_eV": [float(energy[0]), float(energy[-1])],
        "first_nonzero_energy_eV": threshold,
        "maximum_cross_section_m2": float(cross_section.max()),
        "integrated_cross_section_m2_eV": float(np.trapezoid(cross_section, energy)),
        "provenance": {
            "source": request.source,
            "method": request.method,
            "basis": request.basis,
            "mechanism_origin": request.mechanism_origin,
            "discovery_manifest_sha256": request.discovery_manifest_sha256,
            "coverage": request.coverage,
        },
        "scientific_firewall": {
            "known_products_accessible_to_discovery": False,
            "blind_discovery_eligible": request.mechanism_origin == "blind_discovery_artifact",
            "usable_for_rate_without_eedf": False,
            "quantitative_yield_claim_allowed": False,
        },
        "limitations": [
            "This validates one collision cross-section curve; it does not solve an electron-energy distribution.",
            "A predictive plasma mechanism also requires all important collision and heavy-species channels, thermochemistry, electron density and reactor coupling.",
        ],
    }
