"""Read-only first-collision audit for thermalized N(2D) near 298 K.

This does not assign a state branch to the total electron-impact N2
dissociation cross section and does not inject products into gas chemistry.
It asks what removes an N(2D) atom *if* one is formed and thermalized in a
homogeneous bath. The answer is not an HCN yield or a nascent-state model.
"""

from __future__ import annotations

import math
from typing import Any, Mapping

from scipy.constants import Avogadro


HERRON_1999_SOURCE = (
    "Herron, J. Phys. Chem. Ref. Data 28, 1453-1483 (1999), Tables 2 and 3; "
    "DOI 10.1063/1.556043; evaluated 298 K N(2D) bimolecular removal rates "
    "and multiplicative error factors"
)

# Species, k(298 K) [cm3 molecule-1 s-1], multiplicative error factor,
# fate classification. Table 2 covers H2, H2O, NH3, O, O2, N2, NO, N2O;
# Table 3 covers CH4, CO, CO2. Product branching is not inferred from a
# total removal rate. N(2D)+O transfers electronic excitation, but the
# O(3P)/O(1D) product-state split is unresolved; the tabulated temperature
# expression was checked over 300-400 K although k(298 K) is evaluated.
# N2 is the evaluated physical-quenching channel.
_CHANNELS = (
    ("H2", 2.2e-12, 1.3, "reactive_removal"),
    ("H2O", 4.0e-11, 2.0, "products_unresolved"),
    ("NH3", 5.0e-11, 2.0, "products_unresolved"),
    ("O", 1.4e-12, 2.0, "oxygen_product_state_unresolved"),
    ("O2", 5.2e-12, 1.3, "reactive_removal"),
    ("N2", 1.7e-14, 1.5, "physical_quench"),
    ("NO", 6.0e-11, 1.5, "reactive_removal"),
    ("N2O", 2.2e-12, 1.3, "reactive_removal"),
    ("CH4", 4.0e-12, 1.5, "products_unresolved"),
    ("CO", 1.9e-12, 1.3, "removal_or_quench"),
    ("CO2", 3.6e-13, 1.5, "reactive_removal"),
)


def n2d_fate_diagnostic(
    neutral_gas_mol: Mapping[str, float],
    *,
    gas_volume_L: float,
    gas_temperature_K: float,
    active_discharge_duration_s: float | None = None,
) -> dict[str, Any]:
    """Return conditional N(2D) removal hazards without assigning products.

    No fractions are normalized over an incomplete bath. Multiplicative
    rate-factor extrema are sensitivity bounds, not confidence intervals.
    """

    if not math.isfinite(gas_volume_L) or gas_volume_L <= 0:
        raise ValueError("N(2D) fate needs a positive finite gas volume")
    if not math.isfinite(gas_temperature_K) or gas_temperature_K <= 0:
        raise ValueError("N(2D) fate needs a positive finite gas temperature")
    if active_discharge_duration_s is not None and (
        not math.isfinite(active_discharge_duration_s)
        or active_discharge_duration_s <= 0
    ):
        raise ValueError("N(2D) fate needs a positive active duration when supplied")
    if any(not math.isfinite(value) or value < 0 for value in neutral_gas_mol.values()):
        raise ValueError("N(2D) bath inventory must be finite and nonnegative")

    base = {
        "status": "unavailable",
        "source": HERRON_1999_SOURCE,
        "temperature_domain_K": [295.0, 301.0],
        "rate_basis": "evaluated 298 K total removal rates, held constant only near 298 K",
        "state_source_assigned": False,
        "products_injected_into_solver": False,
        "quantitative_hcn_yield_claim_allowed": False,
        "scientific_scope": (
            "Conditional homogeneous-gas first collision of thermalized N(2D); "
            "not nascent-fragment kinetics, product branching, wall/radiative loss, "
            "or a coupled state-resolved HCN yield calculation."
        ),
    }
    if not 295.0 <= gas_temperature_K <= 301.0:
        return {**base, "reason": "outside reviewed near-298 K diagnostic domain"}
    total_mol = sum(neutral_gas_mol.values())
    if total_mol <= 0:
        return {**base, "reason": "empty gas bath"}

    supported = {species for species, *_ in _CHANNELS}
    tolerance_mol = max(total_mol * 1e-12, 1e-18)
    unreviewed = sorted(
        species for species, amount in neutral_gas_mol.items()
        if species not in supported and amount > tolerance_mol
    )
    density_per_mol_cm3 = Avogadro / (gas_volume_L * 1000.0)
    channels = [
        {
            "bath_species": species,
            "kind": kind,
            "rate_coefficient_cm3_molecule_s": coefficient,
            "multiplicative_rate_error_factor_298K": factor,
            "pseudo_first_order_hazard_s_inv": (
                coefficient * neutral_gas_mol.get(species, 0.0) * density_per_mol_cm3
            ),
        }
        for species, coefficient, factor, kind in _CHANNELS
    ]
    total_hazard = sum(item["pseudo_first_order_hazard_s_inv"] for item in channels)
    if unreviewed:
        return {
            **base,
            "reason": "one or more bath partners lack an evaluated N(2D) removal rate",
            "unreviewed_bath_species": unreviewed,
            "unreviewed_bath_amounts_mol": {
                species: float(neutral_gas_mol[species]) for species in unreviewed
            },
            "unreviewed_bath_mole_fraction": (
                sum(neutral_gas_mol[species] for species in unreviewed) / total_mol
            ),
            "reviewed_partial_removal_hazard_s_inv": total_hazard,
            "unreviewed_rate_for_equal_reviewed_hazard_cm3_molecule_s": {
                species: (
                    total_hazard
                    / (float(neutral_gas_mol[species]) * density_per_mol_cm3)
                )
                for species in unreviewed
            },
            "equal_hazard_rate_is_measured_or_bounded": False,
            "branch_fractions_reported": False,
        }
    if total_hazard <= 0:
        return {**base, "reason": "no positive reviewed N(2D) removal hazard"}

    for item in channels:
        hazard = item["pseudo_first_order_hazard_s_inv"]
        factor = item["multiplicative_rate_error_factor_298K"]
        other_low = sum(
            row["pseudo_first_order_hazard_s_inv"]
            / row["multiplicative_rate_error_factor_298K"]
            for row in channels if row is not item
        )
        other_high = sum(
            row["pseudo_first_order_hazard_s_inv"]
            * row["multiplicative_rate_error_factor_298K"]
            for row in channels if row is not item
        )
        lower = hazard / factor
        upper = hazard * factor
        item["conditional_first_collision_fraction"] = hazard / total_hazard
        item["rate_factor_extreme_fraction_range"] = [
            lower / (lower + other_high) if lower + other_high > 0 else 0.0,
            upper / (upper + other_low) if upper + other_low > 0 else 0.0,
        ]

    lifetime_s = 1.0 / total_hazard
    return {
        **base,
        "status": "conditional_first_collision_fate_resolved",
        "bath_complete_for_evaluated_rates": True,
        "branch_fractions_reported": True,
        "total_removal_hazard_s_inv": total_hazard,
        "conditional_thermalized_lifetime_s": lifetime_s,
        "active_duration_to_lifetime_ratio": (
            active_discharge_duration_s / lifetime_s
            if active_discharge_duration_s is not None else None
        ),
        "channels": channels,
        "rate_factor_extremes_are_confidence_intervals": False,
        "product_branching_resolved": False,
    }
