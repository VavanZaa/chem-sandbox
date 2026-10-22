"""Read-only near-room-temperature fate audit for electron-produced O(1D).

These are competing elementary removal channels, not an HCN recipe. No
products or extents from this module enter the predictive inventories yet.
The deliberately narrow temperature and bath-species domain prevents a
published atmospheric-rate table from silently becoming a generic spark
chemistry model.
"""

from __future__ import annotations

import math
from typing import Any, Mapping

from scipy.constants import Avogadro, c, elementary_charge, h


JPL_O1D_SOURCE = (
    "NASA/JPL Evaluation 17, Publication 10-6 (2011), Table 1-1 and notes A2-A10; "
    "https://jpldataeval.jpl.nasa.gov/pdf/JPL%2010-6%20Final%2015June2011_20110620.pdf"
)
DAVIDSON_CO_SOURCE = (
    "Davidson, Schiff, Brown, and Howard, J. Chem. Phys. 69 (1978), "
    "DOI 10.1063/1.436657; "
    "https://ntrs.nasa.gov/citations/19780065438"
)
STREKOWSKI_HCN_SOURCE = (
    "Strekowski, Nicovich, and Wine, ChemPhysChem 11 (2010), 3942-3955, "
    "DOI 10.1002/cphc.201000550; measured total O(1D)+HCN removal and "
    "O(3P) quenching yield over 193-430 K"
)
VRANCKX_HYDROGEN_SOURCE = (
    "Vranckx, Peeters, and Carl, Phys. Chem. Chem. Phys. 12 (2010), 9213-9221, "
    "DOI 10.1039/b923959f; measured O(1D)+H2O/H2 total removal and O(3P) "
    "quenching yields over 227-453 K"
)
NIST_O1D_LEVEL_CM_INV = 15867.862
NIST_O1D_LEVEL_SOURCE = (
    "NIST neutral oxygen energy levels, O I 2s2 2p4 1D2 at 15867.862 cm-1 "
    "above the 3P2 ground level; "
    "https://physics.nist.gov/PhysRefData/Handbook/Tables/oxygentable5_a.htm"
)
NIST_FORMATION_ENTHALPY_KJ_MOL = {
    "CO2": -393.51,
    "CO": -110.53,
    "O": 249.18,
}
NIST_FORMATION_ENTHALPY_SOURCE = (
    "NIST Chemistry WebBook CODATA gas-phase formation enthalpies at 298.15 K: "
    "CO2 https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Mask=2781&Units=SI ; "
    "CO https://webbook.nist.gov/cgi/cbook.cgi?ID=C630080&Mask=FFFFFF ; "
    "O https://webbook.nist.gov/cgi/cbook.cgi?ID=C17778802&Mask=29"
)

# k(298 K) in cm3 molecule-1 s-1. These three channels physically
# deactivate O(1D) to O(3P). H2O/H2 have measured mostly-reactive branches
# below. This audit uses no historical Miller-Urey products, yields, or recipes.
_O1D_CHANNELS = (
    ("N2", 3.10e-11, "physical_quench", {"O": 1, "N2": 1}),
    ("CO2", 1.10e-10, "physical_quench", {"O": 1, "CO2": 1}),
    ("O2", 3.95e-11, "physical_quench", {"O": 1, "O2": 1}),
)

# JPL Evaluation 17, Table 1-1, notes A8, A9 and A11. N2O has two
# explicitly evaluated product branches. CH4 has a measured total removal
# rate, but its 20% methoxy/hydroxymethyl branch is not isomer-resolved here.
_N2O_TO_N2_O2_RATE = 4.95e-11
_N2O_TO_TWO_NO_RATE = 7.75e-11
_N2O_TOTAL_RATE = _N2O_TO_N2_O2_RATE + _N2O_TO_TWO_NO_RATE
_NH3_TO_OH_NH2_RATE = 2.50e-10
_CH4_TOTAL_REMOVAL_RATE = 1.75e-10


def _co_removal_rate_cm3_molecule_s(
    temperature_K: float, *, preexponential: float = 4.7e-11,
    activation_cal_mol: float = 126.0,
) -> float:
    """Measured total O(1D) deactivation by CO, not a product-branch rate.

    Davidson et al. measured (4.7 +/- 0.9)e-11 exp[(126 +/- 33) cal/mol / RT]
    from 113 to 333 K. This audit remains narrower because the other bath
    coefficients are held at their 298 K values.
    """

    gas_constant_cal_mol_K = 1.98720425864083
    return preexponential * math.exp(
        activation_cal_mol / (gas_constant_cal_mol_K * temperature_K)
    )


def _hcn_removal_and_quench_fraction(temperature_K: float) -> tuple[float, float]:
    """Measured O(1D)+HCN total loss and O(3P) yield, not a product recipe."""

    return (
        1.08e-10 * math.exp(105.0 / temperature_K),
        0.150 * math.exp(200.0 / temperature_K),
    )


def o1d_fate_diagnostic(
    neutral_gas_mol: Mapping[str, float],
    *,
    gas_volume_L: float,
    gas_temperature_K: float,
    active_discharge_duration_s: float | None = None,
) -> dict[str, Any]:
    """Audit first-collision fate in a homogeneous gas; never inject products.

    The short-lifetime approximation can only be assessed after *all*
    materially present bath partners are covered. Otherwise branch fractions
    are omitted, rather than renormalized over a partial mixture.
    """

    if not math.isfinite(gas_volume_L) or gas_volume_L <= 0:
        raise ValueError("O(1D) fate needs a positive finite gas volume")
    if not math.isfinite(gas_temperature_K) or gas_temperature_K <= 0:
        raise ValueError("O(1D) fate needs a positive finite gas temperature")
    if active_discharge_duration_s is not None and (
        not math.isfinite(active_discharge_duration_s)
        or active_discharge_duration_s <= 0
    ):
        raise ValueError("O(1D) fate needs a positive active duration when supplied")
    if any(not math.isfinite(value) or value < 0 for value in neutral_gas_mol.values()):
        raise ValueError("O(1D) bath inventory must be finite and nonnegative")

    total_mol = sum(neutral_gas_mol.values())
    base = {
        "status": "unavailable",
        "source": JPL_O1D_SOURCE,
        "temperature_domain_K": [295.0, 301.0],
        "rate_basis": "published 298 K bimolecular removal coefficients, held constant only in this narrow audit",
        "scientific_scope": (
            "Conditional first-collision fate for thermalized O(1D) at near-298 K; "
            "nascent electron-dissociation fragments may be nonthermal. Diagnostic "
            "only: not a full excited-state network, heat balance, or HCN yield."
        ),
        "products_injected_into_solver": False,
    }
    if not 295.0 <= gas_temperature_K <= 301.0:
        return {**base, "reason": "outside reviewed near-298 K diagnostic domain"}
    if total_mol <= 0:
        return {**base, "reason": "empty gas bath"}

    supported = {item[0] for item in _O1D_CHANNELS} | {
        "H2O", "H2", "CO", "HCN", "N2O", "NH3", "CH4",
    }
    tolerance_mol = max(total_mol * 1e-12, 1e-18)
    unreviewed = sorted(
        species for species, amount in neutral_gas_mol.items()
        if species not in supported and amount > tolerance_mol
    )
    if unreviewed:
        molecule_density_per_mol_cm3 = Avogadro / (gas_volume_L * 1000.0)
        known_coefficients = {
            species: coefficient for species, coefficient, _kind, _products in _O1D_CHANNELS
        }
        known_coefficients.update({
            "H2O": 1.70e-10 * math.exp(36.0 / gas_temperature_K),
            "H2": 1.35e-10,
            "CO": _co_removal_rate_cm3_molecule_s(gas_temperature_K),
            "HCN": _hcn_removal_and_quench_fraction(gas_temperature_K)[0],
            "N2O": _N2O_TOTAL_RATE,
            "NH3": _NH3_TO_OH_NH2_RATE,
            "CH4": _CH4_TOTAL_REMOVAL_RATE,
        })
        known_hazard = molecule_density_per_mol_cm3 * sum(
            coefficient * neutral_gas_mol.get(species, 0.0)
            for species, coefficient in known_coefficients.items()
        )
        unreviewed_amount = sum(neutral_gas_mol[species] for species in unreviewed)
        # These are deliberately arbitrary stress coefficients, not measured
        # rate limits. The endpoint remains unavailable at either value.
        rate_stress = (1e-9, 1e-8)
        return {
            **base,
            "reason": "one or more bath partners lack a reviewed O(1D) fate channel",
            "unreviewed_bath_species": unreviewed,
            "unreviewed_bath_amounts_mol": {
                species: float(neutral_gas_mol[species]) for species in unreviewed
            },
            "unreviewed_bath_mole_fraction": (
                unreviewed_amount / total_mol
            ),
            "reviewed_total_removal_hazard_s_inv": known_hazard,
            "unreviewed_rate_stress": [
                {
                    "hypothetical_rate_cm3_molecule_s": rate,
                    "hypothetical_unknown_hazard_fraction": (
                        rate * unreviewed_amount * molecule_density_per_mol_cm3
                        / (known_hazard + rate * unreviewed_amount * molecule_density_per_mol_cm3)
                    ) if known_hazard + rate * unreviewed_amount * molecule_density_per_mol_cm3 > 0 else None,
                }
                for rate in rate_stress
            ],
            "unreviewed_rate_stress_scope": (
                "The two common rates are hypothetical sensitivity points, not measured "
                "values, physical upper bounds, uncertainty limits, or permission to "
                "renormalize the reviewed branches. No product is injected."
            ),
        }

    molecule_density_per_mol_cm3 = Avogadro / (gas_volume_L * 1000.0)
    channels = [
        {
            "bath_species": species,
            "kind": kind,
            "products": products,
            "rate_coefficient_cm3_molecule_s": coefficient,
            "pseudo_first_order_hazard_s_inv": (
                coefficient * neutral_gas_mol.get(species, 0.0)
                * molecule_density_per_mol_cm3
            ),
        }
        for species, coefficient, kind, products in _O1D_CHANNELS
    ]
    for species, coefficient, products, quench_yield, quench_spread in (
        (
            "H2O", 1.70e-10 * math.exp(36.0 / gas_temperature_K),
            {"OH": 2}, 0.010, 0.003,
        ),
        ("H2", 1.35e-10, {"OH": 1, "H": 1}, 0.007, 0.007),
    ):
        channels.append({
            "bath_species": species,
            "kind": "reactive_with_measured_minor_quench",
            "products": products,
            "minor_quench_products": {"O": 1, species: 1},
            "source": VRANCKX_HYDROGEN_SOURCE,
            "rate_coefficient_cm3_molecule_s": coefficient,
            "ground_O_yield_per_O1D_removed": quench_yield,
            "ground_O_yield_reported_plus_minus": quench_spread,
            "pseudo_first_order_hazard_s_inv": (
                coefficient * neutral_gas_mol.get(species, 0.0)
                * molecule_density_per_mol_cm3
            ),
        })
    for species, coefficient, products in (
        ("N2O", _N2O_TOTAL_RATE, {"N2": 1, "O2": 1}),
        ("NH3", _NH3_TO_OH_NH2_RATE, {"OH": 1, "NH2": 1}),
    ):
        item = {
            "bath_species": species,
            "kind": "reactive_with_reviewed_products",
            "products": products,
            "source": JPL_O1D_SOURCE,
            "rate_coefficient_cm3_molecule_s": coefficient,
            "pseudo_first_order_hazard_s_inv": (
                coefficient * neutral_gas_mol.get(species, 0.0)
                * molecule_density_per_mol_cm3
            ),
        }
        if species == "N2O":
            item.update({
                "branch_fraction_to_listed_products": (
                    _N2O_TO_N2_O2_RATE / _N2O_TOTAL_RATE
                ),
                "alternate_products": {"NO": 2},
                "alternate_branch_fraction": (
                    _N2O_TO_TWO_NO_RATE / _N2O_TOTAL_RATE
                ),
            })
        channels.append(item)
    ch4_hazard = (
        _CH4_TOTAL_REMOVAL_RATE * neutral_gas_mol.get("CH4", 0.0)
        * molecule_density_per_mol_cm3
    )
    if ch4_hazard > 0:
        channels.append({
            "bath_species": "CH4",
            "kind": "measured_total_with_unresolved_isomer_branch",
            "products": None,
            "source": JPL_O1D_SOURCE,
            "rate_coefficient_cm3_molecule_s": _CH4_TOTAL_REMOVAL_RATE,
            "pseudo_first_order_hazard_s_inv": ch4_hazard,
            "unresolved_product_scope": (
                "JPL's methoxy-or-hydroxymethyl plus H branch is not split "
                "into the two isomers. Total removal is usable for the hazard "
                "audit, but no CH4 product inventory is projected."
            ),
        })
    co_rate = _co_removal_rate_cm3_molecule_s(gas_temperature_K)
    co_hazard = (
        co_rate
        * neutral_gas_mol.get("CO", 0.0)
        * molecule_density_per_mol_cm3
    )
    if co_hazard > 0:
        channels.append({
            "bath_species": "CO",
            "kind": "unresolved_products",
            "products": None,
            "source": DAVIDSON_CO_SOURCE,
            "rate_coefficient_cm3_molecule_s": co_rate,
            "pseudo_first_order_hazard_s_inv": co_hazard,
        })
    hcn_rate, hcn_ground_o_yield = _hcn_removal_and_quench_fraction(gas_temperature_K)
    hcn_hazard = (
        hcn_rate * neutral_gas_mol.get("HCN", 0.0)
        * molecule_density_per_mol_cm3
    )
    if hcn_hazard > 0:
        channels.append({
            "bath_species": "HCN",
            "kind": "measured_quench_with_unresolved_reactive_products",
            "products": None,
            "source": STREKOWSKI_HCN_SOURCE,
            "rate_coefficient_cm3_molecule_s": hcn_rate,
            "ground_O_yield_per_O1D_removed": hcn_ground_o_yield,
            "pseudo_first_order_hazard_s_inv": hcn_hazard,
        })
    total_hazard = sum(item["pseudo_first_order_hazard_s_inv"] for item in channels)
    if total_hazard <= 0:
        return {**base, "reason": "no covered quencher is present"}
    for item in channels:
        item["first_collision_fraction"] = (
            item["pseudo_first_order_hazard_s_inv"] / total_hazard
        )
    physical_fraction = sum(
        item["first_collision_fraction"]
        for item in channels if item["kind"] == "physical_quench"
    )
    physical_fraction += sum(
        item["first_collision_fraction"] * item["ground_O_yield_per_O1D_removed"]
        for item in channels if item["kind"] == "reactive_with_measured_minor_quench"
    )
    reactive_fraction = sum(
        item["first_collision_fraction"] * (1 - item["ground_O_yield_per_O1D_removed"])
        for item in channels if item["kind"] == "reactive_with_measured_minor_quench"
    )
    other_reactive_fraction = sum(
        item["first_collision_fraction"] for item in channels
        if item["kind"] == "reactive_with_reviewed_products"
    )
    hcn_collision_fraction = hcn_hazard / total_hazard
    hcn_ground_o_fraction = hcn_collision_fraction * hcn_ground_o_yield
    hcn_unresolved_reactive_fraction = (
        hcn_collision_fraction * (1.0 - hcn_ground_o_yield)
    )
    physical_fraction += hcn_ground_o_fraction
    lifetime_s = 1.0 / total_hazard
    result = {
        **base,
        "status": (
            "conditional_fate_partially_resolved"
            if co_hazard > 0 or hcn_hazard > 0 or ch4_hazard > 0
            else "conditional_fate_resolved"
        ),
        "gas_temperature_K": gas_temperature_K,
        "gas_volume_L": gas_volume_L,
        "channels": channels,
        "total_removal_hazard_s_inv": total_hazard,
        "mean_first_collision_lifetime_s": lifetime_s,
        "lifetime_to_active_discharge_ratio": (
            lifetime_s / active_discharge_duration_s
            if active_discharge_duration_s is not None else None
        ),
        "unreviewed_bath_species": [],
    }
    if hcn_hazard > 0 or co_hazard > 0 or ch4_hazard > 0:
        result.update({
            "reason": "Measured removal includes unresolved product branches",
            "known_physical_ground_O_fraction": physical_fraction,
            "known_reactive_OH_H_fraction": reactive_fraction,
            "known_reactive_N2O_NH3_fraction": other_reactive_fraction,
            "unresolved_HCN_reactive_fraction": hcn_unresolved_reactive_fraction,
            "unresolved_CH4_fate_fraction": ch4_hazard / total_hazard,
            "hcn_ground_O_quench_fraction_of_all_collisions": hcn_ground_o_fraction,
            "hcn_branch_scope": (
                "The measured O(3P) yield resolves physical quenching; remaining "
                "O(1D)+HCN reactive products are not assigned. The published "
                "298 K H-atom yield includes a fast secondary contribution and "
                "must not be used as a primary H-product branching fraction."
            ),
        })
    if co_hazard > 0:
        # The measured disappearance rate does not determine whether CO
        # physically deactivates O(1D) or stabilizes a CO2 adduct here.
        co_lower_rate = _co_removal_rate_cm3_molecule_s(
            gas_temperature_K, preexponential=3.8e-11,
            activation_cal_mol=93.0,
        )
        co_upper_rate = _co_removal_rate_cm3_molecule_s(
            gas_temperature_K, preexponential=5.6e-11,
            activation_cal_mol=159.0,
        )
        co_density = neutral_gas_mol["CO"] * molecule_density_per_mol_cm3
        non_co_hazard = total_hazard - co_hazard
        result.update({
            "unresolved_CO_fate_fraction": co_hazard / total_hazard,
            "co_rate_parameter_range_cm3_molecule_s": [co_lower_rate, co_upper_rate],
            "lifetime_range_from_CO_rate_only_s": [
                1.0 / (non_co_hazard + co_upper_rate * co_density),
                1.0 / (non_co_hazard + co_lower_rate * co_density),
            ],
            "co_rate_range_scope": (
                "Extremes of the published A and exponent +/- values, not a confidence "
                "interval or total uncertainty; other channel and thermalization "
                "uncertainties are excluded."
            ),
        })
    if co_hazard <= 0 and hcn_hazard <= 0 and ch4_hazard <= 0:
        result.update({
            "physical_ground_O_fraction": physical_fraction,
            "reactive_OH_H_fraction": reactive_fraction,
            "reactive_N2O_NH3_fraction": other_reactive_fraction,
        })
    return result


def project_thermalized_co2_o1d_events(
    *, events_mol: float, neutral_gas_mol: Mapping[str, float],
    gas_volume_L: float, gas_temperature_K: float,
) -> dict[str, Any]:
    """Project a small CO2 -> CO(X)+O(1D) event batch through reviewed fate.

    This is an isolated conditional source kernel, not an apparatus source.
    It refuses incomplete bath/product coverage or a batch large enough to
    change any consumed reactant by more than 1%, so frozen first-collision
    probabilities are not silently used through substantial depletion.
    """
    if not math.isfinite(events_mol) or events_mol < 0:
        raise ValueError("CO2 O(1D) event amount must be finite and nonnegative")
    fate = o1d_fate_diagnostic(
        neutral_gas_mol, gas_volume_L=gas_volume_L,
        gas_temperature_K=gas_temperature_K,
    )
    if fate["status"] != "conditional_fate_resolved":
        raise ValueError("CO2 O(1D) source needs fully resolved near-room-temperature fate")
    if events_mol > 0.01 * neutral_gas_mol.get("CO2", 0.0):
        raise ValueError("CO2 O(1D) event batch exceeds the 1% frozen-bath CO2 limit")

    probabilities = {
        item["bath_species"]: float(item["first_collision_fraction"])
        for item in fate["channels"]
    }
    water_reactive = events_mol * probabilities["H2O"] * (1 - 0.010)
    hydrogen_reactive = events_mol * probabilities["H2"] * (1 - 0.007)
    n2o_reactive = events_mol * probabilities["N2O"]
    ammonia_reactive = events_mol * probabilities["NH3"]
    for species, consumed in (
        ("H2O", water_reactive), ("H2", hydrogen_reactive),
        ("N2O", n2o_reactive), ("NH3", ammonia_reactive),
    ):
        if consumed > 0.01 * neutral_gas_mol.get(species, 0.0):
            raise ValueError(
                f"CO2 O(1D) event batch exceeds the 1% frozen-bath {species} limit"
            )
    ground_oxygen = events_mol * float(fate["physical_ground_O_fraction"])
    n2o_to_n2_o2 = n2o_reactive * _N2O_TO_N2_O2_RATE / _N2O_TOTAL_RATE
    n2o_to_two_no = n2o_reactive * _N2O_TO_TWO_NO_RATE / _N2O_TOTAL_RATE
    delta = {
        "CO2": -events_mol,
        "CO": events_mol,
        "O": ground_oxygen,
        "H2O": -water_reactive,
        "H2": -hydrogen_reactive,
        "OH": 2 * water_reactive + hydrogen_reactive + ammonia_reactive,
        "H": hydrogen_reactive,
        "N2O": -n2o_reactive,
        "N2": n2o_to_n2_o2,
        "O2": n2o_to_n2_o2,
        "NO": 2 * n2o_to_two_no,
        "NH3": -ammonia_reactive,
        "NH2": ammonia_reactive,
    }
    atom_residual = {
        "C": delta["CO2"] + delta["CO"],
        "H": (
            2 * delta["H2O"] + 2 * delta["H2"] + delta["OH"]
            + delta["H"] + 3 * delta["NH3"] + 2 * delta["NH2"]
        ),
        "O": (
            2 * delta["CO2"] + delta["CO"] + delta["O"] + delta["H2O"]
            + delta["OH"] + delta["N2O"] + 2 * delta["O2"] + delta["NO"]
        ),
        "N": (
            2 * delta["N2O"] + 2 * delta["N2"] + delta["NO"]
            + delta["NH3"] + delta["NH2"]
        ),
    }
    atom_tolerance = max(events_mol * 1e-12, 1e-24)
    if any(abs(value) > atom_tolerance for value in atom_residual.values()):
        raise ValueError("CO2 O(1D) projection failed atom closure")
    excitation_energy_J = events_mol * Avogadro * h * c * NIST_O1D_LEVEL_CM_INV * 100.0
    thermochemical_requirement_J_mol = (
        NIST_FORMATION_ENTHALPY_KJ_MOL["CO"]
        + NIST_FORMATION_ENTHALPY_KJ_MOL["O"]
        - NIST_FORMATION_ENTHALPY_KJ_MOL["CO2"]
    ) * 1000.0 + Avogadro * h * c * NIST_O1D_LEVEL_CM_INV * 100.0
    return {
        "status": "isolated_thermalized_o1d_source_projection_not_apparatus_yield",
        "events_mol": events_mol,
        "signed_gas_change_mol": delta,
        "atom_residual_mol": atom_residual,
        "atom_ledger_closed": True,
        "charge_ledger_closed": True,
        "first_collision_lifetime_s": fate["mean_first_collision_lifetime_s"],
        "first_collision_ground_O_fraction": fate["physical_ground_O_fraction"],
        "first_collision_OH_H_fraction": fate["reactive_OH_H_fraction"],
        "electronic_excitation_energy_initially_carried_J": excitation_energy_J,
        "electronic_excitation_energy_level_source": NIST_O1D_LEVEL_SOURCE,
        "co2_to_co_plus_o1d_standard_enthalpy_requirement_eV_per_event": (
            thermochemical_requirement_J_mol / (Avogadro * elementary_charge)
        ),
        "co2_to_co_plus_o1d_standard_enthalpy_requirement_J_for_batch": (
            events_mol * thermochemical_requirement_J_mol
        ),
        "standard_enthalpy_source": NIST_FORMATION_ENTHALPY_SOURCE,
        "electron_collision_energy_allocation_checked": False,
        "excitation_energy_booked_as_thermal_heat": False,
        "complete_thermal_energy_balance_closed": False,
        "products_injected_into_apparatus": False,
        "science_boundary": (
            "A frozen-bath, thermalized O(1D) first-collision projection only. "
            "The actual electron event extent, nascent translational distribution, "
            "excitation/chemical heat partition and later bath chemistry are not solved."
        ),
    }
