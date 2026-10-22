"""Uncalibrated dry CH4/N2 gas-packet probe for a published DBD geometry.

The feed, geometry, residence, downstream thermometry and average-power basis come
from Kamarinopoulou et al., Science Advances 10 (2024), adl4246. Its HCN
measurements and published HCN mechanism are deliberately not loaded here.
The unknown E/N, electron density and microdischarge exposure remain declared
scenario controls, so this is not yet an experimental yield comparison.
"""
from __future__ import annotations

import argparse
import json
from typing import Any

import numpy as np
from scipy.constants import R

from app.prebiotic.boltzmann import COMMUNITY_CH4_DATASET_ID
from app.prebiotic.coupled_plasma import solve_mixture_plasma_snapshot
from app.prebiotic.heavy_gas import (
    augment_prepared_with_gri,
    gri_ideal_gas_enthalpy_inventory,
    gri_ideal_gas_internal_energy_inventory,
    run_gri_heavy_species_cycle,
)
from app.prebiotic.schema import BoltzmannDischargeProgram


PAPER_DOI = "10.1126/sciadv.adl4246"
REACTOR_VOLUME_L = 0.00157
RESIDENCE_S = 0.47
AVERAGE_DISCHARGE_POWER_W = 28.1
AC_FREQUENCY_HZ = 22_000.0
DOWNSTREAM_PROBE_TEMPERATURE_K = 555.15  # 282 C, 1 cm beyond the reactor in a separate 28.6 W run.
PRESSURE_PA = 101_325.0
AMBIENT_FLOW_REFERENCE_K = 298.15  # Scenario reference; paper does not fix a sccm standard here.


def _gri_inventory(total_feed_mol: float, methane_fraction: float) -> tuple[dict[str, Any], np.ndarray]:
    prepared: dict[str, Any] = {
        "ids": [],
        "index": {},
        "identities": [],
        "transfer": [],
        "stoichiometry": np.zeros((0, 0)),
    }
    augment_prepared_with_gri(prepared)
    gas = np.zeros(len(prepared["ids"]), dtype=float)
    gas[prepared["index"]["CH4"]] = total_feed_mol * methane_fraction
    gas[prepared["index"]["N2"]] = total_feed_mol * (1.0 - methane_fraction)
    return prepared, gas


def run_dry_flow_probe(
    *,
    methane_fraction: float = 0.01,
    nominal_residence_s: float = RESIDENCE_S,
    reduced_electric_field_Td: float = 200.0,
    electron_density_ceiling_m3: float = 1e20,
    discharge_on: bool = True,
    substeps: int = 16,
    source_forcing: str = "impulse_substeps",
    integration_max_step_fraction: float = 0.05,
    flow_reference_temperature_K: float = AMBIENT_FLOW_REFERENCE_K,
    uniform_gas_temperature_K: float = DOWNSTREAM_PROBE_TEMPERATURE_K,
    trace_methane_neutral_pathways: bool = False,
) -> dict[str, Any]:
    if not 0.0 <= methane_fraction <= 1.0:
        raise ValueError("Methane feed fraction must lie in [0, 1]")
    if not 0.05 <= nominal_residence_s <= 5.0:
        raise ValueError("Nominal residence time must lie in [0.05, 5] s")
    if reduced_electric_field_Td <= 0 or electron_density_ceiling_m3 <= 0:
        raise ValueError("Field and electron-density ceiling must be positive")
    if not 1 <= substeps <= 128:
        raise ValueError("Dry-flow substeps must lie in [1, 128]")
    if source_forcing not in {"impulse_substeps", "inventory_limited", "evolving_eedf"}:
        raise ValueError("Unsupported dry-flow source forcing")
    if not 0 < integration_max_step_fraction <= 1:
        raise ValueError("Integration maximum-step fraction must lie in (0, 1]")
    if not 200.0 <= flow_reference_temperature_K <= 350.0:
        raise ValueError("Flow reference temperature must lie in [200, 350] K")
    if not 350.0 <= uniform_gas_temperature_K <= 1000.0:
        raise ValueError("Uniform gas scenario temperature must lie in [350, 1000] K")
    volume_m3 = REACTOR_VOLUME_L * 1e-3
    total_feed_mol = PRESSURE_PA * volume_m3 / (R * uniform_gas_temperature_K)
    nominal_volume_flow_sccm = REACTOR_VOLUME_L * 1000.0 * 60.0 / nominal_residence_s
    reference_flow_mol_s = (
        PRESSURE_PA * nominal_volume_flow_sccm * 1e-6 / 60.0
        / (R * flow_reference_temperature_K)
    )
    hot_ideal_gas_residence_s = total_feed_mol / reference_flow_mol_s
    # A steady-flow packet contains the hot reactor inventory, not the larger
    # reference-temperature inventory. Its exposure time is inventory / molar
    # flow, and its share of average electrical power is P times that time.
    # The paper does not state the sccm reference temperature, so this is an
    # internally consistent scenario calculation, not a matched experiment.
    packet_residence_s = hot_ideal_gas_residence_s
    prepared, inlet = _gri_inventory(total_feed_mol, methane_fraction)
    initial_gas_internal_energy = gri_ideal_gas_internal_energy_inventory(
        inlet, prepared, uniform_gas_temperature_K
    )
    inlet_stream_enthalpy = gri_ideal_gas_enthalpy_inventory(
        inlet, prepared, flow_reference_temperature_K
    )
    deposited_energy_J = (
        AVERAGE_DISCHARGE_POWER_W * packet_residence_s if discharge_on else 0.0
    )
    controls = (
        BoltzmannDischargeProgram.model_validate({
            "collision_set": COMMUNITY_CH4_DATASET_ID,
            "reduced_electric_field_Td": reduced_electric_field_Td,
            "electron_number_density": {
                "value": electron_density_ceiling_m3,
                "unit": "1/m^3",
                "source": "dry-flow sensitivity ceiling; not a measured electron density",
                "method": "bounded by the computed electron-drift power when necessary",
                "basis": "hypothetical",
            },
            "acknowledge_two_term_limit": True,
        })
        if discharge_on else None
    )
    outlet = inlet.copy()
    chemistry: dict[str, Any] = {}
    plasma: dict[str, Any] = {"status": "discharge_off", "source_terms": []}
    allocated_collision_J = 0.0
    source_energy_used_J = 0.0
    active_fragment_channels: set[str] = set()
    source_extents_mol: dict[str, float] = {}
    source_consumption_mol_by_target: dict[str, float] = {}
    methane_neutral_pathways: dict[int, dict[str, Any]] = {}
    methane_pathway_chemical_net_mol = 0.0
    methane_pathway_observed_net_mol = 0.0
    methane_pathway_quadrature_residual_mol = 0.0
    excluded_fragment_channels: dict[str, dict[str, str]] = {}
    plasma_covered_fractions: list[float] = []
    plasma_effective_densities: list[float] = []
    unresolved_n2_outcome_opportunity_mol = 0.0
    unresolved_n2_outcome_threshold_J = 0.0
    conditional_sub14_ground_n_opportunity_mol = 0.0
    n2_excitation_opportunities_mol_by_label: dict[str, float] = {}
    cycle_count = substeps if source_forcing in {"impulse_substeps", "evolving_eedf"} else 1
    step_duration_s = packet_residence_s / cycle_count
    step_energy_J = deposited_energy_J / cycle_count
    for _ in range(cycle_count):
        if controls is not None:
            plasma = solve_mixture_plasma_snapshot(
                controls=controls,
                identities=prepared["identities"],
                gas_amounts_mol=outlet,
                gas_volume_L=REACTOR_VOLUME_L,
                gas_temperature_K=uniform_gas_temperature_K,
                discharge_duration_s=step_duration_s,
                duty_fraction=1.0,
                deposited_energy_J=step_energy_J,
            )
            plasma_covered_fractions.append(float(plasma["covered_gas_fraction"]))
            plasma_effective_densities.append(float(plasma["electron_number_density_m3"]))
            for channel in plasma["source_terms"]:
                if channel["target"] == "N2" and channel.get("role") == "excitation":
                    label = str(channel["product_label"])
                    n2_excitation_opportunities_mol_by_label[label] = (
                        n2_excitation_opportunities_mol_by_label.get(label, 0.0)
                        + float(channel["uncoupled_active_window_opportunity_mol"])
                    )
                if channel["target"] == "N2" and channel.get("product_state_unresolved"):
                    unresolved_n2_outcome_opportunity_mol += float(
                        channel["uncoupled_active_window_opportunity_mol"]
                    )
                    unresolved_n2_outcome_threshold_J += float(
                        channel["threshold_energy_lower_bound_J"]
                    )
                    window = channel.get("ground_atom_energetic_window") or {}
                    fraction = window.get("rate_fraction_below_that_threshold")
                    if fraction is not None:
                        conditional_sub14_ground_n_opportunity_mol += (
                            float(channel["uncoupled_active_window_opportunity_mol"])
                            * float(fraction)
                            * int(window["conditional_minimum_ground_N4S_atoms_per_event_in_window"])
                        )
        outlet, chemistry = run_gri_heavy_species_cycle(
            outlet,
            prepared,
            plasma,
            uniform_gas_temperature_K,
            REACTOR_VOLUME_L,
            step_duration_s,
            step_energy_J,
            source_forcing=(
                "impulse" if source_forcing == "impulse_substeps" else "inventory_limited"
            ),
            include_diagnostic_shadows=False,
            integration_max_step_fraction=integration_max_step_fraction,
            integrated_pathway_species=("CH4",) if trace_methane_neutral_pathways else (),
        )
        if trace_methane_neutral_pathways:
            methane_diagnostic = chemistry["integrated_pathway_diagnostics"]["species"]["CH4"]
            methane_pathway_chemical_net_mol += methane_diagnostic["integrated_chemical_net_mol"]
            methane_pathway_observed_net_mol += methane_diagnostic["observed_chemical_net_mol"]
            methane_pathway_quadrature_residual_mol += methane_diagnostic["quadrature_closure_residual_mol"]
            for item in methane_diagnostic["reaction_contributions"]:
                reaction_index = int(item["reaction_index"])
                if reaction_index not in methane_neutral_pathways:
                    methane_neutral_pathways[reaction_index] = {
                        "reaction_index": reaction_index,
                        "equation": item["equation"],
                        "species_net_mol": 0.0,
                    }
                methane_neutral_pathways[reaction_index]["species_net_mol"] += item["species_net_mol"]
        allocated_collision_J += float(
            chemistry["collision_energy_allocation"]["allocated_collision_energy_J"]
        )
        source_energy_used_J += float(chemistry["source_energy_used_J"])
        active_fragment_channels.update(
            str(item["channel_id"]) for item in chemistry["source_extents"]
        )
        for item in chemistry["source_extents"]:
            channel_id = str(item["channel_id"])
            target = str(item["target"])
            actual_extent_mol = float(
                item.get("actual_integrated_extent_mol", item["accepted_extent_mol"])
            )
            source_extents_mol[channel_id] = (
                source_extents_mol.get(channel_id, 0.0)
                + actual_extent_mol
            )
            source_consumption_mol_by_target[target] = (
                source_consumption_mol_by_target.get(target, 0.0)
                + actual_extent_mol
            )
        for item in chemistry["collision_energy_allocation"]["excluded_fragment_channels"]:
            channel_id = str(item["channel_id"])
            excluded_fragment_channels[channel_id] = {
                "product_label": str(item["product_label"]),
                "reason": str(item["reason"]),
            }
    residual = prepared["atomic"] @ (outlet - inlet)
    final_gas_internal_energy = gri_ideal_gas_internal_energy_inventory(
        outlet, prepared, uniform_gas_temperature_K
    )
    outlet_stream_enthalpy = gri_ideal_gas_enthalpy_inventory(
        outlet, prepared, uniform_gas_temperature_K
    )
    stream_enthalpy_coverage_complete = (
        inlet_stream_enthalpy["all_nonzero_species_covered"]
        and outlet_stream_enthalpy["all_nonzero_species_covered"]
    )
    stream_enthalpy_change_J_per_packet = (
        outlet_stream_enthalpy["covered_enthalpy_J"]
        - inlet_stream_enthalpy["covered_enthalpy_J"]
    )
    stream_enthalpy_rate_W = stream_enthalpy_change_J_per_packet / packet_residence_s
    ideal_gas_internal_energy_change_J = (
        final_gas_internal_energy["covered_internal_energy_J"]
        - initial_gas_internal_energy["covered_internal_energy_J"]
    )
    gas_energy_coverage_complete = (
        initial_gas_internal_energy["all_nonzero_species_covered"]
        and final_gas_internal_energy["all_nonzero_species_covered"]
    )
    maximum_element_residual_mol = float(np.max(np.abs(residual), initial=0.0))
    initial_conserved_amounts = prepared["atomic"] @ inlet
    material_tolerances = np.maximum(np.abs(initial_conserved_amounts) * 1e-12, 1e-12)
    material_closed = bool(np.all(np.abs(residual) <= material_tolerances))
    hcn_mol = float(outlet[prepared["index"]["HCN"]])
    feed_carbon_mol = total_feed_mol * methane_fraction
    outlet_species_mol = {
        species_id: float(outlet[position])
        for position, species_id in enumerate(prepared["ids"])
    }
    outlet_total_mol = sum(outlet_species_mol.values())
    outlet_pressure_Pa = outlet_total_mol * R * uniform_gas_temperature_K / volume_m3
    carbon_atom_counts = prepared["atomic"][prepared["conserved"].index("C")]
    carbon_yield_by_species = {
        species_id: float(carbon_atom_counts[position]) * amount / feed_carbon_mol
        for position, (species_id, amount) in enumerate(outlet_species_mol.items())
        if carbon_atom_counts[position] > 0 and feed_carbon_mol > 0
    }
    inlet_methane_mol = float(inlet[prepared["index"]["CH4"]])
    methane_conversion = (
        (inlet_methane_mol - outlet_species_mol["CH4"]) / inlet_methane_mol
        if inlet_methane_mol > 0 else None
    )
    methane_source_consumption_mol = source_consumption_mol_by_target.get("CH4", 0.0)
    methane_net_neutral_consumption_mol = (
        inlet_methane_mol - outlet_species_mol["CH4"] - methane_source_consumption_mol
    )
    # Paper eqs. S4 and S7 use converted-feed carbon for selectivity and
    # inlet-feed carbon for yield. Keep this all-species model accounting
    # separate from the paper's GC-observed species and deposited coke.
    product_carbon_yield_by_species = {
        species_id: yield_fraction
        for species_id, yield_fraction in carbon_yield_by_species.items()
        if species_id != "CH4"
    }
    product_carbon_selectivity_by_species = (
        {
            species_id: yield_fraction / methane_conversion
            for species_id, yield_fraction in product_carbon_yield_by_species.items()
        }
        if methane_conversion is not None and methane_conversion > 0 else {}
    )
    hydrogen_atom_counts = prepared["atomic"][prepared["conserved"].index("H")]
    consumed_feed_hydrogen_mol = (
        4.0 * inlet_methane_mol * methane_conversion
        if methane_conversion is not None else 0.0
    )
    product_hydrogen_selectivity_by_species = (
        {
            species_id: float(hydrogen_atom_counts[position]) * amount
            / consumed_feed_hydrogen_mol
            for position, (species_id, amount) in enumerate(outlet_species_mol.items())
            if species_id != "CH4" and hydrogen_atom_counts[position] > 0
        }
        if consumed_feed_hydrogen_mol > 0 else {}
    )
    return {
        "status": "dry_flow_uncalibrated_probe_complete",
        "source": {"doi": PAPER_DOI, "measured_hcn_used_as_solver_input": False},
        "conditions": {
            "initial_CH4_mole_fraction": methane_fraction,
            "reactor_idealization": "Lagrangian plug-flow packet; no axial dispersion, radial gradients, or measured streamer occupancy",
            "reactor_volume_L": REACTOR_VOLUME_L,
            "nominal_residence_s": nominal_residence_s,
            "simulated_hot_packet_residence_s": packet_residence_s,
            "measured_steady_downstream_probe_temperature_K": DOWNSTREAM_PROBE_TEMPERATURE_K,
            "probe_distance_above_reactor_cm": 1.0,
            "thermometry_run_power_W": 28.6,
            "thermometry_steady_after_min": 10.0,
            "probe_is_spatially_resolved_reactor_temperature": False,
            "uniform_gas_temperature_K_scenario": uniform_gas_temperature_K,
            "uniform_temperature_equals_downstream_probe_value": (
                uniform_gas_temperature_K == DOWNSTREAM_PROBE_TEMPERATURE_K
            ),
            "pressure_Pa": PRESSURE_PA,
            "average_discharge_power_W": AVERAGE_DISCHARGE_POWER_W if discharge_on else 0.0,
            "nominal_ac_frequency_Hz": AC_FREQUENCY_HZ,
            "nominal_ac_period_s": 1.0 / AC_FREQUENCY_HZ,
            "nominal_ac_cycles_per_hot_packet": packet_residence_s * AC_FREQUENCY_HZ,
            "chemistry_substep_duration_s": step_duration_s,
            "nominal_ac_cycles_per_chemistry_substep": step_duration_s * AC_FREQUENCY_HZ,
            "individual_microdischarge_pulses_resolved": False,
            "power_times_residence_J_per_gas_packet": deposited_energy_J,
            "reduced_electric_field_Td_scenario": reduced_electric_field_Td,
            "electron_density_ceiling_m3_scenario": electron_density_ceiling_m3,
            "microdischarge_exposure_fraction_scenario": 1.0,
            "collision_set": COMMUNITY_CH4_DATASET_ID,
            "discharge_chemistry_substeps": cycle_count,
            "source_forcing": source_forcing,
            "integration_max_step_fraction": integration_max_step_fraction,
        },
        "inlet_total_mol": total_feed_mol,
        "outlet_packet_diagnostics": {
            "scope": "simulated Lagrangian plug-flow packet; not a calibrated continuous-flow outlet assay",
            "outlet_pressure_Pa_from_fixed_volume_inventory": outlet_pressure_Pa,
            "pressure_drift_fraction_from_inlet": (outlet_pressure_Pa - PRESSURE_PA) / PRESSURE_PA,
            "species_mol": outlet_species_mol,
            "species_mole_fraction": {
                species_id: amount / outlet_total_mol
                for species_id, amount in outlet_species_mol.items()
            },
            "methane_conversion_fraction": methane_conversion,
            "methane_electron_fragment_source_consumption_mol": methane_source_consumption_mol,
            "methane_net_neutral_chemistry_consumption_mol": methane_net_neutral_consumption_mol,
            "methane_electron_fragment_source_fraction_of_feed": (
                methane_source_consumption_mol / inlet_methane_mol
                if inlet_methane_mol > 0 else None
            ),
            "methane_net_neutral_chemistry_consumption_fraction_of_feed": (
                methane_net_neutral_consumption_mol / inlet_methane_mol
                if inlet_methane_mol > 0 else None
            ),
            "methane_conversion_partition_interpretation": (
                "Net CH4 conversion equals integrated electron-fragment source "
                "consumption plus net neutral-chemistry consumption. A negative "
                "neutral term means net CH4 regeneration, not a negative reaction rate."
            ),
            "carbon_yield_fraction_by_species": carbon_yield_by_species,
            "carbon_yield_normalization": "carbon atoms in outlet species divided by inlet CH4 carbon atoms",
            "product_carbon_yield_fraction_by_species": product_carbon_yield_by_species,
            "product_carbon_selectivity_fraction_by_species": product_carbon_selectivity_by_species,
            "product_hydrogen_selectivity_fraction_by_species": (
                product_hydrogen_selectivity_by_species
            ),
            "product_selectivity_normalization": (
                "carbon atoms in each modeled non-CH4 outlet species divided by "
                "inlet CH4 carbon atoms consumed; includes model radicals and "
                "omits deposited coke, unlike the paper's GC-based selectivity"
            ),
            "product_hydrogen_selectivity_normalization": (
                "hydrogen atoms in each modeled non-CH4 outlet species divided by "
                "hydrogen atoms in consumed inlet CH4; includes modeled radicals "
                "and omits any unmodeled wall deposits"
            ),
            "measured_product_values_used_as_inputs": False,
        },
        "flow_temperature_basis_diagnostic": {
            "nominal_volume_flow_implied_by_reported_residence_sccm": nominal_volume_flow_sccm,
            "assumed_ambient_flow_reference_K": flow_reference_temperature_K,
            "inferred_feed_molar_flow_mol_s": reference_flow_mol_s,
            "ideal_gas_residence_at_uniform_scenario_temperature_s": hot_ideal_gas_residence_s,
            "hot_to_nominal_residence_ratio": hot_ideal_gas_residence_s / nominal_residence_s,
            "internally_consistent_hot_packet_flow_and_power_basis": (
                abs(total_feed_mol / packet_residence_s - reference_flow_mol_s)
                <= reference_flow_mol_s * 1e-12
                and abs(deposited_energy_J - (
                    AVERAGE_DISCHARGE_POWER_W * packet_residence_s
                    if discharge_on else 0.0
                )) <= 1e-12
            ),
            "comparison_basis_reconciled": False,
            "interpretation": (
                "The packet now uses hot ideal-gas inventory divided by the inferred "
                "reference-flow molar throughput, with the same residence time in "
                "power-times-time energy. The sccm reference temperature and reactor "
                "temperature profile are not established by the accessible paper, so "
                "its yield and energy cost are not yet experimentally matched."
            ),
        },
        "outlet_hcn_mol": hcn_mol,
        "outlet_hcn_per_feed_carbon": hcn_mol / feed_carbon_mol if feed_carbon_mol else None,
        "hcn_nmol_per_declared_discharge_J": (
            hcn_mol * 1e9 / deposited_energy_J if deposited_energy_J else None
        ),
        "plasma": {
            "status": plasma["status"],
            "covered_gas_fraction": plasma.get("covered_gas_fraction"),
            "effective_electron_density_m3": plasma.get("electron_number_density_m3"),
            "density_was_power_capped": (
                plasma.get("drift_power_diagnostic") or {}
            ).get("used_to_scale_source_terms"),
            "reactive_channel_count": plasma.get("reactive_channel_count"),
            "mean_electron_energy_eV": plasma.get("mean_electron_energy_eV"),
            "snapshot_count": len(plasma_covered_fractions),
            "covered_gas_fraction_range": (
                [min(plasma_covered_fractions), max(plasma_covered_fractions)]
                if plasma_covered_fractions else None
            ),
            "effective_electron_density_m3_range": (
                [min(plasma_effective_densities), max(plasma_effective_densities)]
                if plasma_effective_densities else None
            ),
            "unresolved_n2_total_outcome_opportunity_mol": unresolved_n2_outcome_opportunity_mol,
            "unresolved_n2_total_outcome_threshold_energy_J": unresolved_n2_outcome_threshold_J,
            "conditional_sub14_526eV_ground_n_opportunity_mol": conditional_sub14_ground_n_opportunity_mol,
            "n2_excitation_opportunity_mol_by_product_label": (
                n2_excitation_opportunities_mol_by_label
            ),
            "n2_excitation_opportunity_is_propagated": False,
            "n2_excitation_opportunity_interpretation": (
                "Read-only uncoupled electron-impact opportunities at the hypothetical "
                "power-capped density. Spectroscopic bands are not guaranteed unique "
                "metastable states; these are not accepted excited-state inventories, "
                "and their threshold-energy demands overlap the plasma allocation."
            ),
            "unresolved_n2_opportunity_is_propagated": False,
            "unresolved_n2_opportunity_interpretation": (
                "Integral of uncoupled total N2 outcome rates at the power-capped "
                "electron density, before target depletion, product-state branching, "
                "or independent exposure validation. It is not an actual dissociation "
                "count or N-atom yield; the threshold estimate overlaps existing "
                "electronic-loss channels and is not added to the energy ledger."
            ),
            "conditional_ground_n_opportunity_interpretation": (
                "If the tabulated total neutral N2 outcome is valid below 14.526 eV "
                "for ground-state N2, energetics require at least one N(4S) in "
                "each such event. This integrates conditional opportunities only: "
                "near-threshold fragments were not directly state-measured, actual "
                "event counts and streamer exposure are unknown, the partner state "
                "is unresolved, and no N source is propagated."
            ),
            "conditional_n2d_fate_status": (
                (plasma.get("n2d_fate_diagnostic") or {}).get("status")
            ),
            "conditional_n2d_fate_unavailable_reason": (
                (plasma.get("n2d_fate_diagnostic") or {}).get("reason")
            ),
        },
        "energy_ledger": {
            "scope": "electrical collision-source upper bound; not a complete thermal balance",
            "complete_thermal_energy_balance_closed": False,
            "gas_thermochemistry_basis": initial_gas_internal_energy["basis"],
            "gas_internal_energy_inventory_coverage_complete": gas_energy_coverage_complete,
            "initial_ideal_gas_internal_energy_J": initial_gas_internal_energy[
                "covered_internal_energy_J"
            ],
            "final_ideal_gas_internal_energy_J": final_gas_internal_energy[
                "covered_internal_energy_J"
            ],
            "ideal_gas_internal_energy_change_J": ideal_gas_internal_energy_change_J,
            "conditional_unmeasured_non_electrical_energy_into_gas_J": (
                ideal_gas_internal_energy_change_J - deposited_energy_J
                if gas_energy_coverage_complete else None
            ),
            "conditional_exchange_interpretation": (
                "The first-law remainder assumes all declared electrical energy enters "
                "the fixed-volume modeled gas. It combines unmeasured heat, radiation, "
                "wall coupling, and omitted plasma energy storage; it is inferred, "
                "not a measured or closed thermal balance."
            ),
            "steady_flow_enthalpy_scope": (
                "Gas-only inlet-to-outlet stream enthalpy at assumed ambient "
                "feed temperature and modeled uniform gas scenario temperature, "
                "which is not a measured in-reactor or outlet profile; not a closed "
                "reactor heat balance"
            ),
            "steady_flow_enthalpy_coverage_complete": stream_enthalpy_coverage_complete,
            "inlet_stream_enthalpy_J_per_packet": inlet_stream_enthalpy["covered_enthalpy_J"],
            "outlet_stream_enthalpy_J_per_packet": outlet_stream_enthalpy["covered_enthalpy_J"],
            "stream_enthalpy_change_J_per_packet": stream_enthalpy_change_J_per_packet,
            "stream_enthalpy_change_rate_W": stream_enthalpy_rate_W,
            "conditional_electrical_power_minus_stream_enthalpy_rate_W": (
                (AVERAGE_DISCHARGE_POWER_W if discharge_on else 0.0) - stream_enthalpy_rate_W
                if stream_enthalpy_coverage_complete else None
            ),
            "steady_flow_enthalpy_interpretation": (
                "The residual is not measured wall loss: the flow reference "
                "temperature, internal temperature profile, outlet chemistry, "
                "radiation, electrode coupling, and other energy paths are not "
                "independently resolved."
            ),
            "declared_discharge_J": deposited_energy_J,
            "allocated_collision_J": allocated_collision_J,
            "resolved_fragment_J": source_energy_used_J,
            "closed_upper_bound": (
                source_energy_used_J <= deposited_energy_J * (1 + 1e-9)
                and source_energy_used_J <= allocated_collision_J * (1 + 1e-9)
                and allocated_collision_J <= deposited_energy_J * (1 + 1e-9)
            ),
        },
        "material_ledger": {
            "maximum_element_residual_mol": maximum_element_residual_mol,
            "tolerance_mol_by_conserved_quantity": {
                key: float(tolerance)
                for key, tolerance in zip(prepared["conserved"], material_tolerances)
            },
            "residual_mol_by_conserved_quantity": {
                key: float(value)
                for key, value in zip(prepared["conserved"], residual)
            },
            "closed": material_closed,
        },
        "heavy_chemistry": {
            "status": chemistry["status"],
            "active_fragment_channel_count": len(active_fragment_channels),
            "active_fragment_channel_ids": sorted(active_fragment_channels),
            "source_extents_mol_by_channel": source_extents_mol,
            "source_consumption_mol_by_target": source_consumption_mol_by_target,
            "methane_neutral_pathways": (
                {
                    "method": "read-only integrated GRI reaction rates over each packet substep",
                    "integrated_chemical_net_mol": methane_pathway_chemical_net_mol,
                    "observed_chemical_net_mol": methane_pathway_observed_net_mol,
                    "quadrature_closure_residual_mol": methane_pathway_quadrature_residual_mol,
                    "reaction_contributions": sorted(
                        methane_neutral_pathways.values(),
                        key=lambda item: abs(item["species_net_mol"]),
                        reverse=True,
                    ),
                    "not_a_validated_mechanism": True,
                }
                if trace_methane_neutral_pathways else None
            ),
            "excluded_fragment_channels_by_id": excluded_fragment_channels,
            "solver": chemistry["solver"],
            "hcn_pathways_at_outlet": chemistry["final_instantaneous_pathway_diagnostics"].get("HCN"),
        },
        "scientific_firewall": {
            "published_hcn_yield_used_to_choose_source_conditions": False,
            "published_hcn_mechanism_used": False,
            "flow_reactor_spatial_or_microdischarge_model_solved": False,
            "impulse_substep_result_is_timestep_converged": False,
            "inventory_limited_result_is_validated": False,
            "evolving_eedf_result_is_timestep_converged": False,
            "experimental_magnitude_validated": False,
            "wet_miller_urey_reproduction_claimed": False,
            "flow_temperature_residence_basis_reconciled": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--methane-fraction", type=float, default=0.01)
    parser.add_argument("--field-td", type=float, default=200.0)
    parser.add_argument("--electron-density-ceiling-m3", type=float, default=1e20)
    parser.add_argument("--off", action="store_true")
    parser.add_argument("--substeps", type=int, default=16)
    parser.add_argument(
        "--source-forcing", choices=("impulse_substeps", "inventory_limited", "evolving_eedf"),
        default="impulse_substeps",
    )
    parser.add_argument("--integration-max-step-fraction", type=float, default=0.05)
    parser.add_argument("--flow-reference-temperature-k", type=float, default=AMBIENT_FLOW_REFERENCE_K)
    parser.add_argument("--uniform-gas-temperature-k", type=float, default=DOWNSTREAM_PROBE_TEMPERATURE_K)
    arguments = parser.parse_args()
    print(json.dumps(run_dry_flow_probe(
        methane_fraction=arguments.methane_fraction,
        reduced_electric_field_Td=arguments.field_td,
        electron_density_ceiling_m3=arguments.electron_density_ceiling_m3,
        discharge_on=not arguments.off,
        substeps=arguments.substeps,
        source_forcing=arguments.source_forcing,
        integration_max_step_fraction=arguments.integration_max_step_fraction,
        flow_reference_temperature_K=arguments.flow_reference_temperature_k,
        uniform_gas_temperature_K=arguments.uniform_gas_temperature_k,
    ), indent=2))


if __name__ == "__main__":
    main()
