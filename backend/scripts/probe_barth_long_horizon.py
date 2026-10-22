"""Feed-only wet-flask horizon stress test; never loads the HCN assay.

This retains the published initial feed but keeps spark exposure, contact
coefficient, and other probe conditions explicitly hypothetical. It tests
whether the evolving numerical apparatus can reach the reported duration;
its endpoint is not an experimentally validated product yield.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from typing import Any

from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.schema import PrebioticApparatusRequest
from scripts.audit_barth_thermal_capacity import audit_barth_water_thermal_capacity
from scripts.audit_ch_n2_association_kinetics import audit_ch_n2_association
from scripts.audit_ch_water_competition import audit_ch_water_competition
from scripts.audit_hcn_association_kinetics import audit_hcn_association_initial_feed
from scripts.audit_source_thermochemistry import audit_source_thermochemistry
from scripts.probe_barth_feed_pathway import FEED_CARD, make_probe_payload


def _aggregate_species_pathways(
    heavy_steps: list[dict[str, Any]], species: str, mechanism_label: str = "GRI",
) -> dict[str, Any]:
    reaction_totals: dict[int, dict[str, Any]] = {}
    gross_formation = 0.0
    gross_consumption = 0.0
    observed_net = 0.0
    quadrature_net = 0.0
    max_step_residual = 0.0
    method_counts: Counter[str] = Counter()
    for step in heavy_steps:
        trace = step.get("integrated_pathway_diagnostics") or {}
        audit = trace.get("species", {}).get(species)
        if audit is None:
            raise ValueError(f"Full-horizon {species} pathway trace lacks an apparatus substep")
        method = trace.get("method")
        if not isinstance(method, str) or not method:
            raise ValueError(f"Full-horizon {species} pathway trace lacks a quadrature method")
        method_counts[method] += 1
        gross_formation += float(audit["integrated_formation_mol"])
        gross_consumption += float(audit["integrated_consumption_mol"])
        observed_net += float(audit["observed_chemical_net_mol"])
        quadrature_net += float(audit["integrated_chemical_net_mol"])
        max_step_residual = max(
            max_step_residual, abs(float(audit["quadrature_closure_residual_mol"])),
        )
        for item in audit["reaction_contributions"]:
            entry = reaction_totals.setdefault(int(item["reaction_index"]), {
                "reaction_index": int(item["reaction_index"]),
                "equation": item["equation"],
                "formation_mol": 0.0,
                "consumption_mol": 0.0,
            })
            amount = float(item["species_net_mol"])
            if amount > 0:
                entry["formation_mol"] += amount
            else:
                entry["consumption_mol"] -= amount
    entries = list(reaction_totals.values())
    trace_floor_mol = 1e-15
    below_floor = max(gross_formation, gross_consumption, abs(observed_net)) < trace_floor_mol
    # A nearly cancelled radical turnover can have a tiny net despite large,
    # well-resolved opposing fluxes. Check against both scales and expose each
    # component rather than making an absolute net-only percentage claim.
    net_scale_tolerance_mol = 1e-3 * max(abs(observed_net), trace_floor_mol)
    # The BDF state tolerance is 1e-8. Its dense interpolant can differ
    # around fast radical transients even when five-point quadrature is used;
    # 1e-7 of turnover is an audit of pathway ranking, never net accuracy.
    gross_scale_tolerance_mol = 1e-7 * max(gross_formation, gross_consumption)
    closure_tolerance_mol = max(net_scale_tolerance_mol, gross_scale_tolerance_mol)
    closure_residual_mol = quadrature_net - observed_net
    status = (
        "below_numerical_trace_floor" if below_floor
        else "quadrature_reconciled"
        if abs(closure_residual_mol) <= net_scale_tolerance_mol
        else "gross_flux_reconciled_net_ill_conditioned"
        if abs(closure_residual_mol) <= gross_scale_tolerance_mol
        else "quadrature_unreconciled"
    )
    return {
        "species": species,
        "status": status,
        "numerical_trace_floor_mol": trace_floor_mol,
        "quadrature_closure_tolerance_mol": closure_tolerance_mol,
        "net_scale_closure_tolerance_mol": net_scale_tolerance_mol,
        "gross_scale_closure_tolerance_mol": gross_scale_tolerance_mol,
        "closure_error_relative_to_observed_net": (
            abs(closure_residual_mol) / abs(observed_net) if observed_net != 0 else None
        ),
        "closure_error_relative_to_larger_signed_flux": (
            abs(closure_residual_mol) / max(gross_formation, gross_consumption)
            if max(gross_formation, gross_consumption) > 0 else None
        ),
        "opposing_signed_reaction_contributions_over_abs_observed_net": (
            (gross_formation + gross_consumption) / abs(observed_net)
            if abs(observed_net) >= trace_floor_mol else None
        ),
        "opposing_contribution_ratio_is_yield_uncertainty_multiplier": False,
        "method": next(iter(method_counts)) if len(method_counts) == 1 else "mixed per-substep Gauss quadrature",
        "method_substep_counts": dict(method_counts),
        "substeps_traced": len(heavy_steps),
        f"positive_signed_net_reaction_{species.lower()}_mol": gross_formation,
        f"negative_signed_net_reaction_{species.lower()}_mol": gross_consumption,
        f"quadrature_net_gas_chemical_{species.lower()}_mol": quadrature_net,
        f"observed_net_gas_chemical_{species.lower()}_mol": observed_net,
        "total_quadrature_closure_residual_mol": closure_residual_mol,
        "largest_substep_quadrature_closure_residual_mol": max_step_residual,
        "signed_reaction_contributions": sorted(entries, key=lambda item: item["reaction_index"]),
        "top_formation_reactions": [] if below_floor else sorted(
            (item for item in entries if item["formation_mol"] > 0),
            key=lambda item: -item["formation_mol"],
        )[:10],
        "top_consumption_reactions": [] if below_floor else sorted(
            (item for item in entries if item["consumption_mol"] > 0),
            key=lambda item: -item["consumption_mol"],
        )[:10],
        "scope": (
            f"Signed net {mechanism_label} gas-reaction contributions for {species} only, not microscopic "
            "forward/reverse event totals. Electron fragments, "
            "gas-water exchange and aqueous acid-base/speciation are not "
            "individual gas reactions in this trace. This is a "
            "read-only mechanism diagnostic, not a validated yield."
        ),
    }


def _charged_hcnn_net_fate(hcnn_trace: dict[str, Any] | None) -> dict[str, Any] | None:
    """Rank model-internal *net* HCNN sinks, not collision branching odds."""
    if hcnn_trace is None:
        return None
    if hcnn_trace.get("status") not in (
        "quadrature_reconciled", "gross_flux_reconciled_net_ill_conditioned",
    ):
        return {
            "status": "unavailable_or_unreconciled_trace",
            "is_forward_branching_probability": False,
            "independent_kinetic_validation": False,
        }
    total = float(hcnn_trace["negative_signed_net_reaction_hcnn_mol"])
    if total <= 0:
        return {
            "status": "no_resolved_hcnn_consumption",
            "is_forward_branching_probability": False,
            "independent_kinetic_validation": False,
        }
    expected = {
        256: "HCNN + O <=> CO + H + N2",
        257: "HCNN + O <=> HCN + NO",
        260: "H + HCNN <=> CH2 + N2",
    }
    by_index = {
        int(item["reaction_index"]): item
        for item in hcnn_trace["signed_reaction_contributions"]
    }
    routes = []
    for index, equation in expected.items():
        item = by_index.get(index)
        if item is not None and item["equation"] != equation:
            raise ValueError(f"Charged HCNN fate index {index} no longer matches its equation")
        consumption = float(item["consumption_mol"]) if item else 0.0
        routes.append({
            "reaction_index": index,
            "equation": equation,
            "signed_net_hcnn_consumption_mol": consumption,
            "fraction_of_all_signed_net_hcnn_consumption": consumption / total,
        })
    return {
        "status": "model_internal_signed_net_fate_resolved",
        "all_signed_net_hcnn_consumption_mol": total,
        "selected_routes": routes,
        "other_routes_signed_net_hcnn_consumption_mol": max(
            0.0, total - sum(item["signed_net_hcnn_consumption_mol"] for item in routes),
        ),
        "is_forward_branching_probability": False,
        "independent_kinetic_validation": False,
        "scope": (
            "Fractions partition integrated signed net mechanism contributions. "
            "They are not microscopic forward branching ratios, measured rate "
            "constants, a product-yield fraction, or an uncertainty interval."
        ),
    }


def _aggregate_charged_electron_sources(
    snapshots: list[dict[str, Any]] | None,
) -> dict[str, Any] | None:
    if snapshots is None:
        return None
    if not snapshots or any(item.get("electron_sources") is None for item in snapshots):
        return {
            "status": "incomplete_charged_source_trace",
            "is_product_yield": False,
        }
    totals: dict[str, dict[str, Any]] = {}
    for snapshot in snapshots:
        for source in snapshot["electron_sources"]:
            channel_id = str(source["channel_id"])
            metadata = (
                source["target"], tuple(source["mapped_products"]),
                source["source_product_label"],
            )
            amount = float(source["accepted_extent_mol"])
            if not math.isfinite(amount) or amount < 0:
                raise ValueError("Charged electron source extent is invalid")
            entry = totals.setdefault(channel_id, {
                "channel_id": channel_id,
                "target": metadata[0],
                "mapped_products": list(metadata[1]),
                "source_product_label": metadata[2],
                "accepted_extent_mol": 0.0,
            })
            if metadata != (
                entry["target"], tuple(entry["mapped_products"]),
                entry["source_product_label"],
            ):
                raise ValueError("Charged electron channel identity changed across substeps")
            entry["accepted_extent_mol"] += amount
    return {
        "status": "integrated_accepted_electron_source_extents",
        "substeps_traced": len(snapshots),
        "channels": sorted(totals.values(), key=lambda item: -item["accepted_extent_mol"]),
        "is_product_yield": False,
        "scope": (
            "Accepted, energy- and inventory-limited external EEDF source extents "
            "in the charged branch. These do not include gas reaction fluxes or "
            "prove that a downstream product formed through a listed channel."
        ),
    }


def run_long_horizon(
    *, duration_min: float, cycles: int, substeps_per_cycle: int = 48,
    gas_water_kla_per_s: float = 1e-4,
    spark_coupling_factor: float = 1.0,
    trace_hcn_chemistry: bool = False,
    trace_h2cn_chemistry: bool = False,
    trace_charged_hcn_final_rates: bool = False,
    trace_charged_hcn_integrated: bool = False,
    charged_ch_water_rate_factor: float = 1.0,
    charged_ch_n2_rate_factor: float = 1.0,
    charged_hcn_h_association_rate_factor: float = 1.0,
    charged_hcnn_o_to_hcn_rate_factor: float = 1.0,
    charged_hcnn_h_to_ch2_n2_rate_factor: float = 1.0,
    gri_direct_hcn_rate_factor: float = 1.0,
    gri_hcn_h_association_rate_factor: float = 1.0,
    include_literature_h2cn_h_to_hcn: bool = False,
    include_hebrard_h2cn_n_branches: bool = False,
    include_hebrard_h2cn_h_isomer_branches: bool = False,
    include_hebrard_hnc_h_isomerization: bool = False,
    screen_glycolonitrile_endpoint: bool = False,
    formaldehyde_hydration_model: str = "winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    remove_carbon_for_control: bool = False,
    remove_nitrogen_for_control: bool = False,
) -> dict[str, Any]:
    feed = json.loads(FEED_CARD.read_text(encoding="utf-8"))
    reported_duration_min = float(feed["reported_spark_duration_min"])
    if not math.isfinite(duration_min) or not 0 < duration_min <= reported_duration_min:
        raise ValueError("Horizon must be positive and no longer than the published spark duration")
    if not 1 <= cycles <= 100_000 or not 1 <= substeps_per_cycle <= 48:
        raise ValueError("Cycles or substeps are outside the apparatus schema range")
    seconds_per_cycle = duration_min * 60.0 / cycles
    if seconds_per_cycle > 3600:
        raise ValueError("Use enough cycles to keep each cycle at or below one hour")
    if not math.isfinite(spark_coupling_factor) or spark_coupling_factor <= 0:
        raise ValueError("Spark coupling sensitivity factor must be finite and positive")
    for factor in (
        gri_direct_hcn_rate_factor, gri_hcn_h_association_rate_factor,
        charged_ch_water_rate_factor, charged_ch_n2_rate_factor,
        charged_hcn_h_association_rate_factor, charged_hcnn_o_to_hcn_rate_factor,
        charged_hcnn_h_to_ch2_n2_rate_factor,
    ):
        if not math.isfinite(factor) or factor < 0:
            raise ValueError("Diagnostic rate factors must be finite and nonnegative")
    if include_literature_h2cn_h_to_hcn and include_hebrard_h2cn_h_isomer_branches:
        raise ValueError("Choose either Pearce or Hebrard H2CN + H diagnostic branches")
    if include_hebrard_hnc_h_isomerization and not include_hebrard_h2cn_h_isomer_branches:
        raise ValueError("HNC + H isomerization requires the opt-in HNC species branch")
    rate_multipliers = tuple(
        (index, factor)
        for index, factor in (
            (275, gri_direct_hcn_rate_factor),
            (236, gri_hcn_h_association_rate_factor),
        )
        if factor != 1.0
    )
    charged_rate_multipliers = tuple(
        (index, factor)
        for index, factor in (
            (126, charged_ch_water_rate_factor),
            (240, charged_ch_n2_rate_factor),
            (236, charged_hcn_h_association_rate_factor),
            (257, charged_hcnn_o_to_hcn_rate_factor),
            (260, charged_hcnn_h_to_ch2_n2_rate_factor),
        )
        if factor != 1.0
    )
    payload = make_probe_payload(
        duration_s=min(duration_min * 60.0, 3600.0),
        cycles=1, substeps=substeps_per_cycle,
        gas_water_kla_per_s=gas_water_kla_per_s,
        gas_water_solver="coupled_bdf",
        primary_source_forcing="inventory_limited",
        formaldehyde_hydration_model=formaldehyde_hydration_model,
        remove_carbon_for_control=remove_carbon_for_control,
        remove_nitrogen_for_control=remove_nitrogen_for_control,
    )
    payload["id"] = (
        "barth_feed_long_horizon_carbon_removed_control"
        if remove_carbon_for_control else
        "barth_feed_long_horizon_nitrogen_removed_control"
        if remove_nitrogen_for_control else
        "barth_feed_only_long_horizon_stress_not_yield_validation"
    )
    payload["apparatus"]["cycles"] = cycles
    payload["apparatus"]["discharge_duration"] = {
        "value": seconds_per_cycle, "unit": "s",
    }
    coupling = float(payload["spark"]["coupling_efficiency"]) * spark_coupling_factor
    if coupling > 1:
        raise ValueError("Spark coupling sensitivity exceeds physical fraction one")
    payload["spark"]["coupling_efficiency"] = coupling
    payload["spark"]["source"] = (
        "Nominal generator settings with hypothetical coupling sensitivity; "
        "not a measured Barth wet-flask discharge"
    )
    payload["spark"]["basis"] = "hypothetical"
    result = simulate_prebiotic_apparatus(
        PrebioticApparatusRequest.model_validate(payload),
        diagnostic_pathway_species_each_substep=(
            (("HCN",) if trace_hcn_chemistry else ())
            + (("H2CN",) if trace_h2cn_chemistry else ())
            + (("HNC",) if include_hebrard_h2cn_h_isomer_branches else ())
        ),
        diagnostic_gri_reaction_multipliers=rate_multipliers,
        diagnostic_literature_h2cn_h_to_hcn=include_literature_h2cn_h_to_hcn,
        diagnostic_hebrard_h2cn_n_branches=include_hebrard_h2cn_n_branches,
        diagnostic_hebrard_h2cn_h_isomer_branches=include_hebrard_h2cn_h_isomer_branches,
        diagnostic_hebrard_hnc_h_isomerization=include_hebrard_hnc_h_isomerization,
        diagnostic_charged_hcn_final_rate_each_substep=(
            trace_charged_hcn_final_rates
        ),
        diagnostic_charged_hcn_integrated_each_substep=(
            trace_charged_hcn_integrated
        ),
        diagnostic_pavan_reaction_multipliers=charged_rate_multipliers,
    )
    glycolonitrile_endpoint = None
    if screen_glycolonitrile_endpoint:
        from app.prebiotic.glycolonitrile_screen import evaluate_glycolonitrile_endpoint_screen
        glycolonitrile_endpoint = evaluate_glycolonitrile_endpoint_screen(
            result, duration_s=duration_min * 60.0,
        )
    source_energy_audit = audit_source_thermochemistry(
        result, float(payload["apparatus"]["discharge_temperature"]["value"]),
    )
    charged = result["plasma_model"].get("propagated_charged_branch") or {}
    heavy_steps = result["plasma_model"]["heavy_species_cycle_summaries"]
    charged_hcn_final_rate_snapshots = result["plasma_model"].get(
        "diagnostic_charged_hcn_final_rate_snapshots"
    )
    charged_hcn_integrated_snapshots = result["plasma_model"].get(
        "diagnostic_charged_hcn_integrated_snapshots"
    )
    charged_pathway_traces: dict[str, dict[str, Any]] = {}
    charged_electron_sources = None
    if trace_charged_hcn_integrated:
        charged_electron_sources = _aggregate_charged_electron_sources(
            charged_hcn_integrated_snapshots
        )
        if charged_hcn_integrated_snapshots and all(
            item.get("diagnostic") is not None
            for item in charged_hcn_integrated_snapshots
        ):
            charged_steps = [
                {"integrated_pathway_diagnostics": item["diagnostic"]}
                for item in charged_hcn_integrated_snapshots
            ]
            charged_pathway_traces = {
                species: _aggregate_species_pathways(charged_steps, species, "Pavan charged")
                for species in (
                    "HCN", "HCNN", "CH", "CH2", "CH2(S)", "CH2OH", "CH3O",
                    "CH2O", "HCO", "CO",
                )
            }
        else:
            charged_pathway_traces["HCN"] = {
                "status": "incomplete_charged_branch_diagnostic",
                "reason": "One or more charged substeps lacked an integrated trace",
            }
    hcn_trace = _aggregate_species_pathways(heavy_steps, "HCN") if trace_hcn_chemistry else None
    h2cn_trace = _aggregate_species_pathways(heavy_steps, "H2CN") if trace_h2cn_chemistry else None
    hnc_trace = (
        _aggregate_species_pathways(heavy_steps, "HNC")
        if include_hebrard_h2cn_h_isomer_branches else None
    )
    source_extents = [
        item for step in heavy_steps for item in step["source_extents"]
    ]
    co2_o1d_excluded = [
        (item, step["collision_energy_allocation"]["global_collision_energy_scale"])
        for step in heavy_steps
        for item in step["collision_energy_allocation"]["excluded_fragment_channels"]
        if item["product_label"] == "CO(X) + O(1D)"
        and item["reason"] == "unresolved_excited_state"
    ]
    o1d_fates = [
        item["o1d_fate_diagnostic"]
        for item in result["plasma_model"]["cycle_summaries"]
        if item.get("o1d_fate_diagnostic") is not None
    ]
    n2d_fates = [
        item["n2d_fate_diagnostic"]
        for item in result["plasma_model"]["cycle_summaries"]
        if item.get("n2d_fate_diagnostic") is not None
    ]
    n2d_total_outcomes = [
        source
        for item in result["plasma_model"]["cycle_summaries"]
        for source in item.get("source_terms", [])
        if source.get("target") == "N2"
        and source.get("product_state_unresolved")
        and source.get("product_state_energy_support") is not None
    ]
    n2d_total_opportunities_mol = sum(
        float(item["uncoupled_active_window_opportunity_mol"])
        for item in n2d_total_outcomes
    )
    n2d_energy_support = {
        window: (
            sum(
                float(item["uncoupled_active_window_opportunity_mol"])
                * float(item["product_state_energy_support"][
                    f"fraction_of_computed_rate_{window}_measured_window"
                ])
                for item in n2d_total_outcomes
                if item["product_state_energy_support"][
                    f"fraction_of_computed_rate_{window}_measured_window"
                ] is not None
            ) / n2d_total_opportunities_mol
            if n2d_total_opportunities_mol > 0 else None
        )
        for window in ("below", "within", "above")
    }
    n2d_ground_atom_energetic_fraction = (
        sum(
            float(item["uncoupled_active_window_opportunity_mol"])
            * float(item["ground_atom_energetic_window"]["rate_fraction_below_that_threshold"])
            for item in n2d_total_outcomes
        ) / n2d_total_opportunities_mol
        if n2d_total_opportunities_mol > 0 else None
    )
    n2d_analytic_weighted = [
        (
            float(item["uncoupled_active_window_opportunity_mol"])
            * float(item["analytic_fit_diagnostic"]["rate_coefficient_m3_s"])
            / float(item["rate_coefficient_m3_s"]),
            item["analytic_fit_diagnostic"]["product_state_energy_support"],
        )
        for item in n2d_total_outcomes
        if item.get("analytic_fit_diagnostic") is not None
        and float(item["rate_coefficient_m3_s"]) > 0
        and float(item["analytic_fit_diagnostic"]["rate_coefficient_m3_s"]) > 0
    ]
    n2d_analytic_opportunities_mol = sum(
        weight for weight, _ in n2d_analytic_weighted
    )
    n2d_analytic_energy_support = {
        window: (
            sum(
                weight * float(support[f"fraction_of_computed_rate_{window}_measured_window"])
                for weight, support in n2d_analytic_weighted
            ) / n2d_analytic_opportunities_mol
            if n2d_analytic_opportunities_mol > 0 else None
        )
        for window in ("below", "within", "above")
    }
    n2d_analytic_ground_atom_energetic_fraction = (
        sum(
            float(item["uncoupled_active_window_opportunity_mol"])
            * float(item["analytic_fit_diagnostic"]["rate_coefficient_m3_s"])
            / float(item["rate_coefficient_m3_s"])
            * float(item["analytic_fit_diagnostic"]["ground_atom_energetic_window"][
                "rate_fraction_below_that_threshold"
            ])
            for item in n2d_total_outcomes
            if item.get("analytic_fit_diagnostic") is not None
            and float(item["rate_coefficient_m3_s"]) > 0
        ) / n2d_analytic_opportunities_mol
        if n2d_analytic_opportunities_mol > 0 else None
    )
    n2d_water_fractions = [
        float(channel["conditional_first_collision_fraction"])
        for fate in n2d_fates
        if fate["status"] == "conditional_first_collision_fate_resolved"
        for channel in fate["channels"]
        if channel["bath_species"] == "H2O"
    ]
    charged_compartments = charged.get("final_compartments_mol") or {}
    charged_aqueous_cyanide = sum(
        float(charged_compartments.get(place, {}).get("HCN", 0.0))
        + float(charged_compartments.get(place, {}).get("CN-", 0.0))
        for place in ("boiler_liquid", "condensate_trap", "collected_condensate")
    )
    primary_aqueous_cyanide = float(result["hcn_formation"]["final_hcn_plus_cyanide_aqueous_mol"])
    reporting_floor = float(result["hcn_formation"]["numerical_reporting_floor_mol"])
    controls_pass = (
        result["material_ledger"]["closed"] is True
        and (charged.get("material_ledger") or {}).get("closed") is True
        and charged.get("failed_substep_count") == 0
        and source_energy_audit["status"] == "accepted_source_thermochemical_floor_passed"
        and (hcn_trace is None or hcn_trace["status"] != "quadrature_unreconciled")
        and (h2cn_trace is None or h2cn_trace["status"] != "quadrature_unreconciled")
        and (hnc_trace is None or hnc_trace["status"] != "quadrature_unreconciled")
    )
    if remove_carbon_for_control or remove_nitrogen_for_control:
        controls_pass = (
            controls_pass
            and primary_aqueous_cyanide <= reporting_floor
            and charged_aqueous_cyanide <= reporting_floor
        )
    return {
        "status": (
            "long_horizon_numerical_screen_passed_not_yield_validation"
            if controls_pass else "long_horizon_numerical_screen_failed"
        ),
        "carbon_removed_negative_control": remove_carbon_for_control,
        "nitrogen_removed_negative_control": remove_nitrogen_for_control,
        "duration_min": duration_min,
        "reported_spark_duration_min": reported_duration_min,
        "cycles": cycles,
        "substeps_per_cycle": substeps_per_cycle,
        "total_integration_steps": cycles * substeps_per_cycle,
        "formaldehyde_hydration_model": formaldehyde_hydration_model,
        "minimum_modeled_ideal_pH": result["operating_conditions"]["formaldehyde_hydration_minimum_ideal_pH"],
        "pH_extrapolation_used": result["scientific_firewall"]["formaldehyde_hydration_pH_extrapolation_used"],
        "glycolonitrile_endpoint_screen": glycolonitrile_endpoint,
        "primary_aqueous_HCN_plus_CN_minus_mol": primary_aqueous_cyanide,
        "charged_aqueous_HCN_plus_CN_minus_mol": charged_aqueous_cyanide,
        "numerical_reporting_floor_mol": reporting_floor,
        "primary_material_ledger_closed": result["material_ledger"]["closed"],
        "primary_gas_ideal_energy_species_covered_each_step": all(
            bool(step["ideal_gas_internal_energy_audit"]["all_nonzero_gas_species_covered"])
            for step in heavy_steps
        ),
        "complete_thermal_energy_balance_closed": result["energy_ledger"]["complete_thermal_energy_balance_closed"],
        "charged_material_ledger_closed": (charged.get("material_ledger") or {}).get("closed"),
        "charged_failed_substeps": charged.get("failed_substep_count"),
        "hcn_gas_reaction_trace": hcn_trace,
        "h2cn_gas_reaction_trace": h2cn_trace,
        "hnc_gas_reaction_trace": hnc_trace,
        "gri_diagnostic_reaction_rate_factors": {
            "CH3_plus_N_to_HCN_plus_H2": gri_direct_hcn_rate_factor,
            "H_plus_HCN_to_H2CN": gri_hcn_h_association_rate_factor,
        },
        "gri_rate_stress_is_a_fitted_or_validated_mechanism": False,
        "unmodified_gri_h_hcn_initial_feed_reference_audit": (
            audit_hcn_association_initial_feed()
            if not (remove_carbon_for_control or remove_nitrogen_for_control) else None
        ),
        "unmodified_gri_ch_n2_initial_feed_reference_audit": (
            audit_ch_n2_association()
            if not (remove_carbon_for_control or remove_nitrogen_for_control) else None
        ),
        "unmodified_gri_ch_water_initial_feed_reference_audit": (
            audit_ch_water_competition()
            if not (remove_carbon_for_control or remove_nitrogen_for_control) else None
        ),
        "ch_water_reference_audit_matches_active_charged_rate_profile": (
            charged_ch_water_rate_factor == 1.0
            and not (remove_carbon_for_control or remove_nitrogen_for_control)
        ),
        "ch_n2_reference_audit_matches_active_charged_rate_profile": (
            charged_ch_n2_rate_factor == 1.0
            and not (remove_carbon_for_control or remove_nitrogen_for_control)
        ),
        "h_hcn_reference_audit_matches_active_gri_rate_profile": (
            gri_hcn_h_association_rate_factor == 1.0
            and not (remove_carbon_for_control or remove_nitrogen_for_control)
        ),
        "literature_h2cn_h_to_hcn_in_primary_diagnostic": include_literature_h2cn_h_to_hcn,
        "literature_h2cn_h_to_hcn_rate_source": (
            "Pearce, Ayers & Pudritz 2019, DOI 10.1021/acs.jpca.8b11323, "
            "1.8e-11 cm3 molecule^-1 s^-1 over 50-400 K"
            if include_literature_h2cn_h_to_hcn else None
        ),
        "hebrard_h2cn_n_branches_in_primary_diagnostic": include_hebrard_h2cn_n_branches,
        "hebrard_h2cn_n_branch_source": (
            "Hebrard et al. 2012, DOI 10.1051/0004-6361/201218837, "
            "room-temperature estimated N(4S)+H2CN branches: "
            "3e-11 to N2+CH2 and 1e-11 to NH+HCN cm3 molecule^-1 s^-1"
            if include_hebrard_h2cn_n_branches else None
        ),
        "hebrard_h2cn_h_isomer_branches_in_primary_diagnostic": (
            include_hebrard_h2cn_h_isomer_branches
        ),
        "hebrard_h2cn_h_isomer_branch_source": (
            "Hebrard et al. 2012, DOI 10.1051/0004-6361/201218837, "
            "estimated 298 K H2CN+H exits: 6e-11 to HCN+H2 and "
            "1.2e-11 to HNC+H2 cm3 molecule^-1 s^-1; "
            "NIST-grounded 290-310 K HNC thermochemistry"
            if include_hebrard_h2cn_h_isomer_branches else None
        ),
        "hnc_gas_water_fate_reviewed": (
            False if include_hebrard_h2cn_h_isomer_branches else None
        ),
        "hebrard_hnc_h_isomerization_in_primary_diagnostic": (
            include_hebrard_hnc_h_isomerization
        ),
        "charged_hcn_endpoint_rate_trace": charged_hcn_final_rate_snapshots,
        "charged_hcn_integrated_pathway_trace": charged_hcn_integrated_snapshots,
        "charged_hcn_integrated_pathway_summary": charged_pathway_traces.get("HCN"),
        "charged_hcnn_integrated_pathway_summary": charged_pathway_traces.get("HCNN"),
        "charged_ch_integrated_pathway_summary": charged_pathway_traces.get("CH"),
        "charged_ch2_integrated_pathway_summary": charged_pathway_traces.get("CH2"),
        "charged_ch2_s_integrated_pathway_summary": charged_pathway_traces.get("CH2(S)"),
        "charged_ch2oh_integrated_pathway_summary": charged_pathway_traces.get("CH2OH"),
        "charged_ch3o_integrated_pathway_summary": charged_pathway_traces.get("CH3O"),
        "charged_ch2o_integrated_pathway_summary": charged_pathway_traces.get("CH2O"),
        "charged_hco_integrated_pathway_summary": charged_pathway_traces.get("HCO"),
        "charged_co_integrated_pathway_summary": charged_pathway_traces.get("CO"),
        "charged_electron_source_summary": charged_electron_sources,
        "charged_hcnn_net_fate_diagnostic": _charged_hcnn_net_fate(
            charged_pathway_traces.get("HCNN")
        ),
        "charged_diagnostic_reaction_rate_factors": {
            "CH_plus_H2O_to_CH2O_plus_H": charged_ch_water_rate_factor,
            "CH_plus_N2_to_HCNN": charged_ch_n2_rate_factor,
            "H_plus_HCN_to_H2CN": charged_hcn_h_association_rate_factor,
            "HCNN_plus_O_to_HCN_plus_NO": charged_hcnn_o_to_hcn_rate_factor,
            "H_plus_HCNN_to_CH2_plus_N2": charged_hcnn_h_to_ch2_n2_rate_factor,
        },
        "charged_rate_stress_is_a_fitted_or_validated_mechanism": False,
        "charged_mechanism_overlap": charged.get("mechanism_overlap"),
        "hebrard_hnc_h_isomerization_source": (
            "Hebrard et al. 2012, DOI 10.1051/0004-6361/201218837, "
            "estimated relaxed-HNC H+HNC -> H+HCN rate "
            "4.03e-11 exp(-1200/T) cm3 molecule^-1 s^-1; "
            "stated uncertainty factor 10 at 300 K"
            if include_hebrard_hnc_h_isomerization else None
        ),
        "charged_shadow_uses_gri_rate_stress": False,
        "spark_coupling_factor_from_hypothetical_baseline": spark_coupling_factor,
        "declared_electrical_coupling_efficiency": coupling,
        "source_diagnostic": {
            "thermochemical_floor_status": source_energy_audit["status"],
            "ordinary_fragment_enthalpy_deficit_count": len(source_energy_audit["deficits"]),
            "n2_parent_event_capacity_deficit_count": len(source_energy_audit["n2_parent_coverage_deficits"]),
            "n2_parent_minimum_capacity_to_outcome_ratio": source_energy_audit[
                "n2_parent_minimum_capacity_to_outcome_ratio"
            ],
            "accepted_neutral_extent_mol_by_target": {
                target: sum(
                    float(item.get("actual_integrated_extent_mol", item["accepted_extent_mol"]))
                    for item in source_extents if item["target"] == target
                )
                for target in ("N2", "CO2", "O2", "H2", "H2O")
            },
            "excluded_CO2_to_CO_X_plus_O1D_channel_count": len(co2_o1d_excluded),
            "excluded_CO2_to_CO_X_plus_O1D_raw_opportunities_mol": sum(
                float(item["raw_opportunity_mol"]) for item, _ in co2_o1d_excluded
            ),
            "excluded_CO2_to_CO_X_plus_O1D_energy_scaled_opportunities_mol": sum(
                float(item["raw_opportunity_mol"]) * float(scale)
                for item, scale in co2_o1d_excluded
            ),
            "excluded_opportunities_are_product_moles": False,
            "n2d_fate_status_counts": dict(Counter(
                item["status"] for item in n2d_fates
            )),
            "n2d_total_outcome_uncoupled_opportunities_mol": n2d_total_opportunities_mol,
            "n2d_uncoupled_rate_fraction_by_cosby_product_state_energy_window": n2d_energy_support,
            "n2d_ground_N4S_energetic_necessity_window_eV": [0.0, 14.526],
            "n2d_uncoupled_rate_fraction_in_conditional_ground_N4S_window": (
                n2d_ground_atom_energetic_fraction
            ),
            "n2d_analytic_fit_uncoupled_opportunities_mol": n2d_analytic_opportunities_mol,
            "n2d_analytic_fit_uncoupled_rate_fraction_by_cosby_energy_window": (
                n2d_analytic_energy_support
            ),
            "n2d_analytic_fit_uncoupled_rate_fraction_in_conditional_ground_N4S_window": (
                n2d_analytic_ground_atom_energetic_fraction
            ),
            "n2d_ground_N4S_energetic_window_used_as_fragment_source": False,
            "n2d_analytic_fit_used_to_propagate_chemistry": False,
            "n2d_cosby_product_state_measurement_energy_window_eV": [18.5, 148.5],
            "n2d_electron_energy_window_fractions_are_state_branching": False,
            "n2d_water_conditional_first_collision_fraction_range": (
                [min(n2d_water_fractions), max(n2d_water_fractions)]
                if n2d_water_fractions else None
            ),
            "n2d_unreviewed_bath_species": sorted({
                species for item in n2d_fates
                for species in item.get("unreviewed_bath_species", [])
            }),
            "maximum_n2d_unreviewed_bath_mole_fraction": max(
                (float(item.get("unreviewed_bath_mole_fraction") or 0.0)
                 for item in n2d_fates), default=0.0,
            ),
            "maximum_n2d_unreviewed_bath_amount_mol_by_species": {
                species: max(
                    float(item.get("unreviewed_bath_amounts_mol", {}).get(species, 0.0))
                    for item in n2d_fates
                )
                for species in sorted({
                    species for item in n2d_fates
                    for species in item.get("unreviewed_bath_amounts_mol", {})
                })
            },
            "minimum_n2d_unreviewed_rate_for_equal_reviewed_hazard_cm3_molecule_s": {
                species: min(
                    float(item["unreviewed_rate_for_equal_reviewed_hazard_cm3_molecule_s"][species])
                    for item in n2d_fates
                    if species in item.get("unreviewed_rate_for_equal_reviewed_hazard_cm3_molecule_s", {})
                )
                for species in sorted({
                    species for item in n2d_fates
                    for species in item.get("unreviewed_rate_for_equal_reviewed_hazard_cm3_molecule_s", {})
                })
            },
            "n2d_equal_hazard_rates_are_measured_or_bounded": False,
            "n2d_source_state_and_products_propagated": False,
            "o1d_fate_status_counts": dict(Counter(
                item["status"] for item in o1d_fates
            )),
            "o1d_unreviewed_bath_species": sorted({
                species for item in o1d_fates
                for species in item.get("unreviewed_bath_species", [])
            }),
            "maximum_o1d_unreviewed_bath_mole_fraction": max(
                (float(item.get("unreviewed_bath_mole_fraction") or 0.0)
                 for item in o1d_fates), default=0.0,
            ),
            "maximum_o1d_hypothetical_unknown_hazard_fraction_by_rate": {
                f"{rate:g} cm3/molecule/s": max(
                    (
                        float(stress["hypothetical_unknown_hazard_fraction"])
                        for item in o1d_fates
                        for stress in item.get("unreviewed_rate_stress", [])
                        if stress["hypothetical_rate_cm3_molecule_s"] == rate
                        and stress["hypothetical_unknown_hazard_fraction"] is not None
                    ), default=0.0,
                )
                for rate in (1e-9, 1e-8)
            },
            "maximum_o1d_unreviewed_bath_amount_mol_by_species": {
                species: max(
                    float(item.get("unreviewed_bath_amounts_mol", {}).get(species, 0.0))
                    for item in o1d_fates
                )
                for species in sorted({
                    species for item in o1d_fates
                    for species in item.get("unreviewed_bath_amounts_mol", {})
                })
            },
        },
        "deposited_energy_basis": result["energy_ledger"]["deposited_energy_basis"],
        "conditional_water_thermal_capacity_screen": audit_barth_water_thermal_capacity(
            duration_min=duration_min,
        ),
        "declared_electrical_estimate_J": result["energy_ledger"]["deposited_electrical_energy_J"],
        "budgeted_deposited_energy_J": result["energy_ledger"]["deposited_energy_used_J"],
        "source_threshold_energy_lower_bound_J": result["energy_ledger"]["plasma_source_threshold_energy_lower_bound_J"],
        "electrical_source_budget_status": result["energy_ledger"]["status"],
        "experimental_hcn_assay_loaded": False,
        "science_boundary": (
            "Feed-only numerical stress with hypothetical active spark volume, pulse width, "
            "electrical coupling, gas-water kLa, and explicitly reported hydration pH domain. The held-out "
            "cyanide assay is not loaded or used for fitting; this is not a yield validation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--duration-min", type=float, default=120.0)
    parser.add_argument("--cycles", type=int, default=2)
    parser.add_argument("--substeps-per-cycle", type=int, default=48)
    parser.add_argument("--gas-water-kla-per-s", type=float, default=1e-4)
    parser.add_argument("--spark-coupling-factor", type=float, default=1.0)
    parser.add_argument("--trace-hcn-chemistry", action="store_true")
    parser.add_argument("--trace-h2cn-chemistry", action="store_true")
    parser.add_argument("--trace-charged-hcn-final-rates", action="store_true")
    parser.add_argument("--trace-charged-hcn-integrated", action="store_true")
    parser.add_argument("--charged-ch-water-rate-factor", type=float, default=1.0)
    parser.add_argument("--charged-ch-n2-rate-factor", type=float, default=1.0)
    parser.add_argument("--charged-hcn-h-association-rate-factor", type=float, default=1.0)
    parser.add_argument("--charged-hcnn-o-to-hcn-rate-factor", type=float, default=1.0)
    parser.add_argument("--charged-hcnn-h-to-ch2-n2-rate-factor", type=float, default=1.0)
    parser.add_argument("--gri-direct-hcn-rate-factor", type=float, default=1.0)
    parser.add_argument("--gri-hcn-h-association-rate-factor", type=float, default=1.0)
    parser.add_argument("--include-literature-h2cn-h-to-hcn", action="store_true")
    parser.add_argument("--include-hebrard-h2cn-n-branches", action="store_true")
    parser.add_argument("--include-hebrard-h2cn-h-isomer-branches", action="store_true")
    parser.add_argument("--include-hebrard-hnc-h-isomerization", action="store_true")
    parser.add_argument("--screen-glycolonitrile-endpoint", action="store_true")
    parser.add_argument("--remove-carbon-for-control", action="store_true")
    parser.add_argument("--remove-nitrogen-for-control", action="store_true")
    parser.add_argument("--formaldehyde-hydration-model", choices=(
        "winkelman_2002_source_bounded",
        "winkelman_2002_pH4_to_5_extrapolation_diagnostic",
    ), default="winkelman_2002_pH4_to_5_extrapolation_diagnostic")
    args = parser.parse_args()
    report = run_long_horizon(
        duration_min=args.duration_min,
        cycles=args.cycles,
        substeps_per_cycle=args.substeps_per_cycle,
        gas_water_kla_per_s=args.gas_water_kla_per_s,
        spark_coupling_factor=args.spark_coupling_factor,
        trace_hcn_chemistry=args.trace_hcn_chemistry,
        trace_h2cn_chemistry=args.trace_h2cn_chemistry,
        trace_charged_hcn_final_rates=args.trace_charged_hcn_final_rates,
        trace_charged_hcn_integrated=args.trace_charged_hcn_integrated,
        charged_ch_water_rate_factor=args.charged_ch_water_rate_factor,
        charged_ch_n2_rate_factor=args.charged_ch_n2_rate_factor,
        charged_hcn_h_association_rate_factor=args.charged_hcn_h_association_rate_factor,
        charged_hcnn_o_to_hcn_rate_factor=args.charged_hcnn_o_to_hcn_rate_factor,
        charged_hcnn_h_to_ch2_n2_rate_factor=args.charged_hcnn_h_to_ch2_n2_rate_factor,
        gri_direct_hcn_rate_factor=args.gri_direct_hcn_rate_factor,
        gri_hcn_h_association_rate_factor=args.gri_hcn_h_association_rate_factor,
        include_literature_h2cn_h_to_hcn=args.include_literature_h2cn_h_to_hcn,
        include_hebrard_h2cn_n_branches=args.include_hebrard_h2cn_n_branches,
        include_hebrard_h2cn_h_isomer_branches=(
            args.include_hebrard_h2cn_h_isomer_branches
        ),
        include_hebrard_hnc_h_isomerization=args.include_hebrard_hnc_h_isomerization,
        screen_glycolonitrile_endpoint=args.screen_glycolonitrile_endpoint,
        formaldehyde_hydration_model=args.formaldehyde_hydration_model,
        remove_carbon_for_control=args.remove_carbon_for_control,
        remove_nitrogen_for_control=args.remove_nitrogen_for_control,
    )
    print(json.dumps(report, indent=2))
    if report["status"] == "long_horizon_numerical_screen_failed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
