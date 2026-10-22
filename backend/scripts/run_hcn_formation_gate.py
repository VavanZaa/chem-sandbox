"""Unseeded input-control sweep for simulated HCN formation.

The apparatus receives no initial HCN amount or measured HCN-yield target.
HCN is inspected after each run, but the published GRI/Pavan gas mechanisms
already contain HCN species and reactions. Thus the feed/outcome is held out,
not the product identity itself. This is a formation gate, not a de novo
reaction discovery or calibrated experimental yield test.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path
from typing import Any

from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.boltzmann import IST_H2_COMMUNITY_CH4_DATASET_ID
from app.prebiotic.schema import PrebioticApparatusRequest


PROFILE = (
    Path(__file__).resolve().parents[1]
    / "data" / "benchmarks" / "stribling_miller_measurement_profile.json"
)
CASE_NAMES = ("spark_on", "spark_off", "no_methane", "no_nitrogen", "half_methane")


def _gas_amount(payload: dict[str, Any], species: str) -> dict[str, Any]:
    return next(
        item for item in payload["initial_inventory"]
        if item["compartment"] == "gas_loop" and item["species_id"] == species
    )


def _replace_gas_with_hydrogen(payload: dict[str, Any], species: str, fraction: float) -> None:
    source = _gas_amount(payload, species)["amount"]
    hydrogen = _gas_amount(payload, "H2")["amount"]
    removed = float(source["value"]) * fraction
    source["value"] = float(source["value"]) - removed
    hydrogen["value"] = float(hydrogen["value"]) + removed


def make_formation_scenario(*, cycles: int = 1, substeps: int = 2) -> dict[str, Any]:
    """Build the fixed, hypothetical formation scenario without measured HCN data."""

    if not 1 <= cycles <= 168 or not 1 <= substeps <= 48:
        raise ValueError("Cycles must be 1-168 and substeps 1-48")
    base = json.loads(PROFILE.read_text(encoding="utf-8"))
    base["id"] = "hcn_formation_control_scenario"
    base["apparatus"].update({
        "cycles": cycles,
        "operation_mode": "continuous_overlap" if substeps > 1 else "staged",
        "integration_substeps_per_cycle": substeps,
        "primary_source_forcing": "time_averaged",
        "aqueous_acid_base_model": "hcn_nh3_water_equilibrium",
        "condensate_nonwater_return_fraction": 0.0,
        "condenser_phase_transfer_model": "coupled_flow_transport",
        "condenser_contact_water_holdup": {
            "value": 0.1, "unit": "L",
            "source": "Illustrative condenser water contact hold-up for model sensitivity",
            "method": "declared scenario, not a measured reactor property",
            "basis": "hypothetical",
        },
    })
    base["spark"]["boltzmann"]["collision_set"] = IST_H2_COMMUNITY_CH4_DATASET_ID
    return base


def run_gate(
    *, cycles: int = 1, substeps: int = 2,
    case_names: tuple[str, ...] = CASE_NAMES,
) -> dict[str, Any]:
    if not case_names or len(set(case_names)) != len(case_names):
        raise ValueError("Selected cases must be nonempty and unique")
    if set(case_names) - set(CASE_NAMES):
        raise ValueError("Selected cases include an unknown formation control")
    base = make_formation_scenario(cycles=cycles, substeps=substeps)

    cases: dict[str, dict[str, Any]] = {}
    for name in case_names:
        payload = deepcopy(base)
        payload["id"] = f"hcn_control_{name}"
        if name == "spark_off":
            payload["spark"] = {"mode": "off"}
        elif name == "no_methane":
            _replace_gas_with_hydrogen(payload, "CH4", 1.0)
        elif name == "no_nitrogen":
            _replace_gas_with_hydrogen(payload, "N2", 1.0)
        elif name == "half_methane":
            _replace_gas_with_hydrogen(payload, "CH4", 0.5)
        try:
            result = simulate_prebiotic_apparatus(
                PrebioticApparatusRequest.model_validate(payload)
            )
        except Exception as error:
            raise RuntimeError(f"HCN formation control {name} failed") from error
        formation = result["hcn_formation"]
        final_compartments = result["final_compartments_mol"]
        liquid_places = ("boiler_liquid", "condensate_trap", "collected_condensate")
        liquid_carbonate = {
            species: sum(
                float(final_compartments.get(place, {}).get(species, 0.0))
                for place in liquid_places
            )
            for species in ("CO2", "HCO3-", "CO3--")
        }
        liquid_ammonia = {
            species: sum(
                float(final_compartments.get(place, {}).get(species, 0.0))
                for place in liquid_places
            )
            for species in ("NH3", "NH4+")
        }
        liquid_pH = {
            place: result["aqueous_cyanide_speciation"]["compartments"][place].get(
                "pH_ideal_aqueous_model"
            )
            for place in liquid_places
        }
        charged_branch = result["plasma_model"].get("propagated_charged_branch") or {}
        charged_ledger = charged_branch.get("material_ledger")
        charged_compartments = charged_branch.get("final_compartments_mol") or {}
        charged_cyanide_by_compartment = {
            place: float(species.get("HCN", 0.0)) + float(species.get("CN-", 0.0))
            for place, species in charged_compartments.items()
        }
        charged_cyanide_total = (
            sum(charged_cyanide_by_compartment.values())
            if charged_compartments else None
        )
        charged_cyanide_aqueous = (
            sum(
                charged_cyanide_by_compartment.get(place, 0.0)
                for place in ("boiler_liquid", "condensate_trap", "collected_condensate")
            )
            if charged_compartments else None
        )
        cases[name] = {
            "status": formation["status"],
            "net_hcn_generated_mol": formation["net_hcn_generated_mol"],
            "net_hcn_generated_nmol_per_declared_J": formation[
                "net_hcn_generated_nmol_per_declared_J"
            ],
            "final_hcn_aqueous_mol": formation["final_hcn_aqueous_mol"],
            "final_hcn_plus_cyanide_aqueous_mol": formation["final_hcn_plus_cyanide_aqueous_mol"],
            "final_hcn_plus_cyanide_aqueous_nmol_per_declared_J": formation[
                "final_hcn_plus_cyanide_aqueous_nmol_per_declared_J"
            ],
            "final_cyanide_ion_total_mol": formation["final_cyanide_ion_total_mol"],
            "final_liquid_carbonate_family_mol": liquid_carbonate,
            "final_liquid_ammonia_family_mol": liquid_ammonia,
            "final_liquid_pH_ideal_aqueous_model": liquid_pH,
            "final_hcn_total_mol": formation["final_hcn_total_mol"],
            "net_hcn_per_feed_carbon": formation["net_hcn_per_feed_carbon"],
            "numerical_reporting_floor_mol": formation["numerical_reporting_floor_mol"],
            "material_ledger_closed": formation["material_ledger_closed"],
            "material_ledger_maximum_absolute_residual_mol": result[
                "material_ledger"
            ]["maximum_absolute_residual_mol"],
            "charged_branch_material_ledger_closed": (
                charged_ledger["closed"] if charged_ledger is not None else None
            ),
            "charged_branch_final_hcn_plus_cyanide_total_mol": charged_cyanide_total,
            "charged_branch_final_hcn_plus_cyanide_aqueous_mol": charged_cyanide_aqueous,
            "charged_to_primary_final_cyanide_ratio": (
                charged_cyanide_total / formation["final_hcn_plus_cyanide_total_mol"]
                if charged_cyanide_total is not None
                and charged_cyanide_total > formation["numerical_reporting_floor_mol"]
                and formation["final_hcn_plus_cyanide_total_mol"]
                > formation["numerical_reporting_floor_mol"]
                else None
            ),
            "charged_branch_material_ledger_maximum_absolute_residual_mol": (
                max(abs(value) for value in charged_ledger["residual"].values())
                if charged_ledger is not None else None
            ),
            "charged_branch_material_ledger_residual_mol": (
                charged_ledger["residual"] if charged_ledger is not None else None
            ),
            "charged_branch_failed_substep_count": charged_branch.get("failed_substep_count"),
            "charged_branch_failed_substep_examples": charged_branch.get("failed_substep_examples"),
            "charged_branch_conservative_ion_clipping_projection_total_mol": charged_branch.get(
                "conservative_ion_clipping_projection_total_mol"
            ),
            "deposited_energy_basis": formation["deposited_energy_basis"],
            "gas_stage_ideal_gas_internal_energy_change_J": result["energy_ledger"][
                "gas_stage_ideal_gas_internal_energy_change_J"
            ],
            "gas_stage_thermochemistry_coverage_complete": result["energy_ledger"][
                "gas_stage_thermochemistry_coverage_complete"
            ],
            "complete_thermal_energy_balance_closed": result["energy_ledger"][
                "complete_thermal_energy_balance_closed"
            ],
        }

    on = cases.get("spark_on")
    full_gate_executed = set(case_names) == set(CASE_NAMES)
    controls_clear = full_gate_executed and all(
        cases[name]["status"] == "no_hcn_above_numerical_floor"
        for name in ("spark_off", "no_methane", "no_nitrogen")
    )
    charged_controls_clear = full_gate_executed and all(
        cases[name]["charged_branch_final_hcn_plus_cyanide_total_mol"] is not None
        and cases[name]["charged_branch_final_hcn_plus_cyanide_total_mol"]
        <= cases[name]["numerical_reporting_floor_mol"]
        for name in ("no_methane", "no_nitrogen")
    )
    charged_spark_on_formed = bool(
        on is not None
        and on["charged_branch_final_hcn_plus_cyanide_total_mol"] is not None
        and on["charged_branch_final_hcn_plus_cyanide_total_mol"]
        > on["numerical_reporting_floor_mol"]
    )
    charged_failures = sum(
        int(item["charged_branch_failed_substep_count"] or 0)
        for item in cases.values()
    )
    return {
        "status": (
            "unseeded_hcn_formation_gate_passed"
            if on is not None
            and on["status"] == "unseeded_hcn_formed"
            and controls_clear
            and charged_controls_clear
            and charged_spark_on_formed
            and all(item["material_ledger_closed"] for item in cases.values())
            and all(
                item["charged_branch_failed_substep_count"] == 0
                for name, item in cases.items() if name != "spark_off"
            )
            else "partial_diagnostic_with_failed_charged_substeps"
            if not full_gate_executed and charged_failures
            else "partial_formation_diagnostic_only" if not full_gate_executed
            else "formation_gate_not_passed"
        ),
        "scenario": "hypothetical fixed-volume CH4/N2/H2/H2O discharge; not a historical apparatus fit",
        "cycles": cycles,
        "substeps_per_cycle": substeps,
        "executed_cases": list(case_names),
        "cases": cases,
        "methane_response": {
            "half_to_full_hcn_ratio": (
                cases["half_methane"]["net_hcn_generated_mol"] / on["net_hcn_generated_mol"]
                if "half_methane" in cases and on is not None
                and on["net_hcn_generated_mol"] > 0 else None
            ),
            "charged_half_to_full_cyanide_ratio": (
                cases["half_methane"]["charged_branch_final_hcn_plus_cyanide_total_mol"]
                / on["charged_branch_final_hcn_plus_cyanide_total_mol"]
                if "half_methane" in cases
                and on is not None
                and charged_spark_on_formed
                and cases["half_methane"]["charged_branch_final_hcn_plus_cyanide_total_mol"]
                is not None else None
            ),
            "experimental_trend_claimed": False,
        },
        "scientific_firewall": {
            "literature_hcn_yields_used_as_solver_inputs": False,
            "hcn_used_to_select_scenario_conditions": False,
            "all_negative_controls_executed": full_gate_executed,
            "spark_on_case_executed": on is not None,
            "charged_branch_negative_controls_below_floor": charged_controls_clear,
            "charged_branch_spark_on_formed_above_floor": charged_spark_on_formed,
            "preexisting_gri_mechanism_contains_hcn_reactions": True,
            "preexisting_pavan_mechanism_contains_hcn_reactions": True,
            "de_novo_product_identity_discovery_claimed": False,
            "modern_experiment_magnitude_validated": False,
            "historical_reproduction_claimed": False,
        },
    }


def summarize_resolution_sweep(reports: list[dict[str, Any]]) -> dict[str, Any]:
    """Compare identical HCN controls at successively finer apparatus steps."""
    if len(reports) < 2:
        raise ValueError("Resolution sweep requires at least two runs")
    steps = [int(report["substeps_per_cycle"]) for report in reports]
    if any(right <= left for left, right in zip(steps, steps[1:])):
        raise ValueError("Resolution sweep substeps must increase strictly")
    cycles = {int(report["cycles"]) for report in reports}
    if len(cycles) != 1:
        raise ValueError("Resolution sweep runs must have identical cycle counts")
    metrics = (
        "net_hcn_generated_mol",
        "final_hcn_aqueous_mol",
        "final_hcn_plus_cyanide_aqueous_mol",
        "charged_branch_final_hcn_plus_cyanide_total_mol",
        "charged_branch_final_hcn_plus_cyanide_aqueous_mol",
    )
    adjacent = []
    for coarse, fine in zip(reports, reports[1:]):
        comparisons = {}
        for metric in metrics:
            if coarse["cases"]["spark_on"].get(metric) is None or fine["cases"]["spark_on"].get(metric) is None:
                raise ValueError(f"Resolution sweep is missing {metric}")
            coarse_value = float(coarse["cases"]["spark_on"][metric])
            fine_value = float(fine["cases"]["spark_on"][metric])
            comparisons[metric] = {
                "coarse_mol": coarse_value,
                "fine_mol": fine_value,
                "relative_change_from_fine": (
                    abs(fine_value - coarse_value) / abs(fine_value)
                    if fine_value != 0 else None
                ),
            }
        adjacent.append({
            "coarse_substeps": coarse["substeps_per_cycle"],
            "fine_substeps": fine["substeps_per_cycle"],
            "spark_on": comparisons,
        })
    return {
        "status": "apparatus_resolution_sweep_completed",
        "cycles": reports[0]["cycles"],
        "substeps_per_cycle": steps,
        "adjacent_refinements": adjacent,
        "all_input_controls_passed": all(
            report["status"] == "unseeded_hcn_formation_gate_passed"
            for report in reports
        ),
        "scientific_firewall": {
            "full_long_duration_convergence_established": False,
            "experimental_magnitude_validated": False,
            "reference_yield_used_as_solver_input": False,
        },
    }


def run_resolution_sweep(
    *, cycles: int = 1, substeps: tuple[int, ...] = (2, 4, 8, 16)
) -> dict[str, Any]:
    return summarize_resolution_sweep([
        run_gate(cycles=cycles, substeps=step) for step in substeps
    ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cycles", type=int, default=1)
    parser.add_argument("--substeps", type=int, default=2)
    parser.add_argument("--resolution-sweep", action="store_true")
    parser.add_argument("--case", action="append", choices=CASE_NAMES,
                        help="Run selected cases; partial runs are diagnostics, not a passed gate")
    arguments = parser.parse_args()
    if arguments.resolution_sweep and arguments.case is not None:
        parser.error("--case cannot be combined with --resolution-sweep")
    report = (
        run_resolution_sweep(cycles=arguments.cycles)
        if arguments.resolution_sweep
        else run_gate(
            cycles=arguments.cycles, substeps=arguments.substeps,
            case_names=tuple(arguments.case) if arguments.case else CASE_NAMES,
        )
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
