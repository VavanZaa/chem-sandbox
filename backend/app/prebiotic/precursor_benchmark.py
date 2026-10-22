"""Held-out precursor energy-yield comparison for discharge experiments."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


DATA_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "benchmarks"
    / "stribling_miller_1986_precursor_energy_v1.json"
)
PROFILE_PATH = DATA_PATH.with_name("stribling_miller_measurement_profile.json")


def _analytical_amount_mol(species: str, inventory: dict[str, float]) -> float:
    """Match the reported cyanide electrode to free HCN plus CN- equivalents."""
    amount = float(inventory.get(species, 0.0))
    return amount + float(inventory.get("CN-", 0.0)) if species == "HCN" else amount


def stribling_miller_measurement_profile() -> dict[str, Any]:
    """Return the versioned timing/feed profile without invented calorimetry."""

    return json.loads(PROFILE_PATH.read_text(encoding="utf-8"))


def _linear_zero_time_yield(
    cycles: list[dict[str, Any]],
    species: str,
    minimum_time_s: float,
    maximum_time_s: float,
    sample_times_s: tuple[float, ...] | None = None,
) -> dict[str, Any]:
    """Extrapolate simulated aqueous energy yield to zero exposure time.

    Stribling and Miller sampled the aqueous phase from one hour to seven days
    and extrapolated energy yield to zero time. The conference record does not
    specify regression weighting, so this implementation exposes an ordinary
    least-squares linear fit and its intercept standard error as a declared
    sensitivity method rather than silently claiming the authors' exact fit.
    """

    rows: list[tuple[float, float]] = []
    cumulative_energy_J = 0.0
    for cycle in cycles:
        cumulative_energy_J += float(
            cycle.get(
                "deposited_energy_used_J",
                cycle.get("deposited_electrical_energy_J", 0.0),
            )
        )
        time_s = float(cycle.get("end_time_s", 0.0))
        if not (minimum_time_s <= time_s <= maximum_time_s) or cumulative_energy_J <= 0:
            continue
        if sample_times_s is not None and not any(
            math.isclose(time_s, sample_time_s, abs_tol=1e-6)
            for sample_time_s in sample_times_s
        ):
            continue
        compartments = cycle.get("compartments_mol", {})
        aqueous_mol = sum(
            _analytical_amount_mol(species, compartments.get(compartment, {}))
            for compartment in ("boiler_liquid", "condensate_trap", "collected_condensate")
        )
        rows.append((time_s, aqueous_mol * 1e9 / cumulative_energy_J))

    covers_window = bool(
        rows
        and rows[0][0] <= minimum_time_s * (1 + 1e-9)
        and rows[-1][0] >= maximum_time_s * (1 - 1e-9)
    )
    if len(rows) < 3:
        return {
            "available": False,
            "covers_published_sampling_window": covers_window,
            "sample_count": len(rows),
            "intercept_nmol_per_J": None,
            "intercept_standard_error_nmol_per_J": None,
            "slope_nmol_per_J_s": None,
            "r_squared": None,
            "method": "ordinary least-squares linear extrapolation",
        }

    mean_time = sum(row[0] for row in rows) / len(rows)
    mean_yield = sum(row[1] for row in rows) / len(rows)
    sxx = sum((row[0] - mean_time) ** 2 for row in rows)
    if sxx <= 0:
        return {
            "available": False,
            "covers_published_sampling_window": covers_window,
            "sample_count": len(rows),
            "intercept_nmol_per_J": None,
            "intercept_standard_error_nmol_per_J": None,
            "slope_nmol_per_J_s": None,
            "r_squared": None,
            "method": "ordinary least-squares linear extrapolation",
        }
    slope = sum(
        (time_s - mean_time) * (value - mean_yield) for time_s, value in rows
    ) / sxx
    intercept = mean_yield - slope * mean_time
    residual_sum = sum(
        (value - (intercept + slope * time_s)) ** 2 for time_s, value in rows
    )
    total_sum = sum((value - mean_yield) ** 2 for _, value in rows)
    variance = residual_sum / (len(rows) - 2)
    intercept_standard_error = math.sqrt(
        max(variance, 0.0) * (1 / len(rows) + mean_time**2 / sxx)
    )
    r_squared = 1.0 - residual_sum / total_sum if total_sum > 0 else 1.0
    return {
        "available": covers_window and math.isfinite(intercept),
        "covers_published_sampling_window": covers_window,
        "sample_count": len(rows),
        "first_sample_time_s": rows[0][0],
        "last_sample_time_s": rows[-1][0],
        "intercept_nmol_per_J": intercept,
        "intercept_standard_error_nmol_per_J": intercept_standard_error,
        "slope_nmol_per_J_s": slope,
        "r_squared": r_squared,
        "method": "ordinary least-squares linear extrapolation",
        "method_boundary": (
            "The public conference record specifies zero-time extrapolation but not "
            "the original regression weighting; this fit is an explicit sensitivity choice."
        ),
    }


def evaluate_precursor_energy_benchmark(
    initial_mol: dict[str, float],
    final_liquid_mol: dict[str, float],
    deposited_energy_J: float,
    *,
    cycles: list[dict[str, Any]] | None = None,
    energy_basis: str = "declared_electrical_estimate",
    calorimetric_evidence: dict[str, Any] | None = None,
    declared_reference_pressure_Pa: float | None = None,
    simulated_gas_pressure_range_Pa: list[float] | None = None,
) -> dict[str, Any]:
    """Report yields while refusing a score for a composition mismatch."""

    benchmark = json.loads(DATA_PATH.read_text(encoding="utf-8"))
    condition = benchmark["reference_condition"]
    protocol_evidence = benchmark.get("reference_protocol_evidence", {})
    required_protocol_fields = (
        "exact_sample_times_and_regression_weighting",
        "absolute_initial_gas_amounts_and_pressure",
        "discharge_flask_geometry_and_temperature_history",
        "calorimetric_energy_trace_or_per_sample_values",
        "aqueous_collection_and_sampling_recovery",
    )
    missing_reference_protocol_fields = [
        field for field in required_protocol_fields
        if protocol_evidence.get(field) is None
    ]
    reference_protocol_complete = not missing_reference_protocol_fields
    # Source completeness alone is not a match: an eventual full protocol must
    # be checked field-by-field against the simulated apparatus and sample
    # schedule before a historical score can be emitted.
    reference_protocol_match_assessed = False
    h2 = float(initial_mol.get("H2", 0.0))
    ch4 = float(initial_mol.get("CH4", 0.0))
    ratio = h2 / ch4 if ch4 > 0 else None
    required_species = {
        species: float(initial_mol.get(species, 0.0)) > 0
        for species in condition["required_species"]
    }
    excluded_species = {
        species: float(initial_mol.get(species, 0.0)) <= 0
        for species in condition["excluded_species"]
    }
    ratio_matched = ratio is not None and math.isclose(
        ratio,
        float(condition["H2_to_CH4_molar_ratio"]),
        rel_tol=0.05,
        abs_tol=0.0,
    )
    composition_matched = (
        all(required_species.values())
        and all(excluded_species.values())
        and ratio_matched
    )
    energy_available = deposited_energy_J > 0
    sampling_window = condition.get("sampling_window_s", [3600.0, 7 * 24 * 3600.0])
    zero_time = {
        species: _linear_zero_time_yield(
            cycles or [], species, float(sampling_window[0]), float(sampling_window[1])
        )
        for species in benchmark["energy_yields_nmol_per_J"]
    }
    # The accessible journal abstract describes periodic measurements, whereas
    # the conference record specifies only the 1 h–7 d span. This sparse grid
    # is an illustrative post-simulation fit sensitivity, never the scoring fit
    # or a substitute for the original paper's unavailable exact sampling log.
    sparse_times_s = (
        float(sampling_window[0]),
        2 * 24 * 3600.0,
        4 * 24 * 3600.0,
        6 * 24 * 3600.0,
        float(sampling_window[1]),
    )
    sparse_zero_time = {
        species: _linear_zero_time_yield(
            cycles or [], species, float(sampling_window[0]),
            float(sampling_window[1]), sparse_times_s
        )
        for species in benchmark["energy_yields_nmol_per_J"]
    }
    zero_time_yield_basis_matched = bool(
        zero_time and all(item["available"] for item in zero_time.values())
    )
    calorimetric_energy_basis_matched = bool(
        energy_basis == "calorimetric_measurement"
        and calorimetric_evidence is not None
        and calorimetric_evidence.get("basis") in {"experimental", "literature"}
    )
    pressure_self_consistent = (
        all(
            abs(value - declared_reference_pressure_Pa)
            <= 0.1 * declared_reference_pressure_Pa
            for value in simulated_gas_pressure_range_Pa
        )
        if declared_reference_pressure_Pa is not None
        and declared_reference_pressure_Pa > 0
        and simulated_gas_pressure_range_Pa
        else None
    )
    comparable = (
        composition_matched
        and energy_available
        and zero_time_yield_basis_matched
        and calorimetric_energy_basis_matched
        and reference_protocol_complete
        and reference_protocol_match_assessed
        and (declared_reference_pressure_Pa is None or pressure_self_consistent is True)
    )
    final_inventory_diagnostic = {
        species: (
            _analytical_amount_mol(species, final_liquid_mol) * 1e9 / deposited_energy_J
            if energy_available
            else None
        )
        for species in benchmark["energy_yields_nmol_per_J"]
    }
    observed = {
        species: (
            zero_time[species]["intercept_nmol_per_J"]
            if zero_time_yield_basis_matched
            else final_inventory_diagnostic[species]
        )
        for species in benchmark["energy_yields_nmol_per_J"]
    }
    comparison = {}
    for species, reference in benchmark["energy_yields_nmol_per_J"].items():
        value = observed[species]
        comparison[species] = {
            "simulated_nmol_per_J": value,
            "observable_basis": (
                "zero_time_extrapolated_aqueous_energy_yield"
                if zero_time_yield_basis_matched
                else "final_retained_aqueous_inventory_diagnostic"
            ),
            "final_inventory_diagnostic_nmol_per_J": final_inventory_diagnostic[species],
            "reference_nmol_per_J": float(reference),
            "ratio_to_reference": value / float(reference) if comparable and value is not None else None,
            "scored": comparable,
        }
    blockers = []
    if not all(required_species.values()):
        missing = [species for species, present in required_species.items() if not present]
        blockers.append(f"Reference feed species missing: {', '.join(missing)}.")
    if not all(excluded_species.values()):
        present = [species for species, absent in excluded_species.items() if not absent]
        blockers.append(f"Reference feed excludes: {', '.join(present)}.")
    if not ratio_matched:
        blockers.append("The reference CH4 case used H2/CH4 = 1 within the published comparison.")
    if not energy_available:
        blockers.append("A positive deposited-energy measurement is required.")
    if not zero_time_yield_basis_matched:
        blockers.append(
            "The published values are zero-time energy yields from samples spanning 1 hour to 7 days; the current run does not supply that complete extrapolation window."
        )
    if not calorimetric_energy_basis_matched:
        blockers.append(
            "The published spark energy was estimated calorimetrically; the current electrical ledger is a declared deposited-energy scenario."
        )
    if not reference_protocol_complete:
        blockers.append(
            "The accessible reference does not specify a complete experimental "
            "protocol for quantitative reproduction: "
            + ", ".join(missing_reference_protocol_fields) + ". "
            "The 1953 apparatus settings are not substitutes for this 1987 experiment."
        )
    if not reference_protocol_match_assessed:
        blockers.append(
            "No field-by-field historical apparatus and sampling-protocol match "
            "has been implemented; a matching energy basis alone cannot establish reproduction."
        )
    if pressure_self_consistent is False:
        blockers.append(
            "The inventory-derived gas pressure differs by more than 10% from the "
            "declared reference pressure; this fixed-volume scenario is internally "
            "inconsistent for benchmark comparison."
        )
    elif declared_reference_pressure_Pa is not None and pressure_self_consistent is None:
        blockers.append(
            "No inventory-derived gas-pressure range was available to check the "
            "declared reference pressure; benchmark comparison remains unavailable."
        )
    if zero_time_yield_basis_matched:
        blockers.append(
            "The public conference record does not specify the original regression weighting or reference uncertainty; the exposed linear intercept cannot by itself establish a validation pass."
        )
    return {
        "status": (
            "comparable"
            if comparable
            else "diagnostic_only_reference_protocol_incomplete"
            if composition_matched and not reference_protocol_complete
            else "diagnostic_only_reference_protocol_unassessed"
            if composition_matched and not reference_protocol_match_assessed
            else "diagnostic_only_measurement_basis_mismatch"
            if composition_matched
            else "diagnostic_only_composition_mismatch"
        ),
        "benchmark_id": benchmark["benchmark_id"],
        "comparable": comparable,
        "composition_matched": composition_matched,
        "deposited_energy_J": deposited_energy_J,
        "energy_basis": energy_basis,
        "calorimetric_evidence": calorimetric_evidence,
        "reference_protocol_complete": reference_protocol_complete,
        "reference_protocol_match_assessed": reference_protocol_match_assessed,
        "missing_reference_protocol_fields": missing_reference_protocol_fields,
        "pressure_self_consistent": pressure_self_consistent,
        "declared_reference_pressure_Pa": declared_reference_pressure_Pa,
        "simulated_gas_pressure_range_Pa": simulated_gas_pressure_range_Pa,
        "zero_time_extrapolation": zero_time,
        "illustrative_sparse_sampling_sensitivity": {
            "status": "post_simulation_method_sensitivity_not_scoring_basis",
            "requested_sample_times_s": list(sparse_times_s),
            "zero_time_extrapolation": sparse_zero_time,
            "limitations": (
                "The exact experimental sample times and regression weighting are "
                "not available in the accessible conference/journal abstracts; "
                "this illustrative grid must not be called the published fit."
            ),
        },
        "feed_checks": {
            "required_species_present": required_species,
            "excluded_species_absent": excluded_species,
            "H2_to_CH4_molar_ratio": ratio,
            "H2_to_CH4_ratio_matched": ratio_matched,
            "zero_time_yield_basis_matched": zero_time_yield_basis_matched,
            "calorimetric_energy_basis_matched": calorimetric_energy_basis_matched,
            "reference_protocol_complete": reference_protocol_complete,
            "reference_protocol_match_assessed": reference_protocol_match_assessed,
        },
        "species_comparison": comparison,
        "remaining_blockers": blockers,
        "sources": {
            "conference": benchmark["conference_source"],
            "metadata": benchmark["metadata_source"],
            "journal": benchmark["journal_source"],
        },
        "scientific_firewall": {
            "reference_yields_used_for_simulation": False,
            "reference_yields_used_for_fitting": False,
            "mismatched_feed_scored": False,
            "final_inventory_equated_with_zero_time_yield": False,
            "declared_electrical_energy_equated_with_calorimetry": False,
            "pressure_inconsistent_scenario_scored": False,
            "incomplete_reference_protocol_scored": False,
            "unassessed_reference_protocol_scored": False,
            "zero_time_fit_used_for_simulation": False,
            "sparse_sampling_sensitivity_used_for_simulation_or_scoring": False,
            "calorimetric_measurement_used_for_simulation": bool(
                energy_basis == "calorimetric_measurement"
            ),
        },
    }
