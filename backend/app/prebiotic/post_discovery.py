"""One-way, evidence-labelled evaluation of a frozen aqueous discovery artifact."""
from __future__ import annotations

from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq


DATA_PATH = (
    Path(__file__).resolve().parents[2]
    / "data"
    / "kinetics"
    / "post_discovery_strecker_v1.json"
)
WATER_MOLAR_MASS_G = 18.01528
AMMONIA_MOLAR_MASS_G = 17.03052


def freeze_discovery_artifact(discovery: dict[str, Any]) -> tuple[str, bytes]:
    serialized = json.dumps(
        discovery,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return sha256(serialized).hexdigest(), serialized


def _equilibrium_extent(
    carbonyl_M: float,
    ammonia_M: float,
    hcn_M: float,
    equilibrium_constant_M_minus_2: float,
) -> float:
    limit = min(carbonyl_M, ammonia_M, hcn_M)
    if limit <= 0:
        return 0.0

    def residual(extent: float) -> float:
        return extent - equilibrium_constant_M_minus_2 * (
            carbonyl_M - extent
        ) * (ammonia_M - extent) * (hcn_M - extent)

    upper = limit * (1 - 1e-12)
    if residual(upper) <= 0:
        return upper
    return float(brentq(residual, 0.0, upper, xtol=1e-18, rtol=1e-12))


def _equilibrium_increment(
    carbonyl_mol: float,
    ammonia_mol: float,
    hcn_mol: float,
    aminonitrile_mol: float,
    liquid_volume_L: float,
    equilibrium_constant_M_minus_2: float,
) -> float:
    """Positive extent needed to restore the cited equilibrium at one instant."""

    limit = min(carbonyl_mol, ammonia_mol, hcn_mol)
    if limit <= 0:
        return 0.0

    def residual(extent_mol: float) -> float:
        product_M = (aminonitrile_mol + extent_mol) / liquid_volume_L
        reactant_product_M3 = (
            (carbonyl_mol - extent_mol)
            * (ammonia_mol - extent_mol)
            * (hcn_mol - extent_mol)
            / liquid_volume_L**3
        )
        return product_M - equilibrium_constant_M_minus_2 * reactant_product_M3

    if residual(0.0) >= 0:
        # This reference-only projection does not model reverse decomposition.
        return 0.0
    upper = limit * (1 - 1e-12)
    if residual(upper) <= 0:
        return upper
    return float(brentq(residual, 0.0, upper, xtol=1e-18, rtol=1e-12))


def _replay_aqueous_cycles(
    cycle_snapshots: list[dict[str, Any]],
    liquid_volume_L: float,
    equilibrium_constant_M_minus_2: float,
    hydrolysis_rate_s_inv: float,
) -> dict[str, Any]:
    """One-way sidecar replay after discovery is frozen.

    Gas production and partitioning remain untouched.  Within each stored
    apparatus cycle, precursor supply is linearly interpolated, rapid
    aminonitrile equilibration is applied, and hydrolysis is integrated in
    twenty substeps.  This makes the projection cycle resolved without making
    it a validated kinetic mechanism.
    """

    intermediate_mol = 0.0
    product_mol = 0.0
    formed_total_mol = 0.0
    trace: list[dict[str, float]] = []
    previous = {"CH2O": 0.0, "NH3": 0.0, "HCN": 0.0}
    substeps = 20
    for cycle in cycle_snapshots:
        compartments = cycle.get("compartments_mol", {})
        boiler = compartments.get("boiler_liquid", {})
        trap = compartments.get("condensate_trap", {})
        current = {
            species: float(boiler.get(species, 0.0)) + float(trap.get(species, 0.0))
            for species in previous
        }
        start_s = float(cycle.get("start_time_s", 0.0))
        end_s = float(cycle.get("end_time_s", start_s))
        dt = max(0.0, end_s - start_s) / substeps
        for step in range(1, substeps + 1):
            fraction = step / substeps
            supplied = {
                species: previous[species] + fraction * (current[species] - previous[species])
                for species in previous
            }
            available = {
                species: max(0.0, supplied[species] - formed_total_mol)
                for species in previous
            }
            formed = _equilibrium_increment(
                available["CH2O"],
                available["NH3"],
                available["HCN"],
                intermediate_mol,
                liquid_volume_L,
                equilibrium_constant_M_minus_2,
            )
            intermediate_mol += formed
            formed_total_mol += formed
            hydrolyzed = intermediate_mol * (
                1 - math.exp(-hydrolysis_rate_s_inv * dt)
            )
            intermediate_mol -= hydrolyzed
            product_mol += hydrolyzed
        trace.append(
            {
                "cycle": float(cycle.get("cycle", len(trace) + 1)),
                "end_time_s": end_s,
                "aminonitrile_remaining_mol": intermediate_mol,
                "cumulative_aminonitrile_formed_mol": formed_total_mol,
                "candidate_product_mol": product_mol,
            }
        )
        previous = current
    return {
        "rate_s_inv": hydrolysis_rate_s_inv,
        "converted_fraction": (
            product_mol / formed_total_mol if formed_total_mol > 0 else 0.0
        ),
        "candidate_product_mol": product_mol,
        "aminonitrile_remaining_mol": intermediate_mol,
        "cumulative_aminonitrile_formed_mol": formed_total_mol,
        "cycle_trace": trace,
    }


def evaluate_frozen_aqueous_discovery(
    discovery: dict[str, Any],
    liquid_mol: np.ndarray,
    prepared: dict[str, Any],
    liquid_volume_L: float,
    duration_s: float,
    liquid_temperature_range_K: tuple[float, float],
    cycle_snapshots: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Attach reference kinetics after discovery; never feed output back to it."""

    if not math.isfinite(liquid_volume_L) or liquid_volume_L <= 0:
        raise ValueError("Post-discovery evaluation requires a positive liquid volume")

    artifact_sha256, serialized = freeze_discovery_artifact(discovery)
    # The evidence file is intentionally opened only after immutable serialization.
    evidence = json.loads(DATA_PATH.read_text(encoding="utf-8"))
    discovered_smiles = {
        candidate["smiles"]
        for record in discovery.get("formula_records", [])
        for candidate in record.get("candidates", [])
    }
    index = prepared["index"]

    def amount(species: str) -> float:
        return float(liquid_mol[index[species]]) if species in index else 0.0

    water_g = amount("H2O") * WATER_MOLAR_MASS_G
    ammonia_g = amount("NH3") * AMMONIA_MOLAR_MASS_G
    ammonia_weight_percent = (
        100 * ammonia_g / (water_g + ammonia_g) if water_g + ammonia_g > 0 else 0.0
    )
    formation = evidence["formation_equilibria"]["NCC#N"]
    aminonitrile_M = _equilibrium_extent(
        amount("CH2O") / liquid_volume_L,
        amount("NH3") / liquid_volume_L,
        amount("HCN") / liquid_volume_L,
        float(formation["equilibrium_constant_M_minus_2"]),
    )
    aminonitrile_mol = aminonitrile_M * liquid_volume_L
    hydrolysis = evidence["hydrolysis_sensitivity"]["NCC#N"]
    nominal_rate = float(hydrolysis["rate_s_inv"])
    uncertainty = float(hydrolysis["standard_uncertainty_s_inv"])
    rates = {
        "low_1sigma": max(0.0, nominal_rate - uncertainty),
        "nominal": nominal_rate,
        "high_1sigma": nominal_rate + uncertainty,
    }
    if cycle_snapshots:
        sensitivity = {
            label: _replay_aqueous_cycles(
                cycle_snapshots,
                liquid_volume_L,
                float(formation["equilibrium_constant_M_minus_2"]),
                rate,
            )
            for label, rate in rates.items()
        }
        projection_method = "cycle_resolved_post_discovery_sidecar"
    else:
        sensitivity = {
            label: {
                "rate_s_inv": rate,
                "converted_fraction": 1 - math.exp(-rate * duration_s),
                "candidate_product_mol": aminonitrile_mol
                * (1 - math.exp(-rate * duration_s)),
                "aminonitrile_remaining_mol": aminonitrile_mol
                * math.exp(-rate * duration_s),
                "cumulative_aminonitrile_formed_mol": aminonitrile_mol,
                "cycle_trace": [],
            }
            for label, rate in rates.items()
        }
        projection_method = "single_interval_post_discovery_sensitivity"
    temperature_matches = (
        liquid_temperature_range_K[0]
        <= float(hydrolysis["temperature_K"])
        <= liquid_temperature_range_K[1]
    )
    ammonia_matches = math.isclose(
        ammonia_weight_percent,
        float(hydrolysis["ammonia_weight_percent"]),
        rel_tol=0.1,
        abs_tol=0.5,
    )
    product_discovered = hydrolysis["candidate_product_smiles"] in discovered_smiles
    applicability = (
        "within_reported_conditions"
        if product_discovered and temperature_matches and ammonia_matches
        else "sensitivity_only_conditions_outside_evidence"
    )
    return {
        "status": "post_discovery_reference_evaluation",
        "discovery_artifact_sha256": artifact_sha256,
        "frozen_artifact_bytes": len(serialized),
        "candidate_product": {
            "name": hydrolysis["candidate_product_name"] if product_discovered else None,
            "smiles": hydrolysis["candidate_product_smiles"],
            "present_in_frozen_discovery": product_discovered,
        },
        "reference_equilibrium": {
            "intermediate_name": formation["candidate_name"],
            "intermediate_smiles": "NCC#N",
            "estimated_equilibrium_mol_at_298_15_K": aminonitrile_mol,
            "equilibrium_constant_M_minus_2": formation["equilibrium_constant_M_minus_2"],
            "aqueous_volume_L": liquid_volume_L,
            "aqueous_volume_basis": "final water inventory at 18.01528 mL/mol",
            "source": formation["source"],
            "secondary_value_source": formation["secondary_value_source"],
        },
        "hydrolysis_sensitivity": {
            "duration_s": duration_s,
            "projection_method": projection_method,
            "cycles_replayed": len(cycle_snapshots or []),
            "scenarios": sensitivity,
            "evidence_temperature_K": hydrolysis["temperature_K"],
            "evidence_ammonia_weight_percent": hydrolysis["ammonia_weight_percent"],
            "apparatus_liquid_temperature_range_K": list(liquid_temperature_range_K),
            "apparatus_final_ammonia_weight_percent": ammonia_weight_percent,
            "source": hydrolysis["source"],
            "applicability": applicability,
        },
        "limitations": [
            "The equilibrium projection is evaluated at the literature reference temperature, not across the apparatus temperature cycle.",
            "Cycle replay is a one-way post-discovery sidecar: rapid equilibrium and linear precursor supply within each stored cycle are approximations.",
            "Hydrolysis amounts are sensitivity scenarios unless both temperature and ammonia concentration match the cited experiment.",
            "The reference mechanism is loaded only after discovery is frozen and cannot add, remove, or rerank discovery candidates.",
            "No historical Miller-Urey product amount or yield is used in this calculation.",
        ],
        "scientific_firewall": {
            "discovery_frozen_before_evidence_load": True,
            "evaluation_feedback_to_discovery": False,
            "historical_benchmark_yields_used": False,
            "known_product_identity_available_only_post_discovery": True,
            "quantitative_yield_claim_allowed": applicability == "within_reported_conditions",
        },
    }
