"""Local model sensitivity to hypothetical spark exposure, not HCN fitting.

The perturbations are dimensionless factors around the existing formation
scenario. They are not measured uncertainty ranges or experimental estimates.
The held-out Barth cyanide data is never loaded here.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
import math
from typing import Any

from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.schema import PrebioticApparatusRequest
from scripts.run_hcn_formation_gate import make_formation_scenario


PERTURBATIONS = (
    ("baseline", "none", 1.0),
    ("half_zone", "discharge_zone_volume", 0.5),
    ("double_zone", "discharge_zone_volume", 2.0),
    ("half_density", "electron_number_density", 0.5),
    ("double_density", "electron_number_density", 2.0),
    ("half_coupling", "electrical_coupling_efficiency", 0.5),
    ("double_coupling", "electrical_coupling_efficiency", 2.0),
)


def make_sensitivity_payload(
    base: dict[str, Any], *, name: str, axis: str, factor: float,
) -> dict[str, Any]:
    if axis not in {"none", "discharge_zone_volume", "electron_number_density", "electrical_coupling_efficiency"}:
        raise ValueError("Unknown spark-exposure sensitivity axis")
    if not math.isfinite(factor) or factor <= 0 or (axis == "none" and factor != 1):
        raise ValueError("Sensitivity factor must be positive and baseline factor one")
    payload = deepcopy(base)
    payload["id"] = f"hcn_exposure_sensitivity_{name}"
    if axis == "electrical_coupling_efficiency":
        payload["spark"]["coupling_efficiency"] *= factor
        if payload["spark"]["coupling_efficiency"] > 1:
            raise ValueError("Electrical coupling sensitivity exceeds physical fraction one")
        payload["spark"]["source"] = (
            f"Hypothetical {factor:g}x electrical-coupling sensitivity; no HCN fit"
        )
        payload["spark"]["basis"] = "hypothetical"
        return payload
    if axis == "discharge_zone_volume":
        item = payload["apparatus"]["discharge_zone_volume"]
    elif axis == "electron_number_density":
        item = payload["spark"]["boltzmann"]["electron_number_density"]
    else:
        return payload
    item["value"] *= factor
    item["source"] = f"Hypothetical {factor:g}x local sensitivity of the declared formation scenario"
    item["method"] = "dimensionless parameter perturbation; not experimentally constrained"
    item["basis"] = "hypothetical"
    return payload


def _charged_cyanide_total(result: dict[str, Any]) -> float | None:
    branch = result["plasma_model"].get("propagated_charged_branch")
    if branch is None:
        return None
    return sum(
        float(species.get("HCN", 0.0)) + float(species.get("CN-", 0.0))
        for species in branch["final_compartments_mol"].values()
    )


def _log_elasticity(low: float, high: float) -> float | None:
    return math.log(high / low) / math.log(4.0) if low > 0 and high > 0 else None


def run_exposure_sensitivity(
    *, cycles: int = 1, substeps: int = 2,
    perturbations: tuple[tuple[str, str, float], ...] = PERTURBATIONS,
) -> dict[str, Any]:
    if not perturbations or len({name for name, _, _ in perturbations}) != len(perturbations):
        raise ValueError("Sensitivity cases must be nonempty and have unique names")
    base = make_formation_scenario(cycles=cycles, substeps=substeps)
    cases: dict[str, dict[str, Any]] = {}
    for name, axis, factor in perturbations:
        payload = make_sensitivity_payload(base, name=name, axis=axis, factor=factor)
        result = simulate_prebiotic_apparatus(PrebioticApparatusRequest.model_validate(payload))
        branch = result["plasma_model"].get("propagated_charged_branch") or {}
        plasma_steps = result["plasma_model"]["cycle_summaries"]
        cases[name] = {
            "axis": axis,
            "factor_from_declared_baseline": factor,
            "discharge_zone_volume_L": payload["apparatus"]["discharge_zone_volume"]["value"],
            "declared_electron_density_m3": payload["spark"]["boltzmann"]["electron_number_density"]["value"],
            "electrical_coupling_efficiency": payload["spark"]["coupling_efficiency"],
            "net_hcn_plus_cyanide_generated_mol": result["hcn_formation"]["net_hcn_generated_mol"],
            "final_aqueous_hcn_plus_cyanide_mol": result["hcn_formation"]["final_hcn_plus_cyanide_aqueous_mol"],
            "charged_branch_final_hcn_plus_cyanide_total_mol": _charged_cyanide_total(result),
            "declared_electrical_energy_J": result["energy_ledger"]["deposited_electrical_energy_J"],
            "deposited_energy_basis": result["energy_ledger"]["deposited_energy_basis"],
            "effective_electron_density_range_m3": [
                min(step["electron_number_density_m3"] for step in plasma_steps),
                max(step["electron_number_density_m3"] for step in plasma_steps),
            ] if plasma_steps else None,
            "declared_electron_density_power_capped": bool(plasma_steps) and all(
                step["electron_number_density_m3"]
                < payload["spark"]["boltzmann"]["electron_number_density"]["value"]
                * (1 - 1e-9)
                for step in plasma_steps
            ),
            "material_ledger_closed": result["material_ledger"]["closed"],
            "charged_branch_material_ledger_closed": (branch.get("material_ledger") or {}).get("closed"),
            "charged_branch_failed_substeps": branch.get("failed_substep_count"),
            "complete_thermal_energy_balance_closed": result["energy_ledger"]["complete_thermal_energy_balance_closed"],
        }
        if not cases[name]["material_ledger_closed"] or not cases[name]["charged_branch_material_ledger_closed"] or cases[name]["charged_branch_failed_substeps"] != 0:
            raise RuntimeError(f"Spark-exposure sensitivity case {name} failed material/numerical controls")
    elasticities = {}
    for axis, lower, upper in (
        ("discharge_zone_volume", "half_zone", "double_zone"),
        ("electron_number_density", "half_density", "double_density"),
        ("electrical_coupling_efficiency", "half_coupling", "double_coupling"),
    ):
        if lower in cases and upper in cases:
            elasticities[axis] = {
                "final_aqueous_hcn_plus_cyanide": _log_elasticity(
                    cases[lower]["final_aqueous_hcn_plus_cyanide_mol"],
                    cases[upper]["final_aqueous_hcn_plus_cyanide_mol"],
                ),
                "charged_branch_final_hcn_plus_cyanide_total": _log_elasticity(
                    cases[lower]["charged_branch_final_hcn_plus_cyanide_total_mol"],
                    cases[upper]["charged_branch_final_hcn_plus_cyanide_total_mol"],
                ) if cases[lower]["charged_branch_final_hcn_plus_cyanide_total_mol"] is not None
                and cases[upper]["charged_branch_final_hcn_plus_cyanide_total_mol"] is not None else None,
            }
    return {
        "status": "within_model_sensitivity_diagnostic_only",
        "cycles": cycles,
        "substeps_per_cycle": substeps,
        "cases": cases,
        "central_log_elasticities": elasticities,
        "scientific_firewall": {
            "experimental_hcn_loaded": False,
            "hcn_used_to_select_exposure": False,
            "parameter_factors_are_measured_uncertainties": False,
            "absolute_yield_validated": False,
            "de_novo_product_discovery_claimed": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cycles", type=int, default=1)
    parser.add_argument("--substeps", type=int, default=2)
    parser.add_argument(
        "--axis",
        choices=("discharge_zone_volume", "electron_number_density", "electrical_coupling_efficiency"),
        help="Run only the baseline and the selected pair of perturbations",
    )
    arguments = parser.parse_args()
    perturbations = (
        tuple(item for item in PERTURBATIONS if item[1] in {"none", arguments.axis})
        if arguments.axis else PERTURBATIONS
    )
    print(json.dumps(run_exposure_sensitivity(
        cycles=arguments.cycles, substeps=arguments.substeps,
        perturbations=perturbations,
    ), indent=2))


if __name__ == "__main__":
    main()
