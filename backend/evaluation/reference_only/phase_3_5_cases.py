"""Reference-only checks loaded after Phase 3-5 predictions are frozen."""

from __future__ import annotations

import math
from typing import Any


R_J_MOL_K = 8.31446261815324


def evaluate(prediction: dict[str, Any]) -> dict[str, Any]:
    temperature = 500.0
    expected_rate_ratio = math.exp(10_000.0 / (R_J_MOL_K * temperature))
    rate_ratio = prediction["rate"]["ratio"]
    rate_error = abs(rate_ratio - expected_rate_ratio) / expected_rate_ratio

    hcp = 3.3e-4
    liquid_volume = 1e-3
    gas_volume = 1e-3
    total_amount = 1e-3
    gas_factor = R_J_MOL_K * 298.15 / gas_volume
    expected_dissolved = hcp * liquid_volume * gas_factor * total_amount / (
        1.0 + hcp * liquid_volume * gas_factor
    )
    dissolved = prediction["gas_liquid"]["final_dissolved_mol"]
    gas_error = abs(dissolved - expected_dissolved) / expected_dissolved

    cantera_checks = {}
    for reactor, result in prediction["cantera"].items():
        passed = bool(
            result["accepted"]
            and result["conversion"] > 0.1
            and result["water_yield"] > 0.1
            and result["final_temperature_K"] > 1200.0
            and result["material_residual"] < 1e-8
            and result["energy_relative_residual"] < result["energy_relative_tolerance"]
        )
        cantera_checks[reactor] = {
            "passed": passed,
            "criterion": "reactive H2/O2 case ignites, forms water, and closes element/energy balances",
        }

    reaktoro = prediction.get("reaktoro")
    reaktoro_check = {
        "passed": bool(reaktoro and reaktoro.get("accepted") and abs(reaktoro["pH"] - 5.0) < 1e-8),
        "fixed_pH_target": 5.0,
        "absolute_error": None if not reaktoro else abs(reaktoro["pH"] - 5.0),
    }
    checks = {
        "eyring_detailed_balance": {"passed": rate_error < 1e-12, "relative_error": rate_error},
        "closed_headspace_analytic_equilibrium": {"passed": gas_error < 1e-8, "relative_error": gas_error},
        "cantera_reactive_operations": cantera_checks,
        "reaktoro_fixed_pH": reaktoro_check,
    }
    flattened = [checks["eyring_detailed_balance"]["passed"], checks["closed_headspace_analytic_equilibrium"]["passed"]]
    flattened.extend(item["passed"] for item in cantera_checks.values())
    flattened.append(reaktoro_check["passed"])
    return {
        "suite": "phase-3-5-chemical-validation-v1",
        "passed": all(flattened),
        "passed_checks": sum(bool(item) for item in flattened),
        "total_checks": len(flattened),
        "checks": checks,
        "reference_note": (
            "References are analytic identities or fixed test constraints and were imported only after "
            "the raw prediction artifact was frozen. Cantera checks validate reactive behavior and closure, "
            "not experimental mechanism accuracy."
        ),
    }

