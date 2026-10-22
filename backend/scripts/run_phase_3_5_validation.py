from __future__ import annotations

import json
import os
from pathlib import Path

from app.multiphase.gas_liquid import simulate_gas_liquid_transfer
from app.multiphase.reaktoro_bridge import equilibrate_with_reaktoro
from app.multiphase.schema import ElectrolyteEquilibriumRequest, GasLiquidTransferRequest
from app.operations.cantera_runner import simulate_cantera
from app.operations.cantera_schema import CanteraSimulationRequest
from app.pathways.rates import build_rate_artifact
from app.pathways.schema import RateArtifactRequest
from evaluation.protocol import freeze_prediction_artifact, unlock_evaluation
from probe_phase_3_5 import cantera_payload, energy


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "backend" / "data" / "validation" / "phase_3_5_latest"


def jsonable(value):
    return json.loads(json.dumps(value, default=lambda item: item.item()))


def predictions() -> dict:
    rate = build_rate_artifact(RateArtifactRequest.model_validate({
        "reaction_id": "validation_A_to_B",
        "temperature": {"value": 500.0, "unit": "K"},
        "pressure": {"value": 1.0, "unit": "atm"},
        "forward_order": 1,
        "reverse_order": 1,
        "forward_activation_gibbs": energy(80.0),
        "reverse_activation_gibbs": energy(90.0),
        "reaction_gibbs": energy(-10.0),
    }))
    cantera = {}
    for reactor in ("batch", "pfr", "cstr"):
        result = simulate_cantera(CanteraSimulationRequest.model_validate(cantera_payload(reactor)))
        energy_balance = result["energy_balance"]
        cantera[reactor] = {
            "accepted": bool(result["accepted"]),
            "conversion": result["metrics"]["conversion"],
            "water_yield": result["metrics"]["targets"]["H2O"]["yield_fraction"],
            "final_temperature_K": result["samples"][-1]["temperature_K"],
            "material_residual": result["material_balance"]["max_relative_element_residual"],
            "energy_relative_residual": energy_balance.get("relative_residual", 0.0),
            "energy_relative_tolerance": energy_balance.get("relative_tolerance", 0.0),
            "mechanism_sha256": result["provenance"]["mechanism_sha256"],
        }

    gas_request = GasLiquidTransferRequest.model_validate({
        "solute_id": "validation-solute",
        "mode": "closed_headspace",
        "temperature": {"value": 298.15, "unit": "K"},
        "total_pressure": {"value": 2.0, "unit": "atm"},
        "gas_volume": {"value": 1.0, "unit": "liter"},
        "liquid_volume": {"value": 1.0, "unit": "liter"},
        "duration": {"value": 1000.0, "unit": "s"},
        "points": 11,
        "initial_gas_amount": {"value": 0.001, "unit": "mol"},
        "initial_dissolved_amount": {"value": 0.0, "unit": "mol"},
        "henry_hcp": {
            "value": 3.3e-4,
            "unit": "mol/(m^3*Pa)",
            "source": "validation-only measured coefficient fixture",
            "method": "declared Hcp convention",
            "basis": "experimental",
        },
        "kla": {
            "value": 0.1,
            "unit": "1/s",
            "source": "validation-only transfer coefficient fixture",
            "method": "declared dynamic transfer coefficient",
            "basis": "experimental",
        },
    })
    gas = simulate_gas_liquid_transfer(gas_request)

    os.environ.setdefault("CHEM_SANDBOX_REAKTORO_URL", "http://127.0.0.1:8011")
    reaktoro = equilibrate_with_reaktoro(ElectrolyteEquilibriumRequest.model_validate({
        "database_kind": "phreeqc",
        "database_name": "pitzer.dat",
        "aqueous_species": ["H2O", "H+", "OH-", "Na+", "Cl-"],
        "activity_model": "pitzer",
        "temperature": {"value": 298.15, "unit": "K"},
        "pressure": {"value": 1.0, "unit": "atm"},
        "components": [
            {"species": "H2O", "amount": {"value": 1.0, "unit": "kg"}},
            {"species": "Na+", "amount": {"value": 0.1, "unit": "mol"}},
            {"species": "Cl-", "amount": {"value": 0.1, "unit": "mol"}},
        ],
        "fixed_pH": 5.0,
        "source": "post-freeze fixed-pH equilibrium validation request",
    }))
    return jsonable({
        "suite": "phase-3-5-chemical-validation-v1",
        "rate": {
            "certified": rate["certified"],
            "ratio": rate["forward"]["rate_constant"] / rate["reverse"]["rate_constant"],
            "equilibrium_constant": rate["thermodynamics"]["equilibrium_constant_dimensionless"],
        },
        "cantera": cantera,
        "gas_liquid": {
            "accepted": gas["accepted"],
            "final_dissolved_mol": gas["samples"][-1]["dissolved_mol"],
            "mass_balance_residual_mol": gas["mass_balance"]["residual_mol"],
            "equilibrium_residual_mol_m3": gas["final_equilibrium_residual_mol_m3"],
        },
        "reaktoro": {
            "accepted": reaktoro["accepted"],
            "pH": reaktoro["aqueous"]["pH"],
            "request_sha256": reaktoro["request_sha256"],
            "version": reaktoro["solver"]["version"],
        },
    })


def main() -> None:
    raw = predictions()
    frozen_path = OUTPUT / "predictions.frozen.json"
    frozen = freeze_prediction_artifact(raw, frozen_path)
    # This import is intentionally after the immutable prediction write.
    from evaluation.reference_only.phase_3_5_cases import evaluate

    evaluation = evaluate(raw)
    unlocked_path = OUTPUT / "evaluation.unlocked.json"
    unlock_evaluation(frozen_path, evaluation, unlocked_path)
    report = [
        "# Phase 3-5 chemical validation",
        "",
        f"- Overall: {'PASS' if evaluation['passed'] else 'FAIL'}",
        f"- Checks: {evaluation['passed_checks']}/{evaluation['total_checks']}",
        f"- Frozen prediction SHA-256: `{frozen['prediction_sha256']}`",
        "- Reference access: unlocked only after prediction freeze",
        "",
        "## Checks",
        "",
    ]
    for name, result in evaluation["checks"].items():
        if name == "cantera_reactive_operations":
            for reactor, nested in result.items():
                report.append(f"- Cantera {reactor}: {'PASS' if nested['passed'] else 'FAIL'}")
        else:
            report.append(f"- {name}: {'PASS' if result['passed'] else 'FAIL'}")
    report.extend([
        "",
        "The reactive H2/O2 cases validate the operation wrapper and conservation against the bundled mechanism; ",
        "they are not a claim that the mechanism predicts experimental yields for every apparatus.",
        "",
    ])
    (OUTPUT / "REPORT.md").write_text("\n".join(report), encoding="utf-8")
    print(json.dumps({
        "passed": evaluation["passed"],
        "checks": f"{evaluation['passed_checks']}/{evaluation['total_checks']}",
        "prediction_sha256": frozen["prediction_sha256"],
        "report": str(OUTPUT / "REPORT.md"),
    }, indent=2))
    if not evaluation["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
