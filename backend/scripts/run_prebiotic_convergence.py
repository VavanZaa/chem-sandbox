"""Run a declared Miller-Urey apparatus profile at one integration resolution.

Prints compact, machine-readable diagnostics. It does not use benchmark product
amounts to configure the solver or write a discovery artifact.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from app.prebiotic.apparatus import simulate_prebiotic_apparatus
from app.prebiotic.boltzmann import (
    COMMUNITY_CH4_DATASET_ID,
    IST_H2_COMMUNITY_CH4_DATASET_ID,
    REACTIVE_N2_DATASET_ID,
)
from app.prebiotic.schema import PrebioticApparatusRequest


PROFILE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "benchmarks"
    / "stribling_miller_measurement_profile.json"
)


def water_amount(compartments: dict, species: str) -> float:
    return sum(
        float(compartments.get(place, {}).get(species, 0.0))
        for place in ("boiler_liquid", "condensate_trap", "collected_condensate")
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cycles", type=int, default=24)
    parser.add_argument("--substeps", type=int, required=True)
    parser.add_argument("--forcing", choices=("impulse", "time_averaged"), default="impulse")
    parser.add_argument("--nonwater-return", type=float, default=None)
    parser.add_argument(
        "--phase-transfer",
        choices=("step_equilibrium", "flow_limited_capture", "coupled_flow_transport"),
        default=None,
    )
    parser.add_argument("--contact-holdup-L", type=float, default=None)
    parser.add_argument("--discharge-temperature-K", type=float, default=None)
    parser.add_argument("--reduced-field-Td", type=float, default=None)
    parser.add_argument(
        "--collision-set",
        choices=(REACTIVE_N2_DATASET_ID, COMMUNITY_CH4_DATASET_ID, IST_H2_COMMUNITY_CH4_DATASET_ID),
        default=None,
        help="Optional collision-set sensitivity; the profile default remains unchanged.",
    )
    arguments = parser.parse_args()
    if not (1 <= arguments.cycles <= 168):
        parser.error("cycles must be between 1 and 168")
    if not (1 <= arguments.substeps <= 48):
        parser.error("substeps must be between 1 and 48")

    payload = json.loads(PROFILE.read_text(encoding="utf-8"))
    payload["apparatus"]["cycles"] = arguments.cycles
    payload["apparatus"]["operation_mode"] = (
        "staged" if arguments.substeps == 1 else "continuous_overlap"
    )
    payload["apparatus"]["integration_substeps_per_cycle"] = arguments.substeps
    payload["apparatus"]["primary_source_forcing"] = arguments.forcing
    if arguments.collision_set is not None:
        payload["spark"]["boltzmann"]["collision_set"] = arguments.collision_set
    if arguments.nonwater_return is not None:
        payload["apparatus"]["condensate_nonwater_return_fraction"] = arguments.nonwater_return
    if arguments.phase_transfer is not None:
        payload["apparatus"]["condenser_phase_transfer_model"] = arguments.phase_transfer
    if arguments.contact_holdup_L is not None:
        payload["apparatus"]["condenser_contact_water_holdup"] = {
            "value": arguments.contact_holdup_L,
            "unit": "L",
            "source": "hypothetical upper-bound condenser contact-water sensitivity",
            "method": "declared scenario parameter; not measured apparatus hold-up",
            "basis": "hypothetical",
        }
    if arguments.discharge_temperature_K is not None:
        payload["apparatus"]["discharge_temperature"]["value"] = (
            arguments.discharge_temperature_K
        )
    if arguments.reduced_field_Td is not None:
        payload["spark"]["boltzmann"]["reduced_electric_field_Td"] = (
            arguments.reduced_field_Td
        )
    result = simulate_prebiotic_apparatus(
        PrebioticApparatusRequest.model_validate(payload)
    )
    branch = result["plasma_model"]["propagated_charged_branch"]
    if branch is None:
        raise RuntimeError("Expected the propagated charged branch")
    def water_pair(compartments: dict) -> dict[str, float]:
        return {
            "HCN": water_amount(compartments, "HCN"),
            "CH2O": water_amount(compartments, "CH2O"),
        }
    def fits(report: dict) -> dict:
        return {
            species: {
                key: report["zero_time_extrapolation"][species][key]
                for key in (
                    "available",
                    "sample_count",
                    "intercept_nmol_per_J",
                    "intercept_standard_error_nmol_per_J",
                    "r_squared",
                )
            }
            for species in ("HCN", "CH2O")
        }
    summary = {
        "cycles": arguments.cycles,
        "collision_set": payload["spark"]["boltzmann"]["collision_set"],
        "discharge_temperature_K": payload["apparatus"]["discharge_temperature"]["value"],
        "reduced_electric_field_Td": payload["spark"]["boltzmann"]["reduced_electric_field_Td"],
        "substeps": arguments.substeps,
        "operation_mode": result["operating_conditions"]["operation_mode"],
        "primary_source_forcing": result["operating_conditions"]["primary_source_forcing"],
        "condensate_nonwater_return_fraction": result["operating_conditions"]["condensate_nonwater_return_fraction"],
        "condenser_phase_transfer_model": result["operating_conditions"]["condenser_phase_transfer_model"],
        "condenser_contact_water_holdup_L": result["operating_conditions"]["condenser_contact_water_holdup_L"],
        "integration_steps": len(result["plasma_model"]["cycle_summaries"]),
        "total_duration_s": result["total_duration_s"],
        "deposited_energy_J": result["energy_ledger"]["deposited_energy_used_J"],
        "primary_ledger_closed": result["material_ledger"]["closed"],
        "charged_ledger_closed": branch["material_ledger"]["closed"],
        "primary_benchmark_comparable": result["precursor_energy_benchmark"]["comparable"],
        "charged_benchmark_comparable": branch["precursor_energy_benchmark"]["comparable"],
        "pressure_self_consistent": result["precursor_energy_benchmark"][
            "pressure_self_consistent"
        ],
        "charged_pressure_self_consistent": branch["precursor_energy_benchmark"][
            "pressure_self_consistent"
        ],
        "charged_simulated_gas_pressure_range_Pa": branch[
            "simulated_gas_pressure_range_Pa"
        ],
        "primary_final_water_mol": water_pair(result["final_compartments_mol"]),
        "charged_final_water_mol": water_pair(branch["final_compartments_mol"]),
        "first_hour_primary_water_mol": water_pair(result["cycles"][0]["compartments_mol"]),
        "first_hour_charged_water_mol": water_pair(branch["cycle_summaries"][0]["compartments_mol"]),
        "primary_hcn_zero_time_nmol_J": result["precursor_energy_benchmark"]
        ["zero_time_extrapolation"]["HCN"]["intercept_nmol_per_J"],
        "charged_hcn_zero_time_nmol_J": branch["precursor_energy_benchmark"]
        ["zero_time_extrapolation"]["HCN"]["intercept_nmol_per_J"],
        "primary_zero_time_fits": fits(result["precursor_energy_benchmark"]),
        "charged_zero_time_fits": fits(branch["precursor_energy_benchmark"]),
        "primary_illustrative_sparse_zero_time_fits": fits(
            result["precursor_energy_benchmark"]["illustrative_sparse_sampling_sensitivity"]
        ),
        "charged_illustrative_sparse_zero_time_fits": fits(
            branch["precursor_energy_benchmark"]["illustrative_sparse_sampling_sensitivity"]
        ),
        "simulated_gas_pressure_range_Pa": result["operating_conditions"][
            "simulated_gas_pressure_range_Pa"
        ],
        "charged_max_ledger_residual_mol": max(
            (abs(float(value)) for value in branch["material_ledger"]["residual"].values()),
            default=0.0,
        ),
    }
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
