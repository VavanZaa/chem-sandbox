"""Post-run HCN hydrolysis screening from independently measured rates.

The screen consumes frozen aqueous endpoint snapshots. It does not alter the
prediction, invent formamide/formate yields, or imply a closed thermal ledger.
"""

from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Any

from scipy.constants import R


DATA_PATH = (
    Path(__file__).resolve().parents[2]
    / "data" / "kinetics" / "hcn_hydrolysis_murali_2026_v1.json"
)
_AQUEOUS = ("boiler_liquid", "condensate_trap", "collected_condensate")
_MAX_UNAVAILABLE_DETAILS = 20


@lru_cache(maxsize=1)
def load_hcn_hydrolysis_rate_card() -> dict[str, Any]:
    return json.loads(DATA_PATH.read_text(encoding="utf-8"))


def _rate_at_temperature(item: dict[str, Any], temperature_K: float,
                         *, ln_offset: float = 0.0,
                         activation_energy_J_mol: float | None = None) -> float:
    reference_K = load_hcn_hydrolysis_rate_card()["reference_temperature_K"]
    energy = (
        item["activation_energy_J_mol_from_table_1"]
        if activation_energy_J_mol is None else activation_energy_J_mol
    )
    return math.exp(
        item["ln_k_at_reference_M_inv_s_inv"] + ln_offset
        - energy / R * (1.0 / temperature_K - 1.0 / reference_K)
    )


def hcn_hydrolysis_rate_screen(
    temperature_K: float, pH: float, pKa_hcn: float, pKw: float,
) -> dict[str, float]:
    """Return total-cyanide first-order hazard under fast acid-base equilibrium.

    The rate data apply to molecular HCN, so multiply by the equilibrium HCN
    fraction when expressing a hazard against total HCN + CN-.
    """

    if not all(math.isfinite(value) for value in (temperature_K, pH, pKa_hcn, pKw)):
        raise ValueError("HCN hydrolysis inputs must be finite")
    card = load_hcn_hydrolysis_rate_card()
    lower, upper = card["screen_temperature_K"]
    if not lower <= temperature_K <= upper:
        raise ValueError("HCN hydrolysis temperature is outside the screening window")
    acid, base = card["acid"], card["base"]
    hydrogen_M = 10.0 ** -pH
    hydroxide_M = 10.0 ** (pH - pKw)
    molecular_fraction = 1.0 / (1.0 + 10.0 ** (pH - pKa_hcn))

    def effective_hazard(acid_offset: float, base_offset: float,
                         acid_energy: float | None = None) -> float:
        return molecular_fraction * (
            _rate_at_temperature(
                acid, temperature_K, ln_offset=acid_offset,
                activation_energy_J_mol=acid_energy,
            ) * hydrogen_M
            + _rate_at_temperature(base, temperature_K, ln_offset=base_offset)
            * hydroxide_M
        )

    return {
        "nominal_s_inv": effective_hazard(0.0, 0.0),
        "parameter_low_s_inv": effective_hazard(
            -acid["sigma_ln_k_conservative"], -base["sigma_ln_k_conservative"]
        ),
        "parameter_high_s_inv": effective_hazard(
            acid["sigma_ln_k_conservative"], base["sigma_ln_k_conservative"]
        ),
        "abstract_acid_energy_sensitivity_s_inv": effective_hazard(
            0.0, 0.0, card["abstract_acid_activation_energy_J_mol"]
        ),
        "molecular_hcn_fraction": molecular_fraction,
    }


def evaluate_aqueous_hcn_hydrolysis_screen(
    aqueous_exposure: dict[str, Any],
) -> dict[str, Any]:
    """Screen time-varying aqueous loss from frozen, step-endpoint snapshots."""

    card = load_hcn_hydrolysis_rate_card()
    totals = {
        place: {"nominal_mol": 0.0, "parameter_low_mol": 0.0,
                "parameter_high_mol": 0.0, "abstract_acid_energy_mol": 0.0}
        for place in _AQUEOUS
    }
    unavailable: list[dict[str, Any]] = []
    unavailable_count = 0

    def record_unavailable(step: dict[str, Any], place: str, reason: str) -> None:
        nonlocal unavailable_count
        unavailable_count += 1
        if len(unavailable) < _MAX_UNAVAILABLE_DETAILS:
            unavailable.append({
                "cycle": step["cycle"], "substep": step.get("substep"),
                "compartment": place, "reason": reason,
            })

    for step in aqueous_exposure["step_snapshots"]:
        duration = float(step["duration_s"])
        if not math.isfinite(duration) or duration <= 0:
            raise ValueError("Hydrolysis screen requires positive finite step durations")
        for place in _AQUEOUS:
            liquid = step["compartments"][place]
            total_cyanide = (
                float(liquid["molecular_hcn_mol"]) + float(liquid["cyanide_ion_mol"])
            )
            if not math.isfinite(total_cyanide) or total_cyanide < 0:
                raise ValueError("Hydrolysis screen requires finite nonnegative cyanide")
            if total_cyanide == 0:
                continue
            if liquid["pH_status"] != "ideal_limited_species_pH_estimated":
                record_unavailable(step, place, liquid["pH_status"])
                continue
            try:
                rate = hcn_hydrolysis_rate_screen(
                    float(liquid["temperature_K"]),
                    float(liquid["pH_ideal_aqueous_model"]),
                    float(liquid["pKa_HCN"]), float(liquid["pKw"]),
                )
            except ValueError as exc:
                record_unavailable(step, place, str(exc))
                continue
            for output_key, rate_key in (
                ("nominal_mol", "nominal_s_inv"),
                ("parameter_low_mol", "parameter_low_s_inv"),
                ("parameter_high_mol", "parameter_high_s_inv"),
                ("abstract_acid_energy_mol", "abstract_acid_energy_sensitivity_s_inv"),
            ):
                totals[place][output_key] += total_cyanide * (
                    -math.expm1(-rate[rate_key] * duration)
                )
    combined = {
        key: sum(place[key] for place in totals.values())
        for key in ("nominal_mol", "parameter_low_mol", "parameter_high_mol",
                    "abstract_acid_energy_mol")
    }
    return {
        "status": (
            "unavailable_for_one_or_more_cyanide_exposures" if unavailable_count
            else "read_only_step_endpoint_hydrolysis_screen"
        ),
        "dataset_id": card["dataset_id"],
        "source": card["source"],
        "by_compartment": totals,
        "combined_screen_mol": combined,
        "unavailable_exposure_count": unavailable_count,
        "unavailable_exposures": unavailable,
        "prediction_inventory_mutated": False,
        "kinetic_products_injected": False,
        "thermal_ledger_closed_by_this_screen": False,
        "interpretation": (
            "These are right-endpoint quadrature screening estimates, not actual "
            "hydrolysis yields or strict bounds. The parameter range varies only the "
            "published log-rate spreads, omitting pH, salt, source, transport and "
            "within-step timing uncertainty. If any cyanide-bearing water step is "
            "unavailable, combined values are incomplete and must not be used as "
            "a complete loss estimate."
        ),
        "source_discrepancy": card["source_discrepancy"],
    }
