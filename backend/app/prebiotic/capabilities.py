from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version
import hashlib
from pathlib import Path
from typing import Any

from .boltzmann import (
    DATASET_ID,
    HYBRID_DATASET_ID,
    REACTIVE_N2_DATASET_ID,
    COMMUNITY_CH4_DATASET_ID,
    IST_H2_COMMUNITY_CH4_DATASET_ID,
    IST_CO2_O2_DATASET_ID,
    collision_dataset_status,
)


def _cantera_plasma() -> dict[str, Any]:
    try:
        import cantera as ct

        installed = version("Cantera")
        rate_available = hasattr(ct, "ElectronCollisionPlasmaRate")
        example_path: Path | None = None
        pavan_path: Path | None = None
        kamarinopoulou_path = (
            Path(__file__).resolve().parents[2]
            / "data"
            / "reference"
            / "kamarinopoulou_2024"
            / "kamarinopoulou_2024_neutral.yaml"
        )
        for directory in ct.get_data_directories():
            candidate = Path(directory) / "example_data" / "oxygen-plasma-itikawa.yaml"
            if candidate.is_file():
                example_path = candidate
            pavan_candidate = (
                Path(directory) / "example_data" / "methane-plasma-pavan-2023.yaml"
            )
            if pavan_candidate.is_file():
                pavan_path = pavan_candidate
        example_loaded = False
        species = reactions = None
        error = None
        if example_path is not None:
            try:
                # The packaged reference declares ion transport data incompletely for
                # some Windows wheels. Transport is irrelevant to this kinetics/API
                # capability probe, so disable it explicitly.
                phase = ct.Solution(
                    str(example_path),
                    "isotropic-electron-energy-plasma",
                    transport_model=None,
                )
                example_loaded = True
                species, reactions = phase.n_species, phase.n_reactions
            except (ct.CanteraError, ValueError) as caught:
                error = f"{type(caught).__name__}: {caught}"
        pavan_loaded = False
        pavan_species = pavan_reactions = None
        pavan_two_temperature_reactions = pavan_electron_collision_reactions = None
        if pavan_path is not None:
            try:
                pavan = ct.Solution(str(pavan_path), "gri30-plasma")
                pavan_loaded = True
                pavan_species, pavan_reactions = pavan.n_species, pavan.n_reactions
                pavan_rate_types = [type(reaction.rate).__name__ for reaction in pavan.reactions()]
                pavan_two_temperature_reactions = pavan_rate_types.count("TwoTempPlasmaRate")
                pavan_electron_collision_reactions = pavan_rate_types.count(
                    "ElectronCollisionPlasmaRate"
                )
            except (ct.CanteraError, ValueError) as caught:
                error = error or f"{type(caught).__name__}: {caught}"
        kamarinopoulou_loaded = False
        kamarinopoulou_species = kamarinopoulou_reactions = None
        kamarinopoulou_sha256 = None
        if kamarinopoulou_path.is_file():
            try:
                literature_hcn = ct.Solution(str(kamarinopoulou_path))
                kamarinopoulou_loaded = True
                kamarinopoulou_species = literature_hcn.n_species
                kamarinopoulou_reactions = literature_hcn.n_reactions
                kamarinopoulou_sha256 = hashlib.sha256(
                    kamarinopoulou_path.read_bytes()
                ).hexdigest()
            except (ct.CanteraError, ValueError) as caught:
                error = error or f"{type(caught).__name__}: {caught}"
        return {
            "available": rate_available and example_loaded,
            "version": installed,
            "electron_collision_rate": rate_available,
            "reference_example_loaded": example_loaded,
            "reference_example_species": species,
            "reference_example_reactions": reactions,
            "pavan_methane_plasma_loaded": pavan_loaded,
            "pavan_methane_plasma_species": pavan_species,
            "pavan_methane_plasma_reactions": pavan_reactions,
            "pavan_correlated_shadow_connected": pavan_loaded,
            "pavan_charge_balanced_shadow_connected": pavan_loaded,
            "pavan_retained_two_temperature_reactions": pavan_two_temperature_reactions,
            "pavan_disabled_embedded_electron_collision_reactions": pavan_electron_collision_reactions,
            "kamarinopoulou_hcn_shadow_connected": kamarinopoulou_loaded,
            "kamarinopoulou_hcn_species": kamarinopoulou_species,
            "kamarinopoulou_hcn_neutral_reactions": kamarinopoulou_reactions,
            "kamarinopoulou_hcn_fixed_electron_reactions_excluded": 62,
            "kamarinopoulou_hcn_mechanism_sha256": kamarinopoulou_sha256,
            "detail": error or "Cantera plasma API and packaged reference mechanism loaded",
        }
    except (ImportError, PackageNotFoundError) as error:
        return {
            "available": False,
            "version": None,
            "electron_collision_rate": False,
            "reference_example_loaded": False,
            "detail": f"{type(error).__name__}: {error}",
        }


def prebiotic_capabilities() -> dict[str, Any]:
    collision_set = collision_dataset_status(REACTIVE_N2_DATASET_ID)
    collision_sets = [
        collision_dataset_status(DATASET_ID),
        collision_dataset_status(HYBRID_DATASET_ID),
        collision_set,
        collision_dataset_status(COMMUNITY_CH4_DATASET_ID),
        collision_dataset_status(IST_H2_COMMUNITY_CH4_DATASET_ID),
        collision_dataset_status(IST_CO2_O2_DATASET_ID),
    ]
    boltzmann_available = collision_set.get("available", False)
    return {
        "apparatus": {
            "available": True,
            "compartments": ["boiler_liquid", "gas_loop", "condensate_trap", "collected_condensate", "wall"],
            "features": [
                "closed recirculation cycles",
                "evidence-labelled evaporation, condensation and wall loss",
                "gas and liquid mass-action networks",
                "mid-run additions and sampling",
                "element and charge closure",
                "electrical energy ledger",
                "evidence-labelled active discharge volume and whole-loop residence scaling",
                "cycle-resolved gas-inventory coupling to a Boltzmann EEDF",
                "read-only near-298 K O(1D) competing-fate audit with unsupported-bath refusal",
                "conditional primary radical, ion and excited-state source terms",
                "energy-limited neutral gas chemistry with Cantera GRI-Mech",
                "correlated Pavan 2023 methane-plasma heavy-chemistry shadow with mapped N2 excited states",
                "charge-balanced Pavan ion/electron shadow with two-temperature recombination and external EEDF sources",
                "independent charged-Pavan state propagated through apparatus cycles with branch-local EEDF and gas/water transfer",
                "operator-split continuous boil/discharge/condensation with explicit substep and convergence status",
                "time-averaged energy-limited neutral source and separate isolated cold-trap collection scenario",
                "full Kamarinopoulou 2024 CH4/N2 neutral HCN mechanism shadow with fixed-condition electron rates excluded",
                "reviewed Henry-law gas/water partitioning",
                "opt-in explicit aqueous HCN/CN-, NH3/NH4+ and carbonate equilibrium with atom/charge closure",
                "product-agnostic neutral CHNO aqueous candidate enumeration",
                "cycle-resolved post-discovery Strecker sensitivity replay",
                "held-out 1953 material/time benchmark",
                "background four-model formula-local candidate ranking",
                "frozen class-conditional held-out rank evaluation",
                "measured precursor energy-yield mismatch gate",
                "versioned 1-hour-to-7-day Stribling-Miller measurement-timing profile",
            ],
        },
        "reduced_order_discharge": {
            "available": True,
            "quantitative_yield_prediction": False,
            "boundary": "Declared pseudo-rate channels only; no excited-state or EEDF solution.",
        },
        "cantera_cross_section_plasma": {
            **_cantera_plasma(),
            "connected_to_apparatus": True,
            "rate_coefficient_calculation": True,
            "eedf_modes": ["maxwellian", "tabulated", "Boltzmann-two-term"],
            "self_consistent_boltzmann_eedf": boltzmann_available,
            "required_inputs": [
                "electron-impact cross sections",
                "electron-energy distribution",
                "electron number density or explicit electron inventory",
                "species thermochemistry",
                "gas-phase kinetic mechanism",
            ],
        },
        "boltzmann_eedf": {
            "available": boltzmann_available,
            "engine": "BOLOS",
            "version": "0.2" if boltzmann_available else None,
            "collision_set": collision_set,
            "collision_sets": collision_sets,
            "targets": ["H2", "H2O", "CH4", "NH3", "N2"],
            "opt_in_co2_o2_targets": ["CO2", "O2"],
            "mixture_composition_coupled": True,
            "reduced_electric_field_coupled": True,
            "apparatus_cycle_coupled": True,
            "discharge_zone_residence_coupled": True,
            "primary_source_terms": True,
            "heavy_species_inventory_coupled": True,
            "heavy_species_engine": "Cantera gri30.yaml",
            "heavy_species_shadow_engine": "Cantera methane-plasma-pavan-2023.yaml",
            "charge_balanced_plasma_shadow": True,
            "propagated_charged_branch": True,
            "propagated_charged_branch_self_consistent_eedf": True,
            "propagated_charged_branch_gas_water_partition": True,
            "continuous_overlap_operator_splitting": True,
            "automatic_timestep_convergence_claim": False,
            "literature_hcn_shadow_engine": "Cantera Kamarinopoulou 2024 Data S1/S2 import",
            "literature_hcn_shadow_reactions": 189,
            "heavy_species_shadow_is_independent_confidence": False,
            "reviewed_gas_water_partition": True,
            "opt_in_aqueous_hcn_ammonia_ion_equilibrium": True,
            "fully_coupled_pH_transport_and_kinetics": False,
            "aqueous_candidate_enumerator_coupled": True,
            "post_discovery_cycle_replay_coupled": True,
            "held_out_historical_benchmark_coupled": True,
            "background_candidate_energy_ranking": True,
            "frozen_candidate_rank_evaluation": True,
            "measured_precursor_energy_yield_gate": True,
            "achiral_enantiomer_degeneracy_enforced": True,
            "aqueous_reaction_network_coupled": False,
            "quantitative_yield_prediction": False,
        },
        "blind_validation": {
            "known_products_available_to_discovery": False,
            "evaluation_data_connected": True,
            "evaluation_data_available_only_after_discovery_freeze": True,
            "historical_yields_used_for_fitting": False,
        },
    }
