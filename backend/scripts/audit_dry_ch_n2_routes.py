"""Assay-free dry-flow check of thermal CH + N2 route competition.

The GRI direct-HCN and association rates are compared at a declared uniform
gas state. This is not a measured product branch for a hot streamer or a
prediction of HCN through the unrepresented NCN/HNCN network.
"""
from __future__ import annotations

import json
from typing import Any

import cantera as ct

from scripts.run_dry_flow_hcn_probe import DOWNSTREAM_PROBE_TEMPERATURE_K, run_dry_flow_probe


DIRECT_HCN_REACTION_INDEX = 239
DIRECT_HCN_REACTION = "CH + N2 <=> HCN + N"
ASSOCIATION_REACTION_INDEX = 240
ASSOCIATION_REACTION = "CH + N2 (+M) <=> HCNN (+M)"


def audit_dry_ch_n2_routes() -> dict[str, Any]:
    gas = ct.Solution("gri30.yaml")
    gas.TPX = DOWNSTREAM_PROBE_TEMPERATURE_K, ct.one_atm, {"CH4": 0.01, "N2": 0.99}
    if gas.reaction(DIRECT_HCN_REACTION_INDEX).equation != DIRECT_HCN_REACTION:
        raise ValueError("The installed GRI direct CH+N2 reaction index has changed")
    if gas.reaction(ASSOCIATION_REACTION_INDEX).equation != ASSOCIATION_REACTION:
        raise ValueError("The installed GRI CH+N2 association index has changed")
    direct = float(gas.forward_rate_constants[DIRECT_HCN_REACTION_INDEX])
    association = float(gas.forward_rate_constants[ASSOCIATION_REACTION_INDEX])
    if direct < 0 or association <= 0:
        raise ValueError("Expected finite positive CH+N2 association rate")

    probe = run_dry_flow_probe(source_forcing="evolving_eedf", substeps=16)
    hcn_pathways = probe["heavy_chemistry"]["hcn_pathways_at_outlet"]
    outlet_direct = next(
        (row for row in hcn_pathways["formation"]
         if row["reaction_index"] == DIRECT_HCN_REACTION_INDEX), None
    )
    return {
        "status": "thermal_ch_n2_route_applicability_not_hcn_validation",
        "uniform_gas_temperature_K_scenario": DOWNSTREAM_PROBE_TEMPERATURE_K,
        "pressure_Pa_scenario": ct.one_atm,
        "feed_mole_fractions": {"CH4": 0.01, "N2": 0.99},
        "gri_direct_hcn_reaction": DIRECT_HCN_REACTION,
        "gri_association_reaction": ASSOCIATION_REACTION,
        "gri_direct_rate_m3_kmol_s": direct,
        "gri_association_effective_rate_m3_kmol_s": association,
        "direct_over_association_forward_rate_ratio": direct / association,
        "outlet_direct_hcn_net_rate_mol_s": (
            outlet_direct["species_net_rate_mol_s"] if outlet_direct else 0.0
        ),
        "outlet_hcn_mol": probe["outlet_hcn_mol"],
        "outlet_hcn_below_1e_12_mol_reporting_floor": probe["outlet_hcn_mol"] < 1e-12,
        "n2_total_dissociation_source_excluded": (
            probe["heavy_chemistry"]["excluded_fragment_channels_by_id"].get(
                "N2_dissociation_outcome_000", {}
            ).get("reason") == "unresolved_electronic_product_state"
        ),
        "external_thermal_rate_evidence": {
            "source": "Brownsword et al., J. Chem. Soc. Faraday Trans. 92, 1087 (1996)",
            "doi": "10.1039/FT9969201087",
            "temperature_range_K": [202.0, 584.0],
            "bath_and_pressure_caveat": (
                "Argon-supported association measurements extend to about 0.53 bar; "
                "this audit uses 1 atm mostly N2 and a uniform temperature proxy."
            ),
            "reported_qualitative_result": (
                "Association dominates CH+N2 loss; direct HCN+N is too slow "
                "to matter over the study's conditions."
            ),
        },
        "direct_over_association_ratio_is_complete_product_branch": False,
        "hcnn_ncn_hncn_fate_resolved": False,
        "hcn_assay_loaded": False,
        "solver_or_exposure_modified": False,
        "interpretation": (
            "GRI's thermal direct HCN+N channel is dwarfed by its HCNN "
            "association entry at this assumed bulk gas state. The ratio "
            "is not a branching measurement or final HCN yield; hot streamer "
            "CH/N2 state distributions and HCNN/NCN/HNCN fates are unresolved. "
            "An N2 atom source is one possible missing route, not a logical "
            "prerequisite for every HCN-forming route."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit_dry_ch_n2_routes(), indent=2))
