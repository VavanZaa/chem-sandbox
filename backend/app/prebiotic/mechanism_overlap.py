"""Read-only provenance audit for the packaged GRI and Pavan mechanisms."""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Any

import cantera as ct

from .heavy_gas import MECHANISM, PAVAN_MECHANISM, PAVAN_PHASE


def _reaction_key(reaction: ct.Reaction) -> str:
    # A multiset, not an equation-keyed dictionary: GRI has repeated
    # equations with distinct rate records.
    return json.dumps(reaction.input_data, sort_keys=True)


@lru_cache(maxsize=1)
def audit_pavan_gri_overlap() -> dict[str, Any]:
    """Count exact shared reaction records without reading any HCN assay.

    Sharing a kinetic record is not independent mechanism validation even
    though the two coupled reactor states and source terms can differ.
    """

    gri_file = next(
        (
            Path(directory) / MECHANISM
            for directory in ct.get_data_directories()
            if (Path(directory) / MECHANISM).is_file()
        ),
        None,
    )
    packaged = next(
        (
            Path(directory) / PAVAN_MECHANISM
            for directory in ct.get_data_directories()
            if (Path(directory) / PAVAN_MECHANISM).is_file()
        ),
        None,
    )
    if gri_file is None or packaged is None:
        raise ValueError("Cantera's packaged GRI or Pavan mechanism is unavailable")
    gri = ct.Solution(str(gri_file))
    pavan = ct.Solution(str(packaged), PAVAN_PHASE)
    gri_records = Counter(_reaction_key(item) for item in gri.reactions())
    pavan_records = Counter(_reaction_key(item) for item in pavan.reactions())
    shared = gri_records & pavan_records
    gri_hcn = Counter(
        _reaction_key(item) for item in gri.reactions()
        if "HCN" in item.equation
    )
    pavan_hcn = Counter(
        _reaction_key(item) for item in pavan.reactions()
        if "HCN" in item.equation
    )
    shared_hcn = gri_hcn & pavan_hcn
    pavan_equations = set(pavan.reaction_equations())
    return {
        "cantera_version": ct.__version__,
        "gri_mechanism": MECHANISM,
        "gri_file_sha256": hashlib.sha256(gri_file.read_bytes()).hexdigest(),
        "pavan_mechanism": PAVAN_MECHANISM,
        "pavan_file_sha256": hashlib.sha256(packaged.read_bytes()).hexdigest(),
        "gri_reaction_count": gri.n_reactions,
        "pavan_reaction_count": pavan.n_reactions,
        "exact_shared_gri_reaction_records": sum(shared.values()),
        "gri_reaction_records_all_reused": sum(shared.values()) == gri.n_reactions,
        "gri_hcn_reaction_count": sum(gri_hcn.values()),
        "pavan_hcn_reaction_count": sum(pavan_hcn.values()),
        "exact_shared_hcn_reaction_records": sum(shared_hcn.values()),
        "pavan_hcn_reactions_all_reused_from_gri": (
            sum(shared_hcn.values()) == sum(pavan_hcn.values())
        ),
        "pavan_hcnn_species_present": "HCNN" in pavan.species_names,
        "pavan_ncn_species_present": "NCN" in pavan.species_names,
        "pavan_ch_n2_hcnn_pathway_present": (
            "CH + N2 (+M) <=> HCNN (+M)" in pavan_equations
        ),
        "independent_hcn_kinetic_validation": False,
        "scope": (
            "Exact kinetic-record ancestry only. Pavan adds plasma species and "
            "reactions, and its evolving gas state may differ; identical HCN "
            "rates do not make the trajectories identical or validate yield."
        ),
    }
