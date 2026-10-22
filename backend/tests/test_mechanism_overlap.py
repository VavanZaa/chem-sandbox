"""The plasma shadow must not be misrepresented as independent HCN kinetics."""

from app.prebiotic.mechanism_overlap import audit_pavan_gri_overlap


def test_packaged_pavan_reuses_gri_hcn_reaction_records() -> None:
    audit = audit_pavan_gri_overlap()
    assert audit["gri_reaction_count"] == 325
    assert audit["pavan_reaction_count"] == 440
    assert audit["exact_shared_gri_reaction_records"] == 325
    assert audit["gri_reaction_records_all_reused"] is True
    assert audit["gri_hcn_reaction_count"] == 30
    assert audit["pavan_hcn_reaction_count"] == 30
    assert audit["exact_shared_hcn_reaction_records"] == 30
    assert audit["pavan_hcn_reactions_all_reused_from_gri"] is True
    assert audit["independent_hcn_kinetic_validation"] is False
    assert audit["pavan_hcnn_species_present"] is True
    assert audit["pavan_ncn_species_present"] is False
    assert audit["pavan_ch_n2_hcnn_pathway_present"] is True
    assert len(audit["gri_file_sha256"]) == 64
    assert len(audit["pavan_file_sha256"]) == 64
