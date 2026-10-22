from evaluation.reference_only.reaction_suite import ALL_ENGINES, REACTIONS, SPECIES, atom_inventory, build_geometry, validate_suite


def test_reaction_benchmark_suite_is_balanced_and_has_required_categories():
    validate_suite()
    categories = {case.category for case in REACTIONS}
    assert any("industrial" in category for category in categories)
    assert any("organic" in category for category in categories)
    assert any("complex" in category for category in categories)
    assert any(case.catalyst_note for case in REACTIONS)
    assert {"glucose", "sorbitol"} <= set(SPECIES)


def test_reaction_references_and_engine_domains_are_explicit():
    assert all(case.references for case in REACTIONS)
    assert all(set(case.eligible_engines) <= set(ALL_ENGINES) for case in REACTIONS)
    fluorination = next(case for case in REACTIONS if case.id == "hydrogen_fluorination")
    assert "aimnet2_rxn" not in fluorination.eligible_engines
    assert fluorination.reference_kj_mol == -546.60


def test_glucose_hydrogenation_geometry_and_balance():
    glucose, bonds, formula = build_geometry(SPECIES["glucose"])
    assert formula == "C6H12O6"
    assert len(glucose) == 24
    assert bonds
    assert atom_inventory(SPECIES["sorbitol"])["H"] == 14
