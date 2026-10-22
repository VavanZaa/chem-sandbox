from scripts.audit_vulcan_lowt import composition, parse_reaction


def test_vulcan_formula_parser_preserves_atoms_and_excited_suffix() -> None:
    assert composition("CH2_1") == {"C": 1, "H": 2}
    assert composition("CH3NH2") == {"C": 1, "H": 5, "N": 1}
    assert composition("M") == {}


def test_vulcan_audit_detects_balance_and_declared_range() -> None:
    balanced = parse_reaction("1 [ CH4 + H -> CH3 + H2 ] 2.20E-20 3.0 4040.0 300-2500")
    assert balanced is not None
    assert balanced["atom_balanced"] is True
    assert balanced["explicit_temperature_ranges_K"] == [(300, 2500)]
    assert balanced["thermal_rate_triplet_present"] is True
    unbalanced = parse_reaction("2 [ CH4 -> CH3 + H2 ] 1.0E-10 0.0 0.0 ??")
    assert unbalanced is not None
    assert unbalanced["atom_balanced"] is False
    assert unbalanced["unknown_or_estimated_provenance"] is True
