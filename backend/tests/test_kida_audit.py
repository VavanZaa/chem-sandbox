"""The optional broad-network audit must preserve raw reaction integrity."""

from __future__ import annotations

from scripts.audit_kida_v23 import audit, reaction_slots


def fixed_width_line(reactants: list[str], products: list[str]) -> str:
    names = [*reactants, *([""] * (5 - len(reactants))),
             *products, *([""] * (5 - len(products)))]
    return "".join(f" {name:<10}" for name in names) + "  1.0e-10  0.0  0.0"


def test_reaction_slot_parser_preserves_repeated_species() -> None:
    line = fixed_width_line(["H", "H", "M"], ["H2", "M"])
    assert reaction_slots(line) == (["H", "H", "M"], ["H2", "M"])


def test_network_audit_flags_unbalanced_records_without_loading_solver(tmp_path) -> None:
    directory = tmp_path / "V23_KIDA"
    directory.mkdir()
    (directory / "composes.dat").write_text(
        "1 H 1.0 0 1 0 0\n2 H2 2.0 0 2 0 0\n", encoding="utf-8"
    )
    (directory / "readme.txt").write_text(
        "contains 2 chemical reactions (1 reversible and 1 irreversible)",
        encoding="utf-8",
    )
    (directory / "reactions_12.dat").write_text(
        fixed_width_line(["H", "H"], ["H2"]) + "\n", encoding="utf-8"
    )
    (directory / "reactions_10.dat").write_text(
        fixed_width_line(["H2"], ["H"]) + "\n", encoding="utf-8"
    )
    report = audit(directory)
    assert report["reaction_count"] == 2
    assert report["archive_matches_readme_counts"] is True
    assert report["unbalanced_reaction_count"] == 1
    assert report["unbalanced_reactions"][0]["left_CHON"] == (0, 2, 0, 0)
    assert report["not_imported_into_blind_solver"] is True
