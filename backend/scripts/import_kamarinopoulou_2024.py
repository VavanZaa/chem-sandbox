"""Reproducibly import the open Kamarinopoulou et al. 2024 supplement.

The Science Advances supplement embeds Data S1 (CHEMKIN thermochemistry) and
Data S2 (reaction input) as text on PDF pages 22-35.  This importer preserves
the extracted published input and creates a neutral/heavy-chemistry variant
for Cantera.  The paper's fixed-condition electron-impact rates are preserved
but deliberately excluded from that variant because Chem Sandbox supplies
electron source terms from its own mixture EEDF.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

from pypdf import PdfReader


EXPECTED_PDF_SHA256 = "43520da71dc844752958cb586cd251e9b7f8e49b707ff53a1c91e33f3dddc084"
THERMO_PAGE_INDICES = range(21, 28)
REACTION_PAGE_INDICES = range(28, 35)
ELECTRON_MARKER = "! Electron impact reactions from Bolsig+"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _page_text(reader: PdfReader, indices: range) -> str:
    return "\n".join((reader.pages[index].extract_text() or "") for index in indices)


def _extract_blocks(pdf_path: Path) -> tuple[str, str]:
    actual = _sha256(pdf_path)
    if actual != EXPECTED_PDF_SHA256:
        raise ValueError(
            f"Supplement SHA-256 mismatch: expected {EXPECTED_PDF_SHA256}, got {actual}"
        )
    reader = PdfReader(pdf_path)
    if len(reader.pages) < 35:
        raise ValueError(f"Expected at least 35 pages; supplement has {len(reader.pages)}")

    thermo_pages = _page_text(reader, THERMO_PAGE_INDICES)
    thermo_start = thermo_pages.index("THERMO")
    thermo_end = thermo_pages.index("END", thermo_start) + len("END")
    thermo = thermo_pages[thermo_start:thermo_end].strip() + "\n"

    reaction_pages = _page_text(reader, REACTION_PAGE_INDICES)
    reaction_start = reaction_pages.index("REACTIONS")
    reaction_end = reaction_pages.index("END", reaction_start) + len("END")
    reactions = reaction_pages[reaction_start:reaction_end].strip() + "\n"
    return thermo, reactions


def _thermo_records(thermo: str) -> tuple[list[str], list[tuple[str, ...]]]:
    lines = thermo.splitlines()
    records = [
        line
        for line in lines[1:]
        if line.strip() and not line.lstrip().startswith("!") and line.strip() != "END"
    ]
    if len(records) % 4:
        raise ValueError(
            f"Thermo extraction has {len(records)} non-comment lines, not a multiple of four"
        )
    blocks = [tuple(records[index : index + 4]) for index in range(0, len(records), 4)]
    unique: dict[str, tuple[str, ...]] = {}
    for block in blocks:
        name = block[0].split()[0]
        normalized = tuple(" ".join(line.split()) for line in block)
        existing = tuple(" ".join(line.split()) for line in unique.get(name, ()))
        if name in unique and existing != normalized:
            raise ValueError(f"Thermo extraction contains conflicting records for {name}")
        unique.setdefault(name, block)
    return list(unique), list(unique.values())


def _deduplicated_thermo(thermo: str, blocks: list[tuple[str, ...]]) -> str:
    comments = [
        line
        for line in thermo.splitlines()[1:]
        if line.lstrip().startswith("!")
    ]
    lines = ["THERMO ALL", "   100.000  1000.000  5000.000", *comments]
    for block in blocks:
        header_numbers = re.findall(r"[-+]?\d+(?:\.\d+)?", block[0])
        if len(header_numbers) < 4:
            raise ValueError(f"Could not read thermo temperature range: {block[0]}")
        temperature_low, temperature_high, temperature_mid = map(
            float, header_numbers[-4:-1]
        )
        name = block[0].split()[0]
        before_phase = block[0].rsplit(" G ", maxsplit=1)[0]
        composition = {
            element: int(count)
            for element, count in re.findall(r"([CHN])\s*(\d+)", before_phase)
        }
        if not composition:
            raise ValueError(f"Could not read thermo composition: {block[0]}")
        composition_field = "".join(
            f"{element:<2}{count:>3}" for element, count in composition.items()
        ).ljust(20)
        canonical_header = (
            f"{name:<24}"[:24]
            + composition_field
            + "G"
            + f"{temperature_low:10.3f}{temperature_high:10.3f}{temperature_mid:10.3f}"
        ).ljust(79) + "1"
        coefficients = [
            float(value)
            for line in block[1:]
            for value in re.findall(r"[-+]?\d+\.\d+E[-+]\d+", line)
        ]
        if len(coefficients) != 14:
            raise ValueError(
                f"Expected 14 NASA7 coefficients for {name}; found {len(coefficients)}"
            )
        coefficient_rows = (
            coefficients[:5],
            coefficients[5:10],
            coefficients[10:],
        )
        lines.append(canonical_header)
        for marker, values in zip("234", coefficient_rows, strict=True):
            lines.append("".join(f"{value:15.8E}" for value in values).ljust(79) + marker)
    lines.append("END")
    return "\n".join(lines) + "\n"


def _normalize_comments(text: str) -> str:
    equation_keys = [
        line.rsplit(maxsplit=3)[0]
        for line in text.splitlines()
        if "<=>" in line or "=>" in line
    ]
    equation_counts = {key: equation_keys.count(key) for key in set(equation_keys)}
    normalized: list[str] = []
    previous_equation: str | None = None
    for raw_line in text.splitlines():
        line = raw_line.rstrip()
        if "<=>" in line or "=>" in line:
            previous_equation = line.rsplit(maxsplit=3)[0]
        if line == "DUPLICATE" and equation_counts.get(previous_equation or "", 0) < 2:
            normalized.append(
                "! Orphan DUPLICATE marker omitted: the published subset contains no matching equation"
            )
            continue
        recognized = (
            not line
            or line.startswith("!")
            or line in {"REACTIONS", "END", "DUPLICATE"}
            or "<=>" in line
            or "=>" in line
            or line.startswith(("LOW/", "TROE/"))
            or bool(re.fullmatch(r"(?:[A-Za-z0-9()'+-]+/[.0-9+-]+/\s*)+", line))
        )
        if not recognized:
            # PDF line wrapping splits URLs and long bibliographic comments.
            line = "! " + line
        if line.startswith("LOW/"):
            values = line.removeprefix("LOW/").removesuffix("/").split()
            values = [value for value in values if value != "LOW"]
            if len(values) == 2:
                values.append("0")
            line = "LOW/ " + " ".join(values) + "/"
        normalized.append(line)
    return "\n".join(normalized) + "\n"


def _mechanism_input(species: list[str], reaction_block: str) -> str:
    header = [
        "! Extracted from Kamarinopoulou et al., Science Advances 10 (2024), adl4246",
        "ELEMENTS",
        "C H N",
        "END",
        "SPECIES",
    ]
    for offset in range(0, len(species), 8):
        header.append(" ".join(species[offset : offset + 8]))
    header.append("END")
    return "\n".join(header) + "\n" + _normalize_comments(reaction_block)


def _count_equation_records(text: str) -> int:
    return sum("<=>" in line or "=>" in line for line in text.splitlines())


def import_supplement(pdf_path: Path, output_dir: Path) -> dict[str, object]:
    thermo, reactions = _extract_blocks(pdf_path)
    species, thermo_blocks = _thermo_records(thermo)
    before_electron, electron_tail = reactions.split(ELECTRON_MARKER, maxsplit=1)
    neutral_reactions = before_electron.rstrip() + "\nEND\n"
    electron_reactions = "REACTIONS\n" + electron_tail.strip() + "\n"

    output_dir.mkdir(parents=True, exist_ok=True)
    published_thermo_path = output_dir / "kamarinopoulou_2024_published_thermo.dat"
    thermo_path = output_dir / "kamarinopoulou_2024_neutral_thermo.dat"
    full_path = output_dir / "kamarinopoulou_2024_published.inp"
    neutral_path = output_dir / "kamarinopoulou_2024_neutral.inp"
    electron_path = output_dir / "kamarinopoulou_2024_fixed_electron.inp"
    yaml_path = output_dir / "kamarinopoulou_2024_neutral.yaml"

    published_thermo_path.write_text(thermo, encoding="utf-8", newline="\n")
    thermo_path.write_text(
        _deduplicated_thermo(thermo, thermo_blocks), encoding="utf-8", newline="\n"
    )
    full_path.write_text(_mechanism_input(species, reactions), encoding="utf-8", newline="\n")
    neutral_path.write_text(
        _mechanism_input(species, neutral_reactions), encoding="utf-8", newline="\n"
    )
    electron_path.write_text(
        _mechanism_input(species, electron_reactions), encoding="utf-8", newline="\n"
    )

    command = [
        sys.executable,
        "-m",
        "cantera.ck2yaml",
        "--input",
        str(neutral_path),
        "--thermo",
        str(thermo_path),
        "--output",
        str(yaml_path),
        "--permissive",
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode:
        raise RuntimeError(
            "Cantera conversion failed:\n" + completed.stdout + "\n" + completed.stderr
        )
    # ck2yaml writes the wall-clock conversion time into the YAML. Remove that
    # non-scientific field so the generated mechanism has a stable content hash.
    yaml_text = yaml_path.read_text(encoding="utf-8")
    yaml_text = re.sub(r"^date:.*\n\n", "", yaml_text, flags=re.MULTILINE)
    yaml_path.write_text(yaml_text, encoding="utf-8", newline="\n")

    manifest = {
        "source": {
            "title": "Direct HCN synthesis via plasma-assisted conversion of methane and nitrogen",
            "authors": "Kamarinopoulou et al.",
            "journal": "Science Advances 10 (2024), adl4246",
            "doi": "10.1126/sciadv.adl4246",
            "article_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10980260/",
            "supplement_url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10980260/supplementaryFiles",
            "license": "CC BY-NC 4.0",
            "retrieved_utc_date": "2026-09-23",
            "data_s1_pages": "22-28",
            "data_s2_pages": "29-35",
        },
        "source_pdf_sha256": _sha256(pdf_path),
        "thermo_species": len(species),
        "published_reaction_records": _count_equation_records(reactions),
        "neutral_heavy_reaction_records": _count_equation_records(neutral_reactions),
        "fixed_condition_electron_reaction_records": _count_equation_records(
            electron_reactions
        ),
        "generated_files": {
            path.name: _sha256(path)
            for path in (
                published_thermo_path,
                thermo_path,
                full_path,
                neutral_path,
                electron_path,
                yaml_path,
            )
        },
        "exact_duplicate_thermo_records_removed_for_conversion": 1,
        "conversion_normalizations": [
            "Preserved the raw PDF-extracted Data S1 and Data S2 blocks in published files.",
            "Removed one textually identical duplicate CH3CN thermo record from the executable thermo file.",
            "Reconstructed NASA7 fixed-width columns from the printed coefficients without changing numeric values.",
            "Inserted explicit zero activation energies where printed LOW records omit the third CHEMKIN parameter.",
            "Removed orphan DUPLICATE markers only when the published subset contains no second matching equation.",
            "Converted PDF-wrapped URLs and bibliography continuations back to comments.",
        ],
        "scientific_firewall": {
            "paper_experimental_yields_imported_into_mechanism": False,
            "fixed_300_Td_1e11_cm-3_electron_rates_used": False,
            "neutral_heavy_reactions_used_as_shadow_only": True,
        },
    }
    (output_dir / "import_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdf", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    manifest = import_supplement(args.pdf.resolve(), args.output_dir.resolve())
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
