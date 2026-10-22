"""Evaluator-only balanced gas-phase benchmark for installed ML potentials.

The reference values are standard gas-phase reaction enthalpies at 298.15 K
derived from NIST Chemistry WebBook formation enthalpies. Model values are
minimum electronic potential-energy differences without ZPE or thermal terms,
so the raw deviation is diagnostic rather than a publishable accuracy score.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256

import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


@dataclass(frozen=True, slots=True)
class Species:
    id: str
    name: str
    smiles: str


@dataclass(frozen=True, slots=True)
class ReactionCase:
    id: str
    name: str
    category: str
    difficulty: str
    stoichiometry: dict[str, int]  # products positive, reactants negative
    reference_kj_mol: float | None
    reference_note: str
    references: tuple[str, ...]
    eligible_engines: tuple[str, ...]
    catalyst_note: str | None = None


ALL_ENGINES = ("mace_mp", "ani2x", "aimnet2", "aimnet2_rxn")

SPECIES = {
    item.id: item
    for item in (
        Species("h2", "Hydrogen", "[H][H]"),
        Species("f2", "Fluorine", "FF"),
        Species("hf", "Hydrogen fluoride", "[H]F"),
        Species("acetylene", "Acetylene", "C#C"),
        Species("ethylene", "Ethylene", "C=C"),
        Species("ethane", "Ethane", "CC"),
        Species("nitrogen", "Nitrogen", "N#N"),
        Species("ammonia", "Ammonia", "N"),
        Species("carbon_monoxide", "Carbon monoxide", "[C-]#[O+]"),
        Species("water", "Water", "O"),
        Species("carbon_dioxide", "Carbon dioxide", "O=C=O"),
        Species("acetic_acid", "Acetic acid", "CC(=O)O"),
        Species("ethanol", "Ethanol", "CCO"),
        Species("ethyl_acetate", "Ethyl acetate", "CC(=O)OCC"),
        Species("butadiene", "1,3-Butadiene", "C=CC=C"),
        Species("cyclohexene", "Cyclohexene", "C1=CCCCC1"),
        Species("benzene", "Benzene", "c1ccccc1"),
        Species("cyclohexane", "Cyclohexane", "C1CCCCC1"),
        Species("pyridine", "Pyridine", "n1ccccc1"),
        Species("piperidine", "Piperidine", "N1CCCCC1"),
        Species("glucose", "beta-D-Glucopyranose", "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"),
        Species("sorbitol", "D-Sorbitol", "OCC(O)C(O)C(O)C(O)CO"),
    )
}

NIST = "https://webbook.nist.gov/cgi/cbook.cgi?ID="

REACTIONS = (
    ReactionCase(
        "hydrogen_fluorination", "Hydrogen fluorination", "high-energy", "small / strongly exothermic",
        {"h2": -1, "f2": -1, "hf": 2}, -546.60,
        "NIST gas-phase ΔfH°298; H2 and F2 elemental references are zero.",
        (NIST + "C7664393&Mask=1",), ("mace_mp", "ani2x", "aimnet2"),
    ),
    ReactionCase(
        "acetylene_hydrogenation", "Acetylene hydrogenation", "organic", "multiple-bond reduction",
        {"acetylene": -1, "h2": -2, "ethane": 1}, -311.40,
        "Derived from NIST gas-phase ΔfH°298 for acetylene and ethane.",
        (NIST + "C74862&Mask=1", NIST + "C74840&Mask=1"), ALL_ENGINES,
        "Industrially catalyzed; the catalyst is omitted from this thermodynamic endpoint test.",
    ),
    ReactionCase(
        "ethylene_hydrogenation", "Ethylene hydrogenation", "organic", "single π-bond reduction",
        {"ethylene": -1, "h2": -1, "ethane": 1}, -137.15,
        "Derived from NIST gas-phase ΔfH°298 for ethylene and ethane.",
        (NIST + "C74851&Mask=1", NIST + "C74840&Mask=1"), ALL_ENGINES,
        "Normally metal-catalyzed; the catalyst affects the barrier, not the balanced endpoint energy.",
    ),
    ReactionCase(
        "haber_bosch", "Ammonia synthesis", "industrial", "strong triple bond / stoichiometric",
        {"nitrogen": -1, "h2": -3, "ammonia": 2}, -91.80,
        "Derived from NIST gas-phase ΔfH°298 for ammonia; elemental gases are zero.",
        (NIST + "C7664417&Mask=1",), ALL_ENGINES,
        "Haber-Bosch uses an iron-based catalyst. It is omitted because ANI/AIMNet molecular models do not cover the metal surface.",
    ),
    ReactionCase(
        "water_gas_shift", "Water-gas shift", "industrial", "small-molecule rearrangement",
        {"carbon_monoxide": -1, "water": -1, "carbon_dioxide": 1, "h2": 1}, -41.15,
        "Derived from NIST gas-phase ΔfH°298 for CO, H2O and CO2.",
        (NIST + "C630080&Mask=1", NIST + "C7732185&Mask=1", NIST + "C124389&Mask=1"), ALL_ENGINES,
        "Industrial catalysts include Fe/Cr and Cu/Zn systems; no catalyst surface is represented here.",
    ),
    ReactionCase(
        "fischer_esterification", "Fischer esterification", "organic / catalytic", "four oxygenated molecules",
        {"acetic_acid": -1, "ethanol": -1, "ethyl_acetate": 1, "water": 1}, -20.23,
        "Derived from NIST gas-phase ΔfH°298 averages for acetic acid, ethanol, ethyl acetate and water.",
        (NIST + "C64197&Mask=1", NIST + "C64175&Mask=1", NIST + "C141786&Mask=1", NIST + "C7732185&Mask=1"),
        ALL_ENGINES,
        "Acid catalysis changes the pathway and rate. H+ is not included in this endpoint-only benchmark.",
    ),
    ReactionCase(
        "diels_alder", "Butadiene + ethylene Diels-Alder", "organic", "concerted bond rearrangement",
        {"butadiene": -1, "ethylene": -1, "cyclohexene": 1}, -165.59,
        "Derived from NIST gas-phase ΔfH°298: butadiene 108.8, ethylene 52.47, cyclohexene -4.32 kJ/mol.",
        (NIST + "C106990&Mask=1", NIST + "C74851&Mask=1", NIST + "C110838&Mask=1"), ALL_ENGINES,
    ),
    ReactionCase(
        "benzene_hydrogenation", "Benzene hydrogenation", "aromatic / catalytic", "aromatic ring reduction",
        {"benzene": -1, "h2": -3, "cyclohexane": 1}, -206.0,
        "Derived from NIST gas-phase ΔfH°298: benzene 82.9 and cyclohexane -123.1 kJ/mol; H2 is zero.",
        (NIST + "C71432&Mask=181", NIST + "C110827&Mask=1F"), ALL_ENGINES,
        "Industrial hydrogenation uses a metal catalyst. This endpoint test does not model its surface, barrier or selectivity.",
    ),
    ReactionCase(
        "pyridine_hydrogenation", "Pyridine hydrogenation", "heteroaromatic / catalytic", "heteroaromatic ring reduction",
        {"pyridine": -1, "h2": -3, "piperidine": 1}, -193.8,
        "Direct NIST gas-phase reaction enthalpy for pyridine + 3 H2 → piperidine.",
        (NIST + "C110894&Mask=80F&Units=SI",), ALL_ENGINES,
        "Catalyst and adsorption state are omitted; this validates only isolated-molecule endpoint thermochemistry.",
    ),
    ReactionCase(
        "glucose_hydrogenation", "Glucose hydrogenation to sorbitol", "complex / industrial", "24-atom flexible product",
        {"glucose": -1, "h2": -1, "sorbitol": 1}, None,
        "Stress test only: no phase-matched NIST gas-phase ΔH°298 reference is assigned.",
        ("https://pubchem.ncbi.nlm.nih.gov/compound/Glucose", "https://pubchem.ncbi.nlm.nih.gov/compound/Sorbitol"),
        ALL_ENGINES,
        "Commercial hydrogenation uses nickel or ruthenium catalysts; this case tests molecular endpoint robustness, not the catalytic barrier.",
    ),
)


def rdkit_molecule(species: Species) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(species.smiles)
    if molecule is None:
        raise ValueError(f"Invalid benchmark SMILES for {species.id}")
    return Chem.AddHs(molecule)


def atom_inventory(species: Species) -> Counter[str]:
    return Counter(atom.GetSymbol() for atom in rdkit_molecule(species).GetAtoms())


def validate_suite() -> None:
    """Fail if a case is unbalanced, unknown, charged overall or open-shell."""
    for case in REACTIONS:
        balance: Counter[str] = Counter()
        charge = radicals = 0
        for species_id, coefficient in case.stoichiometry.items():
            species = SPECIES[species_id]
            molecule = rdkit_molecule(species)
            balance.update({element: count * coefficient for element, count in atom_inventory(species).items()})
            charge += Chem.GetFormalCharge(molecule) * coefficient
            radicals += sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()) * abs(coefficient)
        nonzero_balance = {element: count for element, count in balance.items() if count}
        if nonzero_balance or charge or radicals:
            raise ValueError(f"Unbalanced or unsupported benchmark case {case.id}: atoms={nonzero_balance}, charge={charge}, radicals={radicals}")


def build_geometry(species: Species, seed: int = 20260922) -> tuple[Atoms, set[tuple[int, int]], str]:
    """Create one deterministic, MMFF/UFF-prepared isolated molecule."""
    molecule = rdkit_molecule(species)
    params = AllChem.ETKDGv3()
    digest = int.from_bytes(sha256(species.id.encode()).digest()[:4], "little")
    params.randomSeed = int((seed + digest) % (2**31 - 1))
    if AllChem.EmbedMolecule(molecule, params) < 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(molecule, params) < 0:
            raise ValueError(f"RDKit could not embed {species.name}")
    if AllChem.MMFFHasAllMoleculeParams(molecule):
        AllChem.MMFFOptimizeMolecule(molecule, maxIters=500)
    elif AllChem.UFFHasAllMoleculeParams(molecule):
        AllChem.UFFOptimizeMolecule(molecule, maxIters=500)
    positions = np.asarray(molecule.GetConformer().GetPositions(), dtype=float)
    atoms = Atoms([atom.GetSymbol() for atom in molecule.GetAtoms()], positions=positions)
    atoms.center(vacuum=6.0)
    atoms.info.update(
        chemistry_mode="general",
        total_charge=int(Chem.GetFormalCharge(molecule)),
        has_ions=False,
        radical_electrons=sum(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()),
        benchmark_species=species.id,
    )
    expected = {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))) for bond in molecule.GetBonds()}
    return atoms, expected, rdMolDescriptors.CalcMolFormula(molecule)


def inferred_connectivity(atoms: Atoms, scale: float = 1.30) -> set[tuple[int, int]]:
    symbols = atoms.get_atomic_numbers()
    distances = atoms.get_all_distances(mic=False)
    bonds: set[tuple[int, int]] = set()
    for left in range(len(atoms)):
        for right in range(left + 1, len(atoms)):
            cutoff = scale * (covalent_radii[symbols[left]] + covalent_radii[symbols[right]])
            if distances[left, right] <= cutoff:
                bonds.add((left, right))
    return bonds


validate_suite()
