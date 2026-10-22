"""RDKit-backed identity metadata for intact, input-derived molecules.

The search engine stores ASE coordinates, which do not contain bond orders or
stereochemical tags.  This module keeps an atom-mapped RDKit template for each
input molecule and reassigns stereochemistry from the current 3D coordinates.
The template is used only while the same mapped atoms retain the same
connectivity; it is never used to name a newly formed product.
"""
from __future__ import annotations

from typing import Any

from ase import Atoms
from rdkit import Chem


def _potential_stereo(molecule: Chem.Mol) -> list[Any]:
    probe = Chem.RemoveHs(Chem.Mol(molecule))
    Chem.AssignStereochemistry(probe, cleanIt=True, force=True)
    return list(Chem.FindPotentialStereo(probe))


def _input_stereo_status(molecule: Chem.Mol) -> tuple[str, int]:
    potential = _potential_stereo(molecule)
    if not potential:
        return "achiral", 0
    unspecified = any(str(item.specified).endswith("Unspecified") for item in potential)
    return ("unassigned" if unspecified else "specified"), len(potential)


def make_instance_metadata(
    molecule: Chem.Mol,
    global_indices: tuple[int, ...],
    *,
    name: str,
    charge: int,
    connectivity_hash: str,
) -> dict[str, Any]:
    """Serialize the chemically authored identity of one input fragment."""
    mapped = Chem.Mol(molecule)
    for atom, index in zip(mapped.GetAtoms(), global_indices, strict=True):
        atom.SetAtomMapNum(index + 1)
    parent = Chem.RemoveHs(Chem.Mol(molecule))
    stereo_status, potential_count = _input_stereo_status(parent)
    return {
        "name": name,
        "charge": int(charge),
        "atom_indices": list(global_indices),
        "connectivity_hash": connectivity_hash,
        "mapped_template_smiles": Chem.MolToSmiles(
            mapped, canonical=False, isomericSmiles=True, allHsExplicit=True
        ),
        "canonical_smiles": Chem.MolToSmiles(parent, canonical=True, isomericSmiles=False),
        "input_isomeric_smiles": Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True),
        "input_stereo_status": stereo_status,
        "potential_stereo_count": potential_count,
        "aromatic": any(atom.GetIsAromatic() for atom in parent.GetAtoms()),
        "aromatic_atom_count": sum(atom.GetIsAromatic() for atom in parent.GetAtoms()),
        "aromatic_bond_count": sum(bond.GetIsAromatic() for bond in parent.GetBonds()),
    }


def _stereo_labels(molecule: Chem.Mol) -> list[str]:
    labels: list[str] = []
    centers = Chem.FindMolChiralCenters(
        molecule, force=True, includeUnassigned=True, useLegacyImplementation=False
    )
    for atom_index, descriptor in centers:
        labels.append(f"atom {atom_index + 1}: {descriptor}")
    for bond in molecule.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        stereo = str(bond.GetStereo())
        descriptor = None
        if stereo in {"STEREOE", "STEREOTRANS"}:
            descriptor = "E"
        elif stereo in {"STEREOZ", "STEREOCIS"}:
            descriptor = "Z"
        if descriptor:
            labels.append(
                f"bond {bond.GetBeginAtomIdx() + 1}-{bond.GetEndAtomIdx() + 1}: {descriptor}"
            )
    return labels


def identity_from_geometry(atoms: Atoms, metadata: dict[str, Any]) -> dict[str, Any]:
    """Reassign R/S and E/Z from current coordinates using a mapped template."""
    parser = Chem.SmilesParserParams()
    parser.removeHs = False
    template = Chem.MolFromSmiles(metadata["mapped_template_smiles"], parser)
    if template is None or template.GetNumAtoms() != len(metadata["atom_indices"]):
        return {
            "stereo_status": "indeterminate",
            "stereochemistry": None,
            "isomeric_smiles": None,
            "identity_key": f"{metadata['canonical_smiles']}|stereo:indeterminate",
        }

    Chem.RemoveStereochemistry(template)
    conformer = Chem.Conformer(template.GetNumAtoms())
    mapped_indices: list[int] = []
    for atom in template.GetAtoms():
        global_index = atom.GetAtomMapNum() - 1
        if global_index < 0 or global_index >= len(atoms):
            return {
                "stereo_status": "indeterminate",
                "stereochemistry": None,
                "isomeric_smiles": None,
                "identity_key": f"{metadata['canonical_smiles']}|stereo:indeterminate",
            }
        mapped_indices.append(global_index)
        x, y, z = atoms.positions[global_index]
        conformer.SetAtomPosition(atom.GetIdx(), (float(x), float(y), float(z)))
    template.RemoveAllConformers()
    template.AddConformer(conformer, assignId=True)
    Chem.AssignStereochemistryFrom3D(template, confId=0, replaceExistingTags=True)
    Chem.AssignStereochemistry(template, cleanIt=True, force=True)

    authored_status = metadata.get("input_stereo_status", "unknown")
    if authored_status == "unassigned":
        return {
            "stereo_status": "unassigned",
            "stereochemistry": "unassigned",
            "isomeric_smiles": metadata["canonical_smiles"],
            "identity_key": f"{metadata['canonical_smiles']}|stereo:unassigned",
        }

    labels = _stereo_labels(template)
    display_labels = [label.rsplit(": ", 1)[-1] for label in labels]
    for atom in template.GetAtoms():
        atom.SetAtomMapNum(0)
    parent = Chem.RemoveHs(template)
    Chem.AssignStereochemistry(parent, cleanIt=True, force=True)
    isomeric_smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)

    if authored_status == "specified" and (
        not labels or any(label.endswith(": ?") for label in labels)
    ):
        return {
            "stereo_status": "indeterminate",
            "stereochemistry": None,
            "isomeric_smiles": isomeric_smiles,
            "identity_key": f"{metadata['canonical_smiles']}|stereo:indeterminate",
        }
    if authored_status == "achiral":
        return {
            "stereo_status": "achiral",
            "stereochemistry": None,
            "isomeric_smiles": isomeric_smiles,
            "identity_key": f"{metadata['canonical_smiles']}|stereo:achiral",
        }
    return {
        "stereo_status": "specified",
        "stereochemistry": ",".join(display_labels),
        "stereo_assignments": labels,
        "isomeric_smiles": isomeric_smiles,
        "identity_key": f"{isomeric_smiles}|stereo:{','.join(display_labels)}",
    }
