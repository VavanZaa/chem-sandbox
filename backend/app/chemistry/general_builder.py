"""SMILES/ETKDG molecular assembly with explicit engine-domain checks."""
import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from .molecule_builder import BuiltSystem
from .species_detector import SpeciesDetector
from .stereochemistry import make_instance_metadata

ANI_ELEMENTS = {"H", "C", "N", "O", "F", "S", "Cl"}
AIMNET_ELEMENTS = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "As", "Se", "Br", "I"}
AIMNET_RXN_ELEMENTS = {"H", "C", "N", "O"}


def parse_molecule(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or not molecule.GetNumAtoms():
        raise ValueError(f"Invalid SMILES: {smiles}")
    if any(atom.GetAtomicNum() == 0 for atom in molecule.GetAtoms()):
        raise ValueError("Wildcard atoms are not supported")
    return Chem.AddHs(molecule)


def validate_molecules(inputs, engine="mace_mp"):
    atom_count = 0
    for item in inputs:
        molecule = parse_molecule(item.smiles)
        atom_count += molecule.GetNumAtoms() * item.copies
        if engine == "ani2x":
            if any(atom.GetSymbol() not in ANI_ELEMENTS or atom.GetFormalCharge() or atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()):
                raise ValueError("ANI-2x is restricted here to neutral, closed-shell HCNOFSCl molecules (no ions, radicals or metals)")
        if engine == "aimnet2":
            if any(atom.GetSymbol() not in AIMNET_ELEMENTS or atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()):
                raise ValueError("AIMNet2 is restricted here to closed-shell H/B/C/N/O/F/Si/P/S/Cl/As/Se/Br/I molecules; sodium and radicals require another model")
        if engine == "aimnet2_rxn":
            if any(atom.GetSymbol() not in AIMNET_RXN_ELEMENTS or atom.GetFormalCharge() or atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()):
                raise ValueError("AIMNet2-RXN supports net-neutral, closed-shell HCNO reaction systems here (no formal ions or radicals)")
    if atom_count > 96:
        raise ValueError("The interactive molecular search is limited to 96 explicit atoms; reduce copies or molecule size")


def build_general(experiment):
    validate_molecules(experiment.molecules, experiment.compute.engine)
    combined, names, tags, catalog, instances = Atoms(), [], [], {}, {}
    notes, x_offset = [], 0.0
    charge, has_ions, radicals = 0, False, 0
    detector = SpeciesDetector()
    molecules = [parse_molecule(item.smiles) for item in experiment.molecules]
    # Keep separate inputs outside the detector's bond cutoff, including large
    # metal/halide radii. A fixed organic-sized gap can falsely join salt ions.
    largest_radius = max(covalent_radii[atom.GetAtomicNum()] for mol in molecules for atom in mol.GetAtoms())
    fragment_gap = max(2.8, 2.4 * largest_radius + .3)
    for item_index, (item, molecule) in enumerate(zip(experiment.molecules, molecules)):
        params = AllChem.ETKDGv3()
        params.randomSeed = (experiment.search.seed + item_index) % (2**31 - 1)
        if AllChem.EmbedMolecule(molecule, params) < 0:
            raise ValueError(f"Could not generate 3D coordinates for {item.name}")
        if AllChem.MMFFHasAllMoleculeParams(molecule):
            AllChem.MMFFOptimizeMolecule(molecule, maxIters=200)
        for copy in range(item.copies):
            # Pack each disconnected SMILES fragment separately: salts must not overlap.
            for fragment in Chem.GetMolFrags(molecule, asMols=True):
                coordinates = np.array(fragment.GetConformer().GetPositions())
                coordinates -= coordinates.mean(axis=0)
                coordinates[:, 0] -= coordinates[:, 0].min()
                coordinates[:, 0] += x_offset
                x_offset = coordinates[:, 0].max() + fragment_gap
                atoms = Atoms([a.GetSymbol() for a in fragment.GetAtoms()], positions=coordinates)
                fragment_charge = Chem.GetFormalCharge(fragment)
                fragment_radicals = sum(a.GetNumRadicalElectrons() for a in fragment.GetAtoms())
                charge += fragment_charge
                has_ions |= bool(fragment_charge)
                radicals += fragment_radicals
                atoms.info["chemistry_mode"] = "general"
                for species in detector.identify(atoms):
                    # Names from user inputs are retained only when the entire fragment matches.
                    if len(species.atom_indices) == len(atoms):
                        key = f"{species.formula}:{species.connectivity_hash}"
                        catalog[key] = {"name": item.name, "charge": fragment_charge}
                        global_indices = tuple(range(len(combined), len(combined) + len(atoms)))
                        instances[",".join(map(str, global_indices))] = make_instance_metadata(
                            fragment,
                            global_indices,
                            name=item.name,
                            charge=fragment_charge,
                            connectivity_hash=species.connectivity_hash,
                        )
                tags.extend([len(names)] * len(atoms))
                names.append(item.name)
                combined += atoms
        notes.append(f"{item.name}: {Chem.MolToSmiles(Chem.RemoveHs(molecule))}; {item.amount_mol:g} mol bulk, {item.copies} formula-unit copies; molar mass {Descriptors.MolWt(molecule):.5f} g/mol.")
    combined.set_tags(tags)
    combined.set_cell([max(x_offset + 8, 18), 20, 20])
    combined.center()
    combined.info.update(chemistry_mode="general", total_charge=charge, has_ions=has_ions,
                         radical_electrons=radicals, species_catalog=catalog,
                         species_instances=instances)
    notes.extend(["RDKit ETKDG creates starting geometries. MMFF, when available, is used only for initial preparation, never as a reaction predictor.",
                  "Bulk amounts do not set molecular copy counts or concentration in the finite cluster.",
                  "Aromatic bond templates and explicitly authored SMILES stereochemistry are retained for intact input-derived molecules. R/S and E/Z are reassigned from each sampled 3D geometry.",
                  "Unspecified stereocenters or double-bond stereochemistry remain unresolved; use @/@@ and / or \\ in SMILES when a stereoisomer is known.",
                  "Unknown general-mode fragments are retained as unassigned connectivity, not invented chemical identities.",
                  "Solvent and pressure are not included in these vacuum neural-potential energies; quantitative solution predictions require separate kinetic evidence."])
    return BuiltSystem(combined, names, notes)
