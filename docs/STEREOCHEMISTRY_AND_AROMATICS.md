# Aromatic and stereochemical identity

## What is now preserved

General-mode SMILES inputs carry an atom-mapped RDKit template alongside the ASE coordinates. While an identified component still contains the same mapped atoms and distance-derived connectivity, Chem Sandbox uses that template for:

- aromatic atoms and bond orders;
- canonical and isomeric SMILES;
- formal charge and the user-supplied name;
- tetrahedral `R`/`S` assignment from the current 3D coordinates;
- alkene `E`/`Z` assignment from the current 3D coordinates.

The reaction-graph signature includes this identity. Enantiomers and E/Z isomers therefore occupy different nodes even when their elemental formula and untyped connectivity graph are identical. Reflecting all coordinates flips an explicitly specified tetrahedral stereoisomer and produces a distinct state signature.

Use `@` or `@@` in SMILES for tetrahedral stereochemistry and `/` or `\` for double-bond stereochemistry. If a potentially stereogenic input omits those marks, the UI reports `stereo ?`, keeps the state unresolved, and does not promote the arbitrary ETKDG conformer to an authored chemical identity.

## Deliberate safety boundary

The template is trusted only for an intact input-derived component. If connectivity changes or an atom set no longer matches, the product is shown as an unassigned sampled fragment. The program does not transfer the reactant's aromaticity, bond orders, name, or stereochemistry onto a new product based only on formula.

Distance-derived ASE connectivity has no bond-order information. Consequently, a newly discovered aromatic product cannot yet be certified automatically. That requires bond-order perception or a mapped electronic-structure workflow followed by sanitization and stereochemical validation.

## Identity is not energy splitting

An `R` molecule and its `S` mirror image are different stereochemical identities. In an achiral vacuum they nevertheless have the same energy, apart from negligible parity-violating effects. Rotation/reflection-invariant neural potentials should reproduce that equality. Different energies or rates require a chiral environment, such as an enzyme, chiral catalyst, receptor, solvent field, or another stereoisomeric reactant that is represented by a model within its domain.

The current models do not establish enantioselectivity merely because the graph distinguishes the products. Reliable enantiomeric excess needs stereospecific transition states, free-energy corrections, solvent/catalyst representation, and validated kinetic evidence.

## Validation added in v0.9

Backend tests cover benzene, naphthalene, pyridine, anisole, lactic-acid enantiomers, alanine enantiomers, limonene-like terpene enantiomers, a multi-center glucose/mirror pair, E/Z 2-butene, unspecified lactic acid, coordinate reflection, graph-node separation, and extended-XYZ round trips.

The evaluator-only aromatic benchmark adds:

- benzene + 3 H₂ → cyclohexane, NIST gas-phase reference −206.0 kJ/mol;
- pyridine + 3 H₂ → piperidine, NIST gas-phase reference −193.8 kJ/mol.

All four CUDA engines converged and retained intended connectivity. The diagnostic endpoint deviations are large and model-dependent: AIMNet2 was closest on this two-case sample (−73.4 and −42.2 kJ/mol), while MACE-MP and ANI-2x were far outside quantitative usefulness. The report is at `backend/data/benchmarks/aromatic_suite_v01/latest.md`. These compare electronic model ΔE with experimental ΔH°298 and do not validate mechanisms, catalyst barriers, rates, selectivity, or yield.

Primary identity behavior follows the current [RDKit stereochemistry and aromaticity documentation](https://www.rdkit.org/docs/RDKit_Book.html), including 3D stereochemistry assignment and CIP labels.
