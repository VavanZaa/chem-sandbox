# Staged laboratory protocols

Chem Sandbox v0.10 adds a sequential protocol layer for experiments that cannot be represented as one simultaneous reaction network. A protocol passes the material inventory from one stage to the next and supports:

- stage-specific reaction subsets;
- reagent or catalyst addition after earlier reactions have run;
- temperature and volume changes;
- nonreactive holds;
- analytical sample removal;
- imperfect transfers with an explicit recovery fraction;
- elementary mass-action, standard-state Eyring mass-action, and Michaelis–Menten rate laws;
- one material ledger across the complete experiment;
- native SciPy/BDF integration cross-checked with libRoadRunner/CVODE;
- an SBML Level 3 network export validated by libSBML.

The UI exposes this under **Protocol**. The included hexose cycle is a hypothetical execution test. It verifies stage ordering, an intermediate hydrogen addition, sample removal, solver agreement, stereochemical metadata, and atom/charge closure; its rates are not a glucose-process prediction.

## Execution order

Each stage starts with condition changes and declared additions. A `reaction` stage then integrates only its `active_reactions`. `hold` advances time with no active chemistry. `sample` removes a homogeneous fraction of material and volume. `transfer` retains `recovery_fraction` and records the rest as process loss.

Every chemical reaction must conserve atoms and formal charge. Additions and removals are ledger events rather than reaction stoichiometry. The final check is:

`initial + additions − samples − transfer losses = final inventory`

for every tracked element and charge.

## Kinetic evidence

`mass_action` accepts a dimension-checked rate constant. A reaction of total reactant order `n` requires units equivalent to `L^(n−1) mol^(1−n) s−1`.

`eyring_mass_action` accepts activation free energy, reference temperature, and optional activation entropy. It uses a 1 mol/L standard-state activity convention. These values still require appropriate solvent, catalyst, pH, and standard-state evidence.

`michaelis_menten` requires an explicitly declared substrate and catalyst pool plus `kcat` and `Km`. The catalyst is not consumed by stoichiometry. This empirical law does not establish a molecular mechanism and should not be extrapolated beyond its evidence range.

## Stereochemistry and enzymes

Small-molecule SMILES retain their explicit stereochemical identity. An enzyme or chiral catalyst can be added in a middle stage, but the software does not infer stereoselectivity from its name. Competing R/S channels need independently supplied stereospecific rate evidence. Composition-only macromolecular species must include an identity note because an elemental formula does not define an enzyme, sequence, fold, or active site.

## Krebs-cycle and metabolic models

A real citric-acid-cycle model needs the enzyme isoforms, cofactors, protonation states, compartments, reversible reactions, inhibition/activation, pH, ionic strength, and condition-specific parameters. The generic protocol engine can execute such a declared network, but the application does not silently insert a textbook Krebs recipe or fabricate missing constants.

SBML import inspection is available at `POST /api/protocols/sbml/inspect`. It validates and summarizes an SBML model without automatically trusting or executing arbitrary rate expressions. The native protocol endpoint is `POST /api/protocols/simulate`.

The exported SBML currently contains the reaction network and rate laws; the ordered laboratory stage schedule remains in the protocol JSON. RoadRunner cross-checking applies that schedule piecewise to the SBML network. SBML itself supports events, and encoding the full schedule as portable SBML Events is a future interoperability step.

## Scientific firewall

Protocol definitions are downstream, user-declared operating knowledge. They never alter candidate generation in blind discovery. Known pathways, expected products, and literature recipes can be used for post-discovery validation or an explicitly declared protocol, but are not available to the internal product search.

Tool references: [libSBML Python API](https://sbml.org/software/libsbml/libsbml-docs/api/python/), [libRoadRunner simulation API](https://sys-bio.github.io/roadrunner/docs-build/PythonAPIReference/cls_RoadRunner.html), and [SciPy `solve_ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html).
