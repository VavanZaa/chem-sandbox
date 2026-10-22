# Operations engine research record

Reviewed 2026-09-22. Primary project documentation was preferred.

| Component | Intended role | Installation decision | Primary source |
|---|---|---|---|
| Cantera 3.2 | Homogeneous kinetics, reactor networks, reacting surfaces, flow and heat transfer | Installed from PyPI in the project venv; Windows wheel supported for Python 3.10–3.14 | [Cantera pip installation](https://www.cantera.org/dev/install/pip.html), [reactor networks](https://cantera.org/stable/reference/reactors/index.html) |
| IDAES 2.13 | Unit operations, thermophysical packages, optimization, and process flowsheets | Installed from PyPI in the project venv; external solvers are a later scoped install | [IDAES optional dependencies](https://idaes-pse.readthedocs.io/en/main/how_to_guides/opt_dependencies.html), [IDAES software](https://idaes.org/software/) |
| RMG-Py | Automated detailed reaction mechanism generation | Deferred to an isolated Linux/WSL or Docker runtime; the official Windows route is Docker or WSL rather than this venv | [RMG installation](https://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/index.html), [kinetics models](https://reactionmechanismgenerator.github.io/RMG-Py/reference/kinetics/index.html) |
| Reaktoro | Chemical equilibrium, aqueous/electrolyte speciation, kinetics, and reactive transport | Deferred to a dedicated Conda/Linux runtime; no compatible direct PyPI package was found | [Reaktoro installation](https://reaktoro.org/installation/index.html), [kinetics tutorial](https://reaktoro.org/v1/tutorials/kinetics/index.html) |
| CatMAP | Heterogeneous catalytic microkinetic modeling | Deferred to a dedicated compatible runtime and only after explicit surface pathways exist | [CatMAP microkinetic tutorial](https://catmap.readthedocs.io/en/latest/tutorials/creating_a_microkinetic_model.html) |
| ASE thermochemistry | Ideal-gas/harmonic corrections and thermochemistry utilities | Already installed through ASE; use only with verified modes and stated assumptions | [ASE thermochemistry](https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html) |
| PLUMED | Enhanced sampling and collective variables | Candidate later integration; requires a compatible molecular dynamics engine | [PLUMED documentation](https://www.plumed.org/doc) |

## Important interpretation

Installing a solver does not create missing science. Cantera needs a kinetic mechanism and thermodynamic data; IDAES needs valid property/reaction packages and solvers; Reaktoro needs an appropriate database/activity model; CatMAP needs surface energetics and an elementary mechanism. Chem Sandbox must carry the origin, units, validity domain, and uncertainty of those inputs.

The current native Windows environment is deliberately limited to packages with a supported, reproducible installation path. Linux-only engines will communicate through versioned artifacts or a local service so dependency changes cannot destabilize discovery and reference data cannot cross the blind-run boundary.
