"""Load both actual engines and require finite GPU energy/forces."""
import json
import sys
import numpy as np
from ase.build import molecule
from app.models.registry import hardware_status, make_engine

device = sys.argv[1] if len(sys.argv) > 1 else 'cuda'
print(json.dumps(hardware_status(), indent=2), flush=True)
atoms = molecule('H2O')
atoms.info.update(chemistry_mode='general', total_charge=0, radical_electrons=0, has_ions=False)
names = sys.argv[2:] or ['mace_mp', 'ani2x', 'aimnet2', 'aimnet2_rxn']
for name in names:
    engine = make_engine(name, device)
    energy, forces = engine.energy(atoms), engine.forces(atoms)
    assert np.isfinite(energy) and np.isfinite(forces).all()
    assert engine.metadata['device'] == device
    print(json.dumps({'engine': name, 'metadata': engine.metadata, 'energy_ev': energy, 'max_force': float(np.linalg.norm(forces, axis=1).max())}), flush=True)
