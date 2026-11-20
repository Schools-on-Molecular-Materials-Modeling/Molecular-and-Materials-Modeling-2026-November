import numpy as np
import torch
from ase.build import molecule
from mace.calculators import mace_off, mace_mp

device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"=== MACE Smoke Test (device={device}) ===")

for label, calc in [
    ("MACE-OFF23 (medium)", mace_off(model="medium", device=device, default_dtype="float64")),
    ("MACE-MP-0 (med+D3)",  mace_mp(model="medium", dispersion=True, device=device, default_dtype="float64")),
]:
    atoms = molecule("H2O")
    atoms.center(vacuum=10.0)
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    print(f"[{label:<19}] PASS | E = {energy:11.6f} eV | max|F| = {np.abs(forces).max():.6f} eV/A")

