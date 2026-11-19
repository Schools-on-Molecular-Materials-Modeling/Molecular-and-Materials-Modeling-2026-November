import os
import shutil
import subprocess
from pathlib import Path
from ase.build import molecule

print("=================================================================")
print("  Exercise 01: Computational Apparatus & ASE Environment Check")
print("=================================================================")

# 1. Check Python Environment & GPU Hardware
import numpy as np
import ase
import torch

print(f"\n[1] Core Python & Hardware")
print(f"  ASE version        : {ase.__version__}")
print(f"  NumPy version      : {np.__version__}")
print(f"  PyTorch version    : {torch.__version__}")
print(f"  CUDA available     : {torch.cuda.is_available()}", end="")
if torch.cuda.is_available():
    print(f" ({torch.cuda.get_device_name(0)})")
else:
    print(" (running on CPU)")

# 2. Check External Quantum Chemistry & GAP Binaries
print(f"\n[2] External Binaries in PATH")
binaries = ["mopac", "nwchem", "mpirun", "gap_fit", "quip"]
for b in binaries:
    path = shutil.which(b)
    status = f"[OK] {path}" if path else "[MISSING]"
    print(f"  {b:<10} : {status}")

# 3. Single-Point Energy & Force Sanity Check on H2O across all Engines
print(f"\n[3] Calculator Smoke Tests (Single-Point H2O)")
results = []

def run_smoke_test(name, setup_fn, cleanup_patterns=()):
    try:
        atoms = molecule("H2O")
        atoms.pbc = False
        calc = setup_fn(atoms)
        atoms.calc = calc
        e = atoms.get_potential_energy()
        f = atoms.get_forces()
        fmax = float(np.abs(f).max())
        results.append((name, "PASS", f"E = {e:10.4f} eV | max|F| = {fmax:.4f} eV/A"))
    except Exception as exc:
        results.append((name, "FAIL", str(exc).splitlines()[0]))
    finally:
        for pat in cleanup_patterns:
            for p in Path(".").glob(pat):
                if p.is_file():
                    p.unlink()

# A. MOPAC (PM7)
def setup_mopac(atoms):
    from ase.calculators.mopac import MOPAC
    return MOPAC(
        label="test_mopac_h2o",
        method="PM7",
        task="1SCF GRADIENTS",
        command="mopac PREFIX.mop > /dev/null 2>&1",
    )
run_smoke_test("MOPAC (PM7)", setup_mopac, cleanup_patterns=("test_mopac_h2o.*",))

# B. GFN2-xTB
def setup_xtb(atoms):
    from xtb.ase.calculator import XTB
    return XTB(method="GFN2-xTB")
run_smoke_test("GFN2-xTB", setup_xtb)

# C. MACE-OFF23 (Organic MLIP)
def setup_mace_off(atoms):
    from mace.calculators import mace_off
    atoms.center(vacuum=10.0)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    return mace_off(model="medium", device=dev, default_dtype="float64")
run_smoke_test("MACE-OFF23 (med)", setup_mace_off)

# D. MACE-MP-0 (Universal Materials MLIP + D3)
def setup_mace_mp(atoms):
    from mace.calculators import mace_mp
    atoms.center(vacuum=10.0)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    return mace_mp(model="medium", dispersion=True, device=dev, default_dtype="float64")
run_smoke_test("MACE-MP-0 (med+D3)", setup_mace_mp)

# E. CHGNet (Universal GNN MLIP)
def setup_chgnet(atoms):
    from chgnet.model.dynamics import CHGNetCalculator
    atoms.center(vacuum=10.0)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    return CHGNetCalculator(use_device=dev)
run_smoke_test("CHGNet", setup_chgnet)

# F. ALIGNN-FF (JARVIS Foundation MLIP)
def setup_alignn(atoms):
    from alignn.ff.ff import AlignnAtomwiseCalculator, default_path
    atoms.center(vacuum=10.0)
    return AlignnAtomwiseCalculator(path=default_path())
run_smoke_test("ALIGNN-FF", setup_alignn)

# G. NWChem (DFT B3LYP/6-31G*)
def setup_nwchem(atoms):
    from ase.calculators.nwchem import NWChem
    return NWChem(
        label="test_nwchem_h2o",
        theory="dft",
        basis="6-31G*",
        dft={"xc": "b3lyp", "mult": 1},
        command="mpirun -np 2 nwchem PREFIX.nwi > PREFIX.nwo",
    )
run_smoke_test("NWChem (B3LYP)", setup_nwchem, cleanup_patterns=("test_nwchem_h2o*",))

# H. QUIP / GAP (soap_turbo check on C2H4 if CH.xml is present)
def setup_quip(atoms):
    from quippy.potential import Potential
    # Check local or relative CH.xml path if available; otherwise verify import only
    candidates = list(Path(".").rglob("CH.xml"))
    if not candidates:
        raise RuntimeError("quippy imported OK (no local CH.xml potential file in cwd)")
    atoms[:] = molecule("C2H4")
    atoms.center(vacuum=10.0)
    return Potential(param_filename=str(candidates[0].resolve()))
run_smoke_test("QUIP / GAP", setup_quip)

for name, status, info in results:
    print(f"  {name:<19} | {status:<4} | {info}")

print("\n=================================================================")

