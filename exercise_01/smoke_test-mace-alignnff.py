#!/usr/bin/env python
"""
smoke_test-mace-alignnff.py

Minimal runtime smoke test for MACE and ALIGNN-FF in the same conda env.
Builds a small periodic Cu cell, evaluates energy + forces with each calculator,
and prints a pass/fail line. Run from the exercise_01 directory:

    conda activate mace_env
    export TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1   # already set via activate.d if you followed earlier step
    python smoke_test-mace-alignnff.py
"""

import sys
import numpy as np
from ase.build import bulk

# ----------------------------------------------------------------------
# Shared test system: 2x2x2 fcc Cu (16 atoms) — small, fast, periodic.
# ----------------------------------------------------------------------
atoms_template = bulk("Cu", "fcc", a=3.6).repeat((2, 2, 2))
print(f"[info] Test system: {len(atoms_template)} atoms, Cu fcc 2x2x2, a=3.6 Å\n")


def report(name, atoms):
    """Print energy and max force magnitude; return True if finite."""
    try:
        e = atoms.get_potential_energy()
        f = atoms.get_forces()
        fmax = np.abs(f).max()
        ok = np.isfinite(e) and np.isfinite(f).all()
        print(f"[{name}] energy = {e: .6f} eV   max|F| = {fmax:.6f} eV/Å   "
              f"{'OK' if ok else 'NON-FINITE'}")
        return ok
    except Exception as exc:
        print(f"[{name}] FAILED: {type(exc).__name__}: {exc}")
        return False


results = {}

# ----------------------------------------------------------------------
# 1) MACE — using the built-in mace_mp helper (downloads/caches model)
# ----------------------------------------------------------------------
print("=" * 70)
print("MACE (mace_mp, small model, CPU)")
print("=" * 70)
try:
    from mace.calculators import mace_mp

    # 'small' is the lightweight MACE-MP-0 checkpoint; use device='cuda' if available
    mace_calc = mace_mp(model="small", device="cpu", default_dtype="float32")
    atoms_mace = atoms_template.copy()
    atoms_mace.calc = mace_calc
    results["MACE"] = report("MACE", atoms_mace)
except Exception as exc:
    print(f"[MACE] import/setup FAILED: {type(exc).__name__}: {exc}")
    results["MACE"] = False

print()

# ----------------------------------------------------------------------
# 2) ALIGNN-FF — uses AlignnAtomwiseCalculator + default_path (v2026.11.1)
# ----------------------------------------------------------------------
print("=" * 70)
print("ALIGNN-FF (AlignnAtomwiseCalculator, default_path, CPU)")
print("=" * 70)
try:
    from alignn.ff.ff import AlignnAtomwiseCalculator, default_path

    alignn_calc = AlignnAtomwiseCalculator(path=default_path(), device="cpu")
    atoms_alignn = atoms_template.copy()
    atoms_alignn.calc = alignn_calc
    results["ALIGNN-FF"] = report("ALIGNN-FF", atoms_alignn)
except Exception as exc:
    print(f"[ALIGNN-FF] import/setup FAILED: {type(exc).__name__}: {exc}")
    results["ALIGNN-FF"] = False

# ----------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
all_ok = True
for name, ok in results.items():
    status = "PASS" if ok else "FAIL"
    print(f"  {name:12s} : {status}")
    all_ok &= ok

sys.exit(0 if all_ok else 1)
