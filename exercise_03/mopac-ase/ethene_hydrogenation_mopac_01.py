from pathlib import Path
from ase.build import molecule
from ase.optimize import BFGS
from ase.units import kJ, mol
from ase.calculators.mopac import MOPAC

method = "PM7"

print("=== MOPAC Calculator Initialization ===")
print(f"Method             : {method} (Semi-Empirical NDDO Hamiltonian)")
print("Energy Definition  : FINAL HEAT OF FORMATION dHf(298.15 K)")


def make_mopac_calc(label):
    # Redirect Fortran binary stdout/stderr to /dev/null so ASE runs quietly
    return MOPAC(
        label=label,
        method=method,
        task="1SCF GRADIENTS PRECISE",
        command="mopac PREFIX.mop > /dev/null 2>&1",
    )


def relax_and_get_dhf(name, fmax=0.0005):
    atoms = molecule(name)
    atoms.pbc = False
    calc = make_mopac_calc(f"mopac_opt_{name}")
    atoms.calc = calc
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=fmax)
    # MOPAC's FINAL HEAT OF FORMATION (298.15 K) is already a formation enthalpy,
    # so no vibrational/thermal analysis is needed here.
    dhf_ev = calc.get_final_heat_of_formation()
    return dhf_ev, opt.nsteps


species = ["H2", "C2H4", "C2H6"]

print("\n=== Geometry Optimization & MOPAC Heat of Formation ===")
dhf = {}  # kJ/mol
for name in species:
    dhf_ev, steps = relax_and_get_dhf(name, fmax=0.0005)
    dhf[name] = dhf_ev / (kJ / mol)
    print(f"[{name}] steps={steps} | dHf_298K={dhf[name]:+8.2f} kJ/mol")

# 1. Raw PM7 reaction enthalpy (affected by the NDDO H2 artifact: dHf(H2) = -133.9 kJ/mol)
dh_raw = dhf["C2H6"] - (dhf["C2H4"] + dhf["H2"])

# 2. Reaction enthalpy using the exact reference state dHf(H2, 298 K) = 0.00 kJ/mol
dh_ref = dhf["C2H6"] - dhf["C2H4"]

print("\n=== Hydrogenation Thermochemistry Summary (PM7): C2H4 + H2 -> C2H6 ===")
print(f"dHf(H2, 298K) raw  : {dhf['H2']:+8.2f} kJ/mol  (Known NDDO H2 artifact; Expt ref: 0.00 kJ/mol)")
print(f"dHf(C2H4, 298K)    : {dhf['C2H4']:+8.2f} kJ/mol  (Expt ref: +52.40 kJ/mol)")
print(f"dHf(C2H6, 298K)    : {dhf['C2H6']:+8.2f} kJ/mol  (Expt ref: -84.00 kJ/mol)")
print(f"dH_298K (raw H2)   : {dh_raw:+8.2f} kJ/mol")
print(f"dH_298K (dHf_H2=0) : {dh_ref:+8.2f} kJ/mol  [Expt ref: -136.30 kJ/mol]")

# Clean up MOPAC scratch files
for ext in ("*.mop", "*.out", "*.arc", "*.aux"):
    for f in Path(".").glob(ext):
        f.unlink()
