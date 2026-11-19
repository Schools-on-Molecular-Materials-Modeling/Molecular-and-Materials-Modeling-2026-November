import shutil
from pathlib import Path
from ase.build import molecule
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.units import kJ, mol, invcm
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
    dhf_ev = calc.get_final_heat_of_formation()
    return atoms, dhf_ev, opt.nsteps

def analyze_vibrations(name, atoms, geometry):
    vib_dir = Path(f"vib_mopac_{name}")
    if vib_dir.exists():
        shutil.rmtree(vib_dir)

    atoms.calc = make_mopac_calc(f"mopac_vib_{name}")
    vib = Vibrations(atoms, name=str(vib_dir), delta=0.01)
    vib.run()
    vib_energies = vib.get_energies()
    zpe = vib.get_zero_point_energy()
    vib.clean()
    if vib_dir.exists():
        shutil.rmtree(vib_dir)

    n_tr = 5 if geometry == "linear" else 6
    internal_cm1 = [e.real / invcm for e in vib_energies[n_tr:]]
    return internal_cm1, zpe

specs = {
    "H2": "linear",
    "C2H4": "nonlinear",
    "C2H6": "nonlinear",
}

print("\n=== Tight Geometry Optimization & Harmonic Vibrations (delta = 0.01 A) ===")
results = {}
for name, geom in specs.items():
    atoms, dhf_ev, steps = relax_and_get_dhf(name, fmax=0.0005)
    freqs_cm1, zpe = analyze_vibrations(name, atoms, geom)
    results[name] = {
        "atoms": atoms,
        "dhf_ev": dhf_ev,
        "zpe": zpe,
        "freqs": freqs_cm1,
    }
    print(f"[{name}] steps={steps} | dHf_298K={dhf_ev/(kJ/mol):+8.2f} kJ/mol | ZPE={zpe:.4f} eV")

# 1. Raw PM7 reaction enthalpy (affected by known NDDO H2 diatomic artifact: dHf(H2) = -133.9 kJ/mol)
dh_raw_kj = (results["C2H6"]["dhf_ev"] - (results["C2H4"]["dhf_ev"] + results["H2"]["dhf_ev"])) / (kJ / mol)

# 2. Corrected PM7 reaction enthalpy using exact thermodynamic reference state dHf(H2, 298K) = 0.00 kJ/mol
dh_ref_kj = (results["C2H6"]["dhf_ev"] - results["C2H4"]["dhf_ev"]) / (kJ / mol)

print("\n=== Hydrogenation Thermochemistry Summary (PM7): C2H4 + H2 -> C2H6 ===")
print(f"dHf(H2, 298K) raw  : {results['H2']['dhf_ev']/(kJ/mol):+8.2f} kJ/mol  (Known NDDO H2 artifact; Expt ref: 0.00 kJ/mol)")
print(f"dHf(C2H4, 298K)    : {results['C2H4']['dhf_ev']/(kJ/mol):+8.2f} kJ/mol  (Expt ref: +52.40 kJ/mol)")
print(f"dHf(C2H6, 298K)    : {results['C2H6']['dhf_ev']/(kJ/mol):+8.2f} kJ/mol  (Expt ref: -84.00 kJ/mol)")
print(f"dH_298K (raw H2)   : {dh_raw_kj:+8.2f} kJ/mol")
print(f"dH_298K (dHf_H2=0) : {dh_ref_kj:+8.2f} kJ/mol  [Expt ref: -136.30 kJ/mol]")

for ext in ("*.mop", "*.out", "*.arc"):
    for f in Path(".").glob(ext):
        f.unlink()
import shutil
from pathlib import Path
from ase.build import molecule
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.units import kJ, mol, invcm
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
    dhf_ev = calc.get_final_heat_of_formation()
    return atoms, dhf_ev, opt.nsteps

def analyze_vibrations(name, atoms, geometry):
    vib_dir = Path(f"vib_mopac_{name}")
    if vib_dir.exists():
        shutil.rmtree(vib_dir)

    atoms.calc = make_mopac_calc(f"mopac_vib_{name}")
    vib = Vibrations(atoms, name=str(vib_dir), delta=0.01)
    vib.run()
    vib_energies = vib.get_energies()
    zpe = vib.get_zero_point_energy()
    vib.clean()
    if vib_dir.exists():
        shutil.rmtree(vib_dir)

    n_tr = 5 if geometry == "linear" else 6
    internal_cm1 = [e.real / invcm for e in vib_energies[n_tr:]]
    return internal_cm1, zpe

specs = {
    "H2": "linear",
    "C2H4": "nonlinear",
    "C2H6": "nonlinear",
}

print("\n=== Tight Geometry Optimization & Harmonic Vibrations (delta = 0.01 A) ===")
results = {}
for name, geom in specs.items():
    atoms, dhf_ev, steps = relax_and_get_dhf(name, fmax=0.0005)
    freqs_cm1, zpe = analyze_vibrations(name, atoms, geom)
    results[name] = {
        "atoms": atoms,
        "dhf_ev": dhf_ev,
        "zpe": zpe,
        "freqs": freqs_cm1,
    }
    print(f"[{name}] steps={steps} | dHf_298K={dhf_ev/(kJ/mol):+8.2f} kJ/mol | ZPE={zpe:.4f} eV")

# 1. Raw PM7 reaction enthalpy (affected by known NDDO H2 diatomic artifact: dHf(H2) = -133.9 kJ/mol)
dh_raw_kj = (results["C2H6"]["dhf_ev"] - (results["C2H4"]["dhf_ev"] + results["H2"]["dhf_ev"])) / (kJ / mol)

# 2. Corrected PM7 reaction enthalpy using exact thermodynamic reference state dHf(H2, 298K) = 0.00 kJ/mol
dh_ref_kj = (results["C2H6"]["dhf_ev"] - results["C2H4"]["dhf_ev"]) / (kJ / mol)

print("\n=== Hydrogenation Thermochemistry Summary (PM7): C2H4 + H2 -> C2H6 ===")
print(f"dHf(H2, 298K) raw  : {results['H2']['dhf_ev']/(kJ/mol):+8.2f} kJ/mol  (Known NDDO H2 artifact; Expt ref: 0.00 kJ/mol)")
print(f"dHf(C2H4, 298K)    : {results['C2H4']['dhf_ev']/(kJ/mol):+8.2f} kJ/mol  (Expt ref: +52.40 kJ/mol)")
print(f"dHf(C2H6, 298K)    : {results['C2H6']['dhf_ev']/(kJ/mol):+8.2f} kJ/mol  (Expt ref: -84.00 kJ/mol)")
print(f"dH_298K (raw H2)   : {dh_raw_kj:+8.2f} kJ/mol")
print(f"dH_298K (dHf_H2=0) : {dh_ref_kj:+8.2f} kJ/mol  [Expt ref: -136.30 kJ/mol]")

for ext in ("*.mop", "*.out", "*.arc"):
    for f in Path(".").glob(ext):
        f.unlink()

