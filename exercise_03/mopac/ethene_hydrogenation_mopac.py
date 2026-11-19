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
print("Energy Definition  : FINAL HEAT OF FORMATION dHf(298.15 K) in eV")

def make_mopac_calc(label):
    return MOPAC(
        label=label,
        method=method,
        task="1SCF GRADIENTS PRECISE",
    )

def relax_and_get_dhf(name, fmax=0.0005):
    atoms = molecule(name)
    atoms.pbc = False
    calc = make_mopac_calc(f"mopac_opt_{name}")
    atoms.calc = calc
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=fmax)
    # Explicitly read FINAL HEAT OF FORMATION (298.15 K) in eV
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
    print(f"\n[{name}] (BFGS steps={steps}, dHf_298K={dhf_ev:+.6f} eV [{dhf_ev/(kJ/mol):+.2f} kJ/mol], ZPE={zpe:.4f} eV)")
    print("  Internal frequencies (cm^-1):", ", ".join(f"{f:.1f}" for f in freqs_cm1))

d_hh = results["H2"]["atoms"].get_distance(0, 1)
d_cc_ene = results["C2H4"]["atoms"].get_distance(0, 1)
d_cc_ane = results["C2H6"]["atoms"].get_distance(0, 1)

print("\n=== Relaxed Bond Lengths vs. Experiment ===")
print(f"H-H  in H2         : {d_hh:.4f} A  (Expt: 0.741 A)")
print(f"C=C  in C2H4       : {d_cc_ene:.4f} A  (Expt: 1.339 A)")
print(f"C-C  in C2H6       : {d_cc_ane:.4f} A  (Expt: 1.535 A)")

# Direct reaction enthalpy at 298.15 K from PM7 Heats of Formation (no double-counted ZPE)
dH_298_ev = results["C2H6"]["dhf_ev"] - (results["C2H4"]["dhf_ev"] + results["H2"]["dhf_ev"])
dZPE_ev = results["C2H6"]["zpe"] - (results["C2H4"]["zpe"] + results["H2"]["zpe"])

print("\n=== Hydrogenation Thermochemistry Summary (PM7): C2H4 + H2 -> C2H6 ===")
print(f"dHf(H2, 298K)      : {results['H2']['dhf_ev']/(kJ/mol):+8.2f} kJ/mol  (Expt ref:   0.00 kJ/mol)")
print(f"dHf(C2H4, 298K)    : {results['C2H4']['dhf_ev']/(kJ/mol):+8.2f} kJ/mol  (Expt ref: +52.40 kJ/mol)")
print(f"dHf(C2H6, 298K)    : {results['C2H6']['dhf_ev']/(kJ/mol):+8.2f} kJ/mol  (Expt ref: -84.00 kJ/mol)")
print(f"Harmonic dZPE      : {dZPE_ev/(kJ/mol):+8.2f} kJ/mol  (Informational only; already implicit in PM7 dHf)")
print(f"Enthalpy dH_298K   : {dH_298_ev:8.4f} eV  ({dH_298_ev/(kJ/mol):8.2f} kJ/mol)  [Expt ref: -1.41 eV (-136.3 kJ/mol)]")

for ext in ("*.mop", "*.out", "*.arc"):
    for f in Path(".").glob(ext):
        f.unlink()

