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
print("Reaction           : CH4(g) + 2 O2(g) -> CO2(g) + 2 H2O(g)")
print("Energy Definition  : FINAL HEAT OF FORMATION dHf(298.15 K)")

def make_mopac_calc(label, is_triplet=False):
    # Ground-state O2 is a spin triplet (2 unpaired electrons -> TRIPLET UHF in MOPAC)
    spin_kw = " TRIPLET UHF" if is_triplet else ""
    return MOPAC(
        label=label,
        method=method,
        task=f"1SCF GRADIENTS PRECISE{spin_kw}",
        command="mopac PREFIX.mop > /dev/null 2>&1",
    )

def relax_and_get_dhf(name, is_triplet=False, fmax=0.0005):
    atoms = molecule(name)
    atoms.pbc = False
    calc = make_mopac_calc(f"mopac_opt_{name}", is_triplet=is_triplet)
    atoms.calc = calc
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=fmax)
    dhf_ev = calc.get_final_heat_of_formation()
    return atoms, dhf_ev, opt.nsteps

def analyze_vibrations(name, atoms, geometry, is_triplet=False):
    vib_dir = Path(f"vib_mopac_{name}")
    if vib_dir.exists():
        shutil.rmtree(vib_dir)

    atoms.calc = make_mopac_calc(f"mopac_vib_{name}", is_triplet=is_triplet)
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
    "CH4": ("nonlinear", False),
    "O2":  ("linear",    True),
    "CO2": ("linear",    False),
    "H2O": ("nonlinear", False),
}

print("\n=== Tight Geometry Optimization & Harmonic Vibrations (delta = 0.01 A) ===")
results = {}
for name, (geom, is_triplet) in specs.items():
    atoms, dhf_ev, steps = relax_and_get_dhf(name, is_triplet=is_triplet, fmax=0.0005)
    freqs_cm1, zpe = analyze_vibrations(name, atoms, geom, is_triplet=is_triplet)
    results[name] = {
        "atoms": atoms,
        "dhf_ev": dhf_ev,
        "zpe": zpe,
        "freqs": freqs_cm1,
    }
    state_label = "triplet" if is_triplet else "singlet"
    print(f"[{name:<3} ({state_label})] steps={steps:<2} | dHf_298K={dhf_ev/(kJ/mol):+8.2f} kJ/mol | ZPE={zpe:.4f} eV")

# Key relaxed bond lengths
d_ch = results["CH4"]["atoms"].get_distance(0, 1)
d_oo = results["O2"]["atoms"].get_distance(0, 1)
d_co = results["CO2"]["atoms"].get_distance(0, 1)
d_oh = results["H2O"]["atoms"].get_distance(0, 1)

print("\n=== Relaxed Bond Lengths vs. Experiment ===")
print(f"C-H in CH4         : {d_ch:.4f} A  (Expt: 1.087 A)")
print(f"O=O in O2 (triplet): {d_oo:.4f} A  (Expt: 1.208 A)")
print(f"C=O in CO2         : {d_co:.4f} A  (Expt: 1.160 A)")
print(f"O-H in H2O         : {d_oh:.4f} A  (Expt: 0.958 A)")

# Reaction enthalpy: CH4(g) + 2 O2(g) -> CO2(g) + 2 H2O(g)
dh_ch4 = results["CH4"]["dhf_ev"] / (kJ / mol)
dh_o2  = results["O2"]["dhf_ev"] / (kJ / mol)
dh_co2 = results["CO2"]["dhf_ev"] / (kJ / mol)
dh_h2o = results["H2O"]["dhf_ev"] / (kJ / mol)

dc_H_raw = (dh_co2 + 2.0 * dh_h2o) - (dh_ch4 + 2.0 * dh_o2)
dc_H_ref = (dh_co2 + 2.0 * dh_h2o) - (dh_ch4 + 2.0 * 0.0)

print("\n=== Gas-Phase Methane Combustion Summary (PM7): CH4 + 2 O2 -> CO2 + 2 H2O ===")
print(f"dHf(CH4, 298K)     : {dh_ch4:+8.2f} kJ/mol  (Expt ref:  -74.87 kJ/mol)")
print(f"dHf(O2, 298K) raw  : {dh_o2:+8.2f} kJ/mol  (Expt ref:    0.00 kJ/mol)")
print(f"dHf(CO2, 298K)     : {dh_co2:+8.2f} kJ/mol  (Expt ref: -393.51 kJ/mol)")
print(f"dHf(H2O, 298K)     : {dh_h2o:+8.2f} kJ/mol  (Expt ref: -241.83 kJ/mol)")
print(f"dcH_298K (raw O2)  : {dc_H_raw:+8.2f} kJ/mol  ({dc_H_raw * (kJ / mol):+.4f} eV)")
print(f"dcH_298K (dHf_O2=0): {dc_H_ref:+8.2f} kJ/mol  ({dc_H_ref * (kJ / mol):+.4f} eV)  [Expt ref: -802.30 kJ/mol]")

# Clean up temporary MOPAC scratch files
for ext in ("*.mop", "*.out", "*.arc"):
    for f in Path(".").glob(ext):
        f.unlink()

