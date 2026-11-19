import shutil
import torch
from pathlib import Path
from ase.build import molecule
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.units import kJ, mol, invcm
from ase.calculators.mopac import MOPAC
from xtb.ase.calculator import XTB
from mace.calculators import mace_off

# Species specifications: (ase_build_name, symmetrynumber)
# Note: trans-buta-1,3-diene is C2h (sym=2), C2H4 is D2h (sym=4), half-chair cyclohexene is C2 (sym=2)
SPECIES = {
    "C4H6": ("trans-buta-1,3-diene", 2),
    "C2H4": ("C2H4", 4),
    "C6H10": ("cyclohexene", 2),
}

def get_thermochemistry(atoms, e_pot, sym, vib_prefix, use_thermo=True):
    vib_dir = Path(vib_prefix)
    if vib_dir.exists():
        shutil.rmtree(vib_dir)

    vib = Vibrations(atoms, name=str(vib_dir), delta=0.01)
    vib.run()
    vib_energies = vib.get_energies()
    zpe = vib.get_zero_point_energy()
    vib.clean()
    if vib_dir.exists():
        shutil.rmtree(vib_dir)

    internal_energies = vib_energies[6:]  # All three species are nonlinear
    if not use_thermo:
        return zpe, e_pot  # MOPAC PM7 dHf is already parameterized at 298.15 K

    thermo = IdealGasThermo(
        vib_energies=internal_energies,
        potentialenergy=e_pot,
        atoms=atoms,
        geometry="nonlinear",
        symmetrynumber=sym,
        spin=0,
    )
    h_298 = thermo.get_enthalpy(temperature=298.15, verbose=False)
    return zpe, h_298

def run_method(method_name, make_calc_fn, is_mopac=False):
    print(f"\n--- Running {method_name} ---")
    data = {}
    for label, (mol_name, sym) in SPECIES.items():
        atoms = molecule(mol_name)
        atoms.pbc = False
        if not is_mopac:
            atoms.center(vacuum=12.0)

        calc = make_calc_fn(f"tmp_{method_name}_{label}")
        atoms.calc = calc
        opt = BFGS(atoms, logfile=None)
        opt.run(fmax=0.001)

        if is_mopac:
            e_val = calc.get_final_heat_of_formation()
        else:
            e_val = atoms.get_potential_energy()

        zpe, h_298 = get_thermochemistry(
            atoms,
            e_val,
            sym,
            f"vib_{method_name}_{label}",
            use_thermo=(not is_mopac),
        )
        data[label] = {"e_opt": e_val, "zpe": zpe, "h_298": h_298}
        print(f"  {label:<6} ({mol_name:<20}) | steps={opt.nsteps:<3} | ZPE={zpe:.4f} eV")

    dE_kj = (data["C6H10"]["e_opt"] - (data["C4H6"]["e_opt"] + data["C2H4"]["e_opt"])) / (kJ / mol)
    dZPE_kj = (data["C6H10"]["zpe"] - (data["C4H6"]["zpe"] + data["C2H4"]["zpe"])) / (kJ / mol)
    dH_kj = (data["C6H10"]["h_298"] - (data["C4H6"]["h_298"] + data["C2H4"]["h_298"])) / (kJ / mol)
    return dE_kj, dZPE_kj, dH_kj

# 1. MOPAC (PM7)
res_mopac = run_method(
    "MOPAC_PM7",
    lambda lbl: MOPAC(
        label=lbl,
        method="PM7",
        task="1SCF GRADIENTS PRECISE",
        command="mopac PREFIX.mop > /dev/null 2>&1",
    ),
    is_mopac=True,
)

# 2. GFN2-xTB
res_xtb = run_method(
    "GFN2_xTB",
    lambda lbl: XTB(method="GFN2-xTB"),
    is_mopac=False,
)

# 3. MACE-OFF23 (medium)
device = "cuda" if torch.cuda.is_available() else "cpu"
mace_calc = mace_off(model="medium", device=device, default_dtype="float64")
res_mace = run_method(
    "MACE_OFF23",
    lambda lbl: mace_calc,
    is_mopac=False,
)

# Clean up MOPAC scratch files
for ext in ("*.mop", "*.out", "*.arc"):
    for f in Path(".").glob(ext):
        f.unlink()

# Summary Comparison Table
print("\n=== Diels-Alder Reaction Enthalpy Summary: C4H6 + C2H4 -> C6H10 (298.15 K) ===")
print(f"{'Method':<18} | {'dE_elec (kJ/mol)':>16} | {'dZPE (kJ/mol)':>14} | {'dH_298K (kJ/mol)':>16}")
print("-" * 73)
print(f"{'MOPAC (PM7)':<18} | {'(implicit)':>16} | {res_mopac[1]:+14.2f} | {res_mopac[2]:+16.2f}")
print(f"{'GFN2-xTB':<18} | {res_xtb[0]:+16.2f} | {res_xtb[1]:+14.2f} | {res_xtb[2]:+16.2f}")
print(f"{'MACE-OFF23 (med)':<18} | {res_mace[0]:+16.2f} | {res_mace[1]:+14.2f} | {res_mace[2]:+16.2f}")
print(f"{'Experiment (NIST)':<18} | {'—':>16} | {'~ +22.0':>14} | {-166.70:+16.2f}")

