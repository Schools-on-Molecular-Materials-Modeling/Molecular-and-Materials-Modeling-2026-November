import shutil
import torch
from pathlib import Path
from ase import Atoms
from ase.build import molecule
from ase.collections import g2
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.units import kJ, mol
from ase.calculators.mopac import MOPAC
from xtb.ase.calculator import XTB
from mace.calculators import mace_off

def build_species(name):
    if name == "C4H6":
        # ASE G2 collection name for trans-1,3-butadiene is "butadiene"
        return molecule("butadiene")
    elif name == "C2H4":
        return molecule("C2H4")
    elif name == "C6H10":
        if "cyclohexene" in g2.names:
            return molecule("cyclohexene")
        # Explicit half-chair cyclohexene (C6H10) starting coordinates (in Angstroms)
        symbols = ["C"] * 6 + ["H"] * 10
        positions = [
            [ 0.000,  0.668,  1.258],  # C1 (sp2)
            [ 0.000, -0.668,  1.258],  # C2 (sp2)
            [-0.145, -1.495, -0.002],  # C3 (sp3)
            [ 0.258, -0.748, -1.275],  # C4 (sp3)
            [-0.258,  0.748, -1.275],  # C5 (sp3)
            [ 0.145,  1.495, -0.002],  # C6 (sp3)
            [ 0.032,  1.225,  2.193],  # H on C1
            [ 0.032, -1.225,  2.193],  # H on C2
            [ 0.453, -2.411,  0.088],  # H on C3
            [-1.193, -1.825, -0.082],  # H on C3
            [ 1.351, -0.822, -1.352],  # H on C4
            [-0.145, -1.262, -2.155],  # H on C4
            [-1.351,  0.822, -1.352],  # H on C5
            [ 0.145,  1.262, -2.155],  # H on C5
            [ 1.193,  1.825, -0.082],  # H on C6
            [-0.453,  2.411,  0.088],  # H on C6
        ]
        return Atoms(symbols=symbols, positions=positions)
    raise ValueError(f"Unknown species: {name}")

# Species specifications: (display_name, symmetrynumber)
SPECIES = {
    "C4H6":  ("1,3-butadiene", 2),
    "C2H4":  ("ethene",        4),
    "C6H10": ("cyclohexene",   2),
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
    for label, (disp_name, sym) in SPECIES.items():
        atoms = build_species(label)
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
        print(f"  {label:<6} ({disp_name:<14}) | steps={opt.nsteps:<3} | ZPE={zpe:.4f} eV")

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
print(f"{'MOPAC (PM7)':<18} | {'(implicit)':>16} | {res_mopac:+14.2f} | {res_mopac:+16.2f}")
print(f"{'GFN2-xTB':<18} | {res_xtb:+16.2f} | {res_xtb:+14.2f} | {res_xtb:+16.2f}")
print(f"{'MACE-OFF23 (med)':<18} | {res_mace:+16.2f} | {res_mace:+14.2f} | {res_mace:+16.2f}")
print(f"{'Experiment (NIST)':<18} | {'—':>16} | {'~ +22.0':>14} | {-166.70:+16.2f}")

