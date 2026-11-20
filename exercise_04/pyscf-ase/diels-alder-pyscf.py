import shutil
from pathlib import Path

import numpy as np
from ase import Atoms, units
from ase.build import molecule
from ase.calculators.calculator import Calculator, all_changes
from ase.collections import g2
from ase.io import read, write
from ase.optimize import BFGS
from ase.thermochemistry import IdealGasThermo
from ase.units import kJ, mol, invcm
from ase.vibrations import Vibrations

from pyscf import gto, dft, lib

# D3(BJ) dispersion needs the separate pyscf-dispersion package; fail early with a clear message
try:
    from pyscf.dispersion import dftd3  # noqa: F401
except ImportError as err:
    raise SystemExit("D3 dispersion backend missing: pip install pyscf-dispersion") from err

# 1. Configure DFT functional, basis set, threads
# NOTE: in PySCF the plain libxc "wb97x-d" contains only the xc functional part. The
# "-d3bj" name adds an explicit D3(BJ) dispersion term (needs: pip install pyscf-dispersion).
xc_func = "wb97x-d3bj"      # range-separated hybrid + D3(BJ); alternatives: "wb97x-d", "pbe0"
basis_set = "def2-svp"      # quick test; switch to "def2-tzvp" for production-quality energies
density_fit = True          # RI-JK speeds up hybrid-DFT gradients a lot (error ~ microHartree)
nthreads = 8
lib.num_threads(nthreads)

# Optional: start from the PM7-optimized geometries to save expensive DFT optimization steps
mopac_geom_dir = Path("../mopac-ase/geometries")

out_dir = Path("geometries_pyscf")
out_dir.mkdir(exist_ok=True)

print("=== PySCF DFT Calculator Initialization ===")
print(f"Functional / Basis : {xc_func.upper()} / {basis_set}  (density fitting: {density_fit})")
print(f"OpenMP threads     : {nthreads}")
print("Reaction           : C4H6(g, 1,3-butadiene) + C2H4(g) -> C6H10(g, cyclohexene)")
print("Energy Definition  : H(298.15 K) = E_elec + ZPE + thermal (ASE IdealGasThermo)")


# 2. Minimal ASE calculator wrapping PySCF RKS (energy + analytic forces)
class PySCFCalculator(Calculator):
    implemented_properties = ["energy", "forces"]
    _announced = False

    def __init__(self, xc, basis, density_fit=True, grid_level=4, conv_tol=1e-10, **kwargs):
        super().__init__(**kwargs)
        self.xc = xc
        self.basis = basis
        self.density_fit = density_fit
        self.grid_level = grid_level
        self.conv_tol = conv_tol
        self._dm0 = None  # density from previous geometry as initial guess

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        atoms = self.atoms

        m = gto.M(
            atom=[(s, tuple(p)) for s, p in zip(atoms.get_chemical_symbols(), atoms.positions)],
            basis=self.basis,
            unit="Angstrom",
            charge=0,
            spin=0,
            verbose=0,
        )

        mf = dft.RKS(m, xc=self.xc)
        if self.density_fit:
            mf = mf.density_fit()
        mf.grids.level = self.grid_level
        mf.conv_tol = self.conv_tol
        mf.conv_tol_grad = 1e-7
        mf.max_cycle = 200

        dm0 = self._dm0 if (self._dm0 is not None and self._dm0.shape[-1] == m.nao) else None
        mf.kernel(dm0=dm0)

        # Dispersion is resolved when the SCF starts, so report it after the first kernel call
        if not PySCFCalculator._announced:
            print(f"  [PySCF] xc={mf.xc}, dispersion={getattr(mf, 'disp', None)}, "
                  f"DF={self.density_fit}, grid level={self.grid_level}")
            PySCFCalculator._announced = True

        if not mf.converged:
            raise RuntimeError("PySCF SCF did not converge")
        self._dm0 = mf.make_rdm1()

        self.results["energy"] = mf.e_tot * units.Hartree

        if "forces" in properties:
            g = mf.nuc_grad_method()
            g.grid_response = True        # consistent gradients for finite-difference Hessians
            grad = g.kernel()             # Hartree/Bohr
            self.results["forces"] = -np.asarray(grad) * (units.Hartree / units.Bohr)


def make_pyscf_calc():
    return PySCFCalculator(xc=xc_func, basis=basis_set, density_fit=density_fit)


# 3. Starting geometries
def build_species(name):
    if name == "C4H6":
        return molecule("butadiene")
    elif name == "C2H4":
        return molecule("C2H4")
    elif name == "C6H10":
        if "cyclohexene" in g2.names:
            return molecule("cyclohexene")
        # Half-chair cyclohexene (C6H10) starting coordinates (in Angstroms)
        symbols = ["C"] * 6 + ["H"] * 10
        positions = [
            [ 0.000,  0.668,  1.258],
            [ 0.000, -0.668,  1.258],
            [-0.145, -1.495, -0.002],
            [ 0.258, -0.748, -1.275],
            [-0.258,  0.748, -1.275],
            [ 0.145,  1.495, -0.002],
            [ 0.032,  1.225,  2.193],
            [ 0.032, -1.225,  2.193],
            [ 0.453, -2.411,  0.088],
            [-1.193, -1.825, -0.082],
            [ 1.351, -0.822, -1.352],
            [-0.145, -1.262, -2.155],
            [-1.351,  0.822, -1.352],
            [ 0.145,  1.262, -2.155],
            [ 1.193,  1.825, -0.082],
            [-0.453,  2.411,  0.088],
        ]
        return Atoms(symbols=symbols, positions=positions)
    raise ValueError(f"Unknown species: {name}")


def start_geometry(label):
    pm7_file = mopac_geom_dir / f"{label}_final.xyz"
    if pm7_file.exists():
        print(f"  ({label}: starting from PM7 geometry {pm7_file})")
        return read(str(pm7_file))
    return build_species(label)


# 4. Relax, write geometries, harmonic vibrations, thermochemistry
def relax_and_analyze(label, display_name, symmetrynumber, fmax=0.001):
    atoms = start_geometry(label)
    atoms.pbc = False
    atoms.calc = make_pyscf_calc()

    traj_file = out_dir / f"{label}_opt.traj"
    opt = BFGS(atoms, logfile=None, trajectory=str(traj_file))
    opt.run(fmax=fmax)
    e_opt = atoms.get_potential_energy()

    write(out_dir / f"{label}_opt_traj.xyz", read(str(traj_file), index=":"), format="extxyz")
    write(
        out_dir / f"{label}_final.xyz",
        atoms,
        format="xyz",
        comment=f"{display_name} {xc_func}/{basis_set} E = {e_opt:.6f} eV",
    )

    vib_dir = Path(f"vib_pyscf_{label}")
    if vib_dir.exists():
        shutil.rmtree(vib_dir)
    vib = Vibrations(atoms, name=str(vib_dir), delta=0.01)
    vib.run()
    vib_energies = vib.get_energies()
    vib.clean()
    if vib_dir.exists():
        shutil.rmtree(vib_dir)

    # Drop the 6 translational/rotational modes (all species are nonlinear)
    internal = vib_energies[6:]
    n_imag = int(sum(1 for e in internal if abs(np.imag(e)) > 0))
    internal_cm1 = [np.real(e) / invcm for e in internal]
    zpe = 0.5 * float(np.sum(np.real(internal)))

    thermo = IdealGasThermo(
        vib_energies=internal,
        potentialenergy=e_opt,
        atoms=atoms,
        geometry="nonlinear",
        symmetrynumber=symmetrynumber,
        spin=0,
    )
    h_298 = thermo.get_enthalpy(temperature=298.15, verbose=False)

    print(f"\n[{label:<5} ({display_name:<14})] steps={opt.nsteps:<3} | E_opt={e_opt:14.6f} eV | "
          f"ZPE={zpe:.4f} eV | imaginary modes={n_imag}")
    print(f"  Lowest internal modes (cm^-1): {', '.join(f'{f:.1f}' for f in internal_cm1[:4])}")
    return atoms, e_opt, zpe, h_298


# (label, display name, rotational symmetry number)
species_list = [
    ("C4H6",  "1,3-butadiene", 2),   # s-trans, C2h
    ("C2H4",  "ethene",        4),   # D2h
    ("C6H10", "cyclohexene",   2),   # half-chair, C2
]

print("\n=== Geometry Optimization & Harmonic Vibrations (delta = 0.01 A) ===")
results = {}
for label, disp, sigma in species_list:
    atoms, e_opt, zpe, h_298 = relax_and_analyze(label, disp, sigma, fmax=0.001)
    results[label] = {"atoms": atoms, "e_opt": e_opt, "zpe": zpe, "h_298": h_298}

# 5. Key double-bond lengths (atom pairs as in the MOPAC script)
d_cc_diene = results["C4H6"]["atoms"].get_distance(0, 1)
d_cc_ene   = results["C2H4"]["atoms"].get_distance(0, 1)
d_cc_cyhex = results["C6H10"]["atoms"].get_distance(0, 1)

print("\n=== Relaxed C=C Double Bond Lengths ===")
print(f"C=C in 1,3-butadiene : {d_cc_diene:.4f} A  (Expt: 1.341 A)")
print(f"C=C in ethene        : {d_cc_ene:.4f} A  (Expt: 1.339 A)")
print(f"C=C in cyclohexene   : {d_cc_cyhex:.4f} A  (Expt: 1.335 A)")

# 6. Reaction energetics: C4H6 + C2H4 -> C6H10
r = results
dE_elec = r["C6H10"]["e_opt"] - (r["C4H6"]["e_opt"] + r["C2H4"]["e_opt"])
dZPE = r["C6H10"]["zpe"] - (r["C4H6"]["zpe"] + r["C2H4"]["zpe"])
dE_0K = dE_elec + dZPE
dH_298 = r["C6H10"]["h_298"] - (r["C4H6"]["h_298"] + r["C2H4"]["h_298"])

print(f"\n=== Diels-Alder Reaction Enthalpy Summary (PySCF {xc_func.upper()}/{basis_set}) ===")
print(f"Electronic dE_rxn  : {dE_elec:8.4f} eV  ({dE_elec / (kJ / mol):8.2f} kJ/mol)")
print(f"Zero-Point dZPE    : {dZPE:+8.4f} eV  ({dZPE / (kJ / mol):+8.2f} kJ/mol)")
print(f"ZPE-Corr dE_0K     : {dE_0K:8.4f} eV  ({dE_0K / (kJ / mol):8.2f} kJ/mol)")
print(f"Enthalpy drH_298K  : {dH_298:8.4f} eV  ({dH_298 / (kJ / mol):8.2f} kJ/mol)  [Expt ref: -166.70 kJ/mol]")
