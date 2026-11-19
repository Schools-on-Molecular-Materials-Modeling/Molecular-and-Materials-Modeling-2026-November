import numpy as np
from ase.build import molecule
from ase.optimize import BFGS
from pyscf import gto, dft
from pyscf.ase.ase import PySCF
from pyscf.prop.magnetizability import rks

# 1. Configure DFT functional and augmented triple-zeta basis set
xc_func = "b3lyp"
basis_set = "aug-cc-pVTZ"

# Exact CODATA unit conversion factors from atomic units (e^2 a0^2 / m_e)
AU_TO_J_T2 = 7.8910366e-29      # 1 a.u. = 78.910366 x 10^-30 J/T^2
AU_TO_CM3_MOL = 4.751991e-6     # 1 a.u. = 4.751991 x 10^-6 cm^3/mol (cgs-emu molar susceptibility)

print("=== ASE + PySCF GIAO Molecular Magnetizability Calculation ===")
print(f"Functional / Basis : {xc_func.upper()} / {basis_set}")
print("Target Molecule    : H2 (simplest closed-shell diamagnetic molecule)")

# 2. Build H2 in ASE and optimize geometry with ASE BFGS + PySCF calculator
h2_atoms = molecule("H2")
h2_atoms.pbc = False

mol = gto.Mole()
mol.atom = [(sym, pos) for sym, pos in zip(h2_atoms.get_chemical_symbols(), h2_atoms.get_positions())]
mol.basis = basis_set
mol.charge = 0
mol.spin = 0
mol.verbose = 0
mol.build()

mf = dft.RKS(mol)
mf.xc = xc_func
mf.grids.level = 4

h2_atoms.calc = PySCF(method=mf)
opt = BFGS(h2_atoms, logfile=None)
opt.run(fmax=0.001)

e_opt = h2_atoms.get_potential_energy()
d_hh = h2_atoms.get_distance(0, 1)

# 3. Rebuild PySCF Mole at the ASE-relaxed equilibrium geometry and compute GIAO magnetizability
mol_opt = gto.Mole()
mol_opt.atom = [(sym, pos) for sym, pos in zip(h2_atoms.get_chemical_symbols(), h2_atoms.get_positions())]
mol_opt.basis = basis_set
mol_opt.charge = 0
mol_opt.spin = 0
mol_opt.verbose = 0
mol_opt.build()

mf_opt = dft.RKS(mol_opt)
mf_opt.xc = xc_func
mf_opt.grids.level = 4
mf_opt.kernel()

# GIAO magnetizability tensor (3x3 matrix in atomic units)
chi_tensor_au = rks.Magnetizability(mf_opt).kernel()
chi_iso_au = float(np.trace(chi_tensor_au) / 3.0)

# Principal components (sorted eigenvalues: two degenerate perpendicular, one parallel)
eigs = np.sort(np.linalg.eigvalsh(chi_tensor_au))
if abs(eigs[0] - eigs[1]) < abs(eigs[1] - eigs[2]):
    chi_perp_au = 0.5 * (eigs[0] + eigs[1])
    chi_par_au = eigs[2]
else:
    chi_perp_au = 0.5 * (eigs[1] + eigs[2])
    chi_par_au = eigs[0]
aniso_au = float(chi_par_au - chi_perp_au)

chi_iso_si = chi_iso_au * AU_TO_J_T2 * 1e30        # in 10^-30 J/T^2
chi_iso_cgs = chi_iso_au * AU_TO_CM3_MOL * 1e6     # in 10^-6 cm^3/mol

print(f"\n[H2] BFGS steps={opt.nsteps} | E_opt={e_opt:.6f} eV | r(H-H)={d_hh:.4f} A")
print("\n=== GIAO Molecular Magnetizability Summary for H2 ===")
print(f"Isotropic chi_iso (a.u.)     : {chi_iso_au:+8.4f} a.u.             (Expt: -0.847 a.u.)")
print(f"Isotropic chi_iso (SI)       : {chi_iso_si:+8.3f} x 10^-30 J/T^2   (Expt: -6.68 x 10^-30 J/T^2)")
print(f"Molar Susceptibility chi_M   : {chi_iso_cgs:+8.3f} x 10^-6 cm^3/mol (Expt: -4.02 x 10^-6 cm^3/mol)")
print(f"Magnetizability Anisotropy   : {aniso_au:+8.4f} a.u.             (Expt: +0.114 a.u.)")

