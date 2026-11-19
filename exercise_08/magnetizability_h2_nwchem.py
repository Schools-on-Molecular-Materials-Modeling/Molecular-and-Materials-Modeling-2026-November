import re
import subprocess
from pathlib import Path
import numpy as np
from ase.build import molecule
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem

# 1. Configure DFT functional, augmented triple-zeta basis set, and MPI ranks
xc_func = "b3lyp"
basis_set = "aug-cc-pVTZ"   # Augmented basis functions are essential for magnetic response properties
nprocs = 6
nwchem_cmd = f"mpirun -np {nprocs} nwchem PREFIX.nwi > PREFIX.nwo"

# Exact CODATA unit conversion factors from atomic units (e^2 a0^2 / m_e)
AU_TO_J_T2 = 7.8910366e-29      # 1 a.u. = 78.910366 x 10^-30 J/T^2
AU_TO_CM3_MOL = 4.751991e-6     # 1 a.u. = 4.751991 x 10^-6 cm^3/mol (cgs-emu molar susceptibility)

print("=== NWChem GIAO Molecular Magnetizability Calculation ===")
print(f"Functional / Basis : {xc_func.upper()} / {basis_set}")
print(f"MPI Launch Command : {nwchem_cmd}")
print("Target Molecule    : H2 (simplest closed-shell diamagnetic molecule)")

def make_nwchem_calc(label, run_prop=False):
    params = {
        "label": label,
        "theory": "dft",
        "basis": basis_set,
        "dft": {
            "xc": xc_func,
            "mult": 1,
            "grid": "fine",
            "convergence": {
                "energy": 1e-8,
                "density": 1e-6,
                "gradient": 5e-6,
            },
        },
        "command": nwchem_cmd,
    }
    if run_prop:
        # NWChem property module keyword for molecular magnetizability (susceptibility) tensor
        params["property"] = {"susceptibility": None}
        params["task"] = "property"
    return NWChem(**params)

# 2. Optimize H2 geometry with ASE BFGS
h2_atoms = molecule("H2")
h2_atoms.pbc = False
h2_atoms.calc = make_nwchem_calc("nwchem_opt_H2", run_prop=False)
opt = BFGS(h2_atoms, logfile=None)
opt.run(fmax=0.001)
e_opt = h2_atoms.get_potential_energy()
d_hh = h2_atoms.get_distance(0, 1)

# 3. Compute GIAO molecular magnetizability tensor via NWChem property module
prop_label = "nwchem_mag_H2"
prop_calc = make_nwchem_calc(prop_label, run_prop=True)
prop_calc.write_input(h2_atoms)
subprocess.run(nwchem_cmd.replace("PREFIX", prop_label), shell=True, check=True)

nwo_text = Path(f"{prop_label}.nwo").read_text(errors="ignore")

# Parse isotropic magnetizability (in a.u.) and principal tensor components from NWChem output
iso_match = re.search(r"isotropic\s*=\s*([+-]?\d+\.\d+)", nwo_text, re.IGNORECASE)
if not iso_match:
    raise RuntimeError(f"Could not parse isotropic magnetizability from {prop_label}.nwo")

chi_iso_au = float(iso_match.group(1))
chi_iso_si = chi_iso_au * AU_TO_J_T2 * 1e30        # in 10^-30 J/T^2
chi_iso_cgs = chi_iso_au * AU_TO_CM3_MOL * 1e6     # in 10^-6 cm^3/mol

# Parse principal components (eigenvalues of the 3x3 susceptibility tensor) if printed
eig_match = re.search(
    r"principal\s+components.*?\n\s*([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)\s+([+-]?\d+\.\d+)",
    nwo_text,
    re.IGNORECASE | re.DOTALL,
)

print(f"\n[H2] BFGS steps={opt.nsteps} | E_opt={e_opt:.6f} eV | r(H-H)={d_hh:.4f} A")
print("\n=== Molecular Magnetizability Summary for H2 ===")
print(f"Isotropic chi_iso (a.u.)     : {chi_iso_au:+8.4f} a.u.          (Expt: -0.847 a.u.)")
print(f"Isotropic chi_iso (SI)       : {chi_iso_si:+8.3f} x 10^-30 J/T^2 (Expt: -6.68 x 10^-30 J/T^2)")
print(f"Molar Susceptibility chi_M   : {chi_iso_cgs:+8.3f} x 10^-6 cm^3/mol (Expt: -4.02 x 10^-6 cm^3/mol)")

if eig_match:
    eigs = sorted([float(eig_match.group(i)) for i in (1, 2, 3)])
    # For linear H2 aligned along one axis, two components are perpendicular (degenerate) and one is parallel
    chi_perp = 0.5 * (eigs[0] + eigs[1]) if abs(eigs[0] - eigs[1]) < abs(eigs[1] - eigs[2]) else 0.5 * (eigs[1] + eigs[2])
    chi_par = eigs[2] if abs(eigs[0] - eigs[1]) < abs(eigs[1] - eigs[2]) else eigs[0]
    aniso_au = chi_par - chi_perp
    print(f"Magnetizability Anisotropy   : {aniso_au:+8.4f} a.u.          (Expt: +0.114 a.u.)")

# Clean up NWChem scratch files
for pattern in ("nwchem_opt_*", "nwchem_mag_*"):
    for f in Path(".").glob(pattern):
        if f.is_file():
            f.unlink()

