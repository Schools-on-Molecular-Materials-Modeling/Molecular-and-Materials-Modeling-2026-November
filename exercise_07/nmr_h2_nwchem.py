import re
import subprocess
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem

# 1. Configure DFT functional, triple-zeta basis set, and MPI ranks
xc_func = "b3lyp"           # Try "b3lyp", "pbe0", or "kt2" (specialized for NMR shifts)
basis_set = "6-311++G**"    # Triple-zeta + polarization + diffuse functions for GIAO NMR
nprocs = 6
nwchem_cmd = f"mpirun -np {nprocs} nwchem PREFIX.nwi > PREFIX.nwo"

print("=== NWChem GIAO NMR Chemical Shift Calculation ===")
print(f"Functional / Basis : {xc_func.upper()} / {basis_set}")
print(f"MPI Launch Command : {nwchem_cmd}")
print("Target Molecule    : H2 (simplest molecule, 2 atoms)")
print("NMR Reference      : TMS = Si(CH3)4 (IUPAC 1H standard, 29 atoms)")

def build_tms():
    """Construct Td tetramethylsilane Si(CH3)4 (29 atoms) in Angstroms."""
    symbols = ["Si"]
    positions = [[0.0, 0.0, 0.0]]

    # Tetrahedral Si-C vectors (r_SiC = 1.875 A -> 1.875 / sqrt(3) = 1.0825 A)
    d = 1.875 / np.sqrt(3.0)
    c_dirs = np.array([
        [ 1.0,  1.0,  1.0],
        [ 1.0, -1.0, -1.0],
        [-1.0,  1.0, -1.0],
        [-1.0, -1.0,  1.0],
    ])
    c_coords = d * c_dirs

    # Staggered methyl hydrogens (r_CH = 1.093 A, tetrahedral angle 109.47 deg)
    r_ch = 1.093
    cos_t = np.cos(np.deg2rad(109.4712))  # -1/3 along C->Si axis
    sin_t = np.sin(np.deg2rad(109.4712))  # 2*sqrt(2)/3 perpendicular

    for c_pos in c_coords:
        symbols.append("C")
        positions.append(c_pos.tolist())

        # Unit vector pointing from C away from Si
        u_axis = c_pos / np.linalg.norm(c_pos)
        # Orthonormal basis perpendicular to u_axis
        ref = np.array([1.0, 0.0, 0.0]) if abs(u_axis) < 0.9 else np.array([0.0, 1.0, 0.0])
        v1 = np.cross(u_axis, ref)
        v1 /= np.linalg.norm(v1)
        v2 = np.cross(u_axis, v1)

        for k in range(3):
            phi = np.deg2rad(120.0 * k)
            h_dir = (-cos_t) * u_axis + sin_t * (np.cos(phi) * v1 + np.sin(phi) * v2)
            symbols.append("H")
            positions.append((c_pos + r_ch * h_dir).tolist())

    return Atoms(symbols=symbols, positions=positions, pbc=False)

def make_nwchem_calc(label, run_nmr=False):
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
    if run_nmr:
        params["property"] = {"shielding": None}
        params["task"] = "property"
    return NWChem(**params)

def optimize_and_compute_nmr(name, atoms, fmax=0.001):
    # Step 1: Geometry optimization with ASE BFGS
    opt_label = f"nwchem_opt_{name}"
    atoms.calc = make_nwchem_calc(opt_label, run_nmr=False)
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=fmax)
    e_opt = atoms.get_potential_energy()

    # Step 2: GIAO NMR shielding tensor calculation via NWChem property module
    nmr_label = f"nwchem_nmr_{name}"
    nmr_calc = make_nwchem_calc(nmr_label, run_nmr=True)
    nmr_calc.write_input(atoms)

    cmd = nwchem_cmd.replace("PREFIX", nmr_label)
    subprocess.run(cmd, shell=True, check=True)

    nwo_text = Path(f"{nmr_label}.nwo").read_text(errors="ignore")
    iso_values = [float(val) for val in re.findall(r"isotropic\s*=\s*([+-]?\d+\.\d+)", nwo_text)]

    symbols = atoms.get_chemical_symbols()
    if len(iso_values) != len(symbols):
        raise RuntimeError(f"Expected {len(symbols)} shielding values for {name}, parsed {len(iso_values)}")

    h_shieldings = [iso for sym, iso in zip(symbols, iso_values) if sym == "H"]
    sigma_h_mean = float(np.mean(h_shieldings))
    sigma_h_std = float(np.std(h_shieldings))

    print(
        f"[{name:<4}] steps={opt.nsteps:<2} | E_opt={e_opt:12.6f} eV | "
        f"sigma_iso(1H) = {sigma_h_mean:7.3f} +/- {sigma_h_std:.3f} ppm"
    )
    return atoms, sigma_h_mean, iso_values

print("\n=== Step 1 & 2: Geometry Optimization + GIAO NMR Shielding ===")
h2_atoms = molecule("H2")
h2_atoms.pbc = False
h2_opt, sigma_h2, _ = optimize_and_compute_nmr("H2", h2_atoms, fmax=0.001)

tms_atoms = build_tms()
tms_opt, sigma_tms, tms_all = optimize_and_compute_nmr("TMS", tms_atoms, fmax=0.002)

# Calculate 1H chemical shift relative to TMS: delta = sigma_ref - sigma_sample
delta_h2_ppm = sigma_tms - sigma_h2
d_hh = h2_opt.get_distance(0, 1)

print("\n=== 1H NMR Chemical Shift Summary for H2 (vs. TMS) ===")
print(f"Relaxed H-H Bond Length   : {d_hh:.4f} A  (Expt: 0.741 A)")
print(f"sigma_iso(1H) in TMS      : {sigma_tms:7.3f} ppm  (Expt gas-phase: ~30.58 ppm)")
print(f"sigma_iso(1H) in H2       : {sigma_h2:7.3f} ppm  (Expt gas-phase:  26.26 ppm)")
print(f"Chemical Shift delta(1H)  : {delta_h2_ppm:+7.3f} ppm")
print("Reference Values          : +4.32 ppm (gas phase, 300 K) | +4.47 to +4.62 ppm (solution)")

# Clean up NWChem scratch files
for pattern in ("nwchem_opt_*", "nwchem_nmr_*"):
    for f in Path(".").glob(pattern):
        if f.is_file():
            f.unlink()

