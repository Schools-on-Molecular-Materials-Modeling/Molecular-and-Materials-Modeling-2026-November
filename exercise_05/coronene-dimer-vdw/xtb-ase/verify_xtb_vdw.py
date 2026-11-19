import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from ase.units import kJ, mol
from xtb.ase.calculator import XTB

print("=================================================================")
print("  GFN2-xTB Van der Waals (D4 Dispersion) Verification: Coronene Dimer")
print("=================================================================")
print("Method             : GFN2-xTB (Extended Tight-Binding + Analytical D4)")
print("System             : Coronene Dimer (C24H12 ... C24H12, 72 atoms)")
print("Benchmark Ref      : CCSD(T)/CBS & FN-DMC (~ -81.6 to -85.8 kJ/mol slip-stacked)")

def build_coronene():
    """Build planar D6h coronene (C24H12, 36 atoms) in the xy-plane."""
    r_c_inner = 1.425
    r_c_spoke = 2.845
    r_c_rim   = 3.685
    r_h       = 4.765

    symbols = []
    positions = []

    # 1. Inner benzene core (6 C atoms)
    for k in range(6):
        ang = np.deg2rad(60 * k)
        symbols.append("C")
        positions.append([r_c_inner * np.cos(ang), r_c_inner * np.sin(ang), 0.0])

    # 2. Spoke-junction carbons (6 C atoms)
    for k in range(6):
        ang = np.deg2rad(60 * k)
        symbols.append("C")
        positions.append([r_c_spoke * np.cos(ang), r_c_spoke * np.sin(ang), 0.0])

    # 3. Outer CH carbons (12 C atoms, +/- 19.1 deg around each spoke)
    for k in range(6):
        for offset in (-19.1, +19.1):
            ang = np.deg2rad(60 * k + offset)
            symbols.append("C")
            positions.append([r_c_rim * np.cos(ang), r_c_rim * np.sin(ang), 0.0])

    # 4. Outer hydrogens (12 H atoms)
    for k in range(6):
        for offset in (-19.1, +19.1):
            ang = np.deg2rad(60 * k + offset)
            symbols.append("H")
            positions.append([r_h * np.cos(ang), r_h * np.sin(ang), 0.0])

    atoms = Atoms(symbols=symbols, positions=positions)
    atoms.pbc = False
    return atoms

calc = XTB(method="GFN2-xTB")

# 1. Relax isolated coronene monomer (36 atoms)
monomer = build_coronene()
monomer.calc = calc
opt = BFGS(monomer, trajectory="gfn2_xtb_monomer_opt.traj", logfile=None)
opt.run(fmax=0.005)
e_mono = monomer.get_potential_energy()
write("gfn2_xtb_monomer_opt.xyz", monomer)
print(f"\nRelaxed Monomer E  : {e_mono:.6f} eV ({len(monomer)} atoms, BFGS steps={opt.nsteps})")

# 2. Scan vertical interplanar distance z (slip-stacked by slip_x = 1.42 A, AB-graphite offset)
slip_x = 1.42
distances = [3.0, 3.2, 3.3, 3.4, 3.5, 3.6, 3.8, 4.0, 4.5, 5.0, 6.5]
print(f"\n{'z (A)':>7} | {'E_int (eV)':>12} | {'E_int (kJ/mol)':>15}")
print("-" * 41)

scan_frames = []
min_z, min_e_kj, best_dimer = None, 1e9, None

for z in distances:
    top = monomer.copy()
    top.translate([slip_x, 0.0, z])
    dimer = monomer + top
    dimer.pbc = False
    dimer.calc = calc

    e_dimer = dimer.get_potential_energy()
    e_int_ev = e_dimer - 2.0 * e_mono
    e_int_kj = e_int_ev / (kJ / mol)
    print(f"{z:7.2f} | {e_int_ev:+12.4f} | {e_int_kj:+15.2f}")

    frame = dimer.copy()
    frame.info["z_angstrom"] = z
    frame.info["E_int_eV"] = e_int_ev
    frame.info["E_int_kJ_mol"] = e_int_kj
    scan_frames.append(frame)

    if e_int_kj < min_e_kj:
        min_e_kj = e_int_kj
        min_z = z
        best_dimer = frame.copy()

# 3. Save full z-scan trajectory and equilibrium dimer geometry
scan_traj = "gfn2_xtb_dimer_zscan.traj"
eq_xyz = f"gfn2_xtb_dimer_eq_z{min_z:.2f}A.xyz"
write(scan_traj, scan_frames)
write(eq_xyz, best_dimer)

print(f"\n=> Equilibrium minimum : z ≈ {min_z:.2f} A, E_vdW = {min_e_kj:+.2f} kJ/mol ({len(best_dimer)} atoms)")
print(f"=> Saved trajectories  : gfn2_xtb_monomer_opt.traj, {scan_traj}, {eq_xyz}")
print("=> Reference L7/FN-DMC : z_eq ≈ 3.38 A, E_vdW ≈ -81.6 to -85.8 kJ/mol")

