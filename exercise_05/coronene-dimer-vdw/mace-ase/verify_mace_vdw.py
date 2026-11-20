import numpy as np
import torch
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from ase.units import kJ, mol
from mace.calculators import mace_off, mace_mp

device = "cuda" if torch.cuda.is_available() else "cpu"
print("=================================================================")
print("  MACE Van der Waals (Dispersion) Verification: Coronene Dimer")
print("=================================================================")
print(f"Compute Device     : {device} (PyTorch {torch.__version__})")
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
    atoms.center(vacuum=15.0)
    return atoms

def evaluate_vdw_curve(label, tag, calc, slip_x=1.42):
    print(f"\n--- {label} ---")
    # Relax isolated monomer and save trajectory + final geometry
    monomer = build_coronene()
    monomer.calc = calc
    mono_traj = f"{tag}_monomer_opt.traj"
    opt = BFGS(monomer, trajectory=mono_traj, logfile=None)
    opt.run(fmax=0.005)
    e_mono = monomer.get_potential_energy()
    write(f"{tag}_monomer_opt.xyz", monomer)
    print(f"Relaxed Monomer E  : {e_mono:.6f} eV ({len(monomer)} atoms, BFGS steps={opt.nsteps})")

    # Scan vertical interplanar distance z (slip-stacked by slip_x = 1.42 A, AB-graphite offset)
    distances = [3.0, 3.2, 3.3, 3.4, 3.45, 3.5, 3.55, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
    print(f"{'z (A)':>7} | {'E_int (eV)':>12} | {'E_int (kJ/mol)':>15}")
    print("-" * 41)

    scan_frames = []
    min_z, min_e_kj, best_dimer = None, 1e9, None

    for z in distances:
        top = monomer.copy()
        top.translate([slip_x, 0.0, z])
        dimer = monomer + top
        dimer.center(vacuum=15.0)
        dimer.calc = calc

        e_dimer = dimer.get_potential_energy()
        e_int_ev = e_dimer - 2.0 * e_mono
        e_int_kj = e_int_ev / (kJ / mol)
        print(f"{z:7.2f} | {e_int_ev:+12.4f} | {e_int_kj:+15.2f}")

        # Attach metadata and store frame for trajectory export
        frame = dimer.copy()
        frame.info["z_angstrom"] = z
        frame.info["E_int_eV"] = e_int_ev
        frame.info["E_int_kJ_mol"] = e_int_kj
        scan_frames.append(frame)

        if e_int_kj < min_e_kj:
            min_e_kj = e_int_kj
            min_z = z
            best_dimer = frame.copy()

    # Save full z-scan trajectory and equilibrium dimer geometry
    scan_traj = f"{tag}_dimer_zscan.traj"
    eq_xyz = f"{tag}_dimer_eq_z{min_z:.2f}A.xyz"
    write(scan_traj, scan_frames)
    write(eq_xyz, best_dimer)

    print(f"=> Equilibrium minimum : z ≈ {min_z:.2f} A, E_vdW = {min_e_kj:+.2f} kJ/mol ({len(best_dimer)} atoms)")
    print(f"=> Saved trajectories  : {mono_traj}, {scan_traj}, {eq_xyz}")
    return min_z, min_e_kj

# 1. MACE-OFF23 (trained on wB97M-D3(BJ) -> includes dispersion natively)
calc_off = mace_off(model="medium", device=device, default_dtype="float64")
res_off = evaluate_vdw_curve("MACE-OFF23 (medium, native wB97M-D3BJ)", "mace_off23", calc_off)

# 2. MACE-MP-0 WITHOUT D3 dispersion (pure PBE -> misses long-range London dispersion)
calc_mp_nod3 = mace_mp(model="medium", dispersion=False, device=device, default_dtype="float64")
res_mp_nod3 = evaluate_vdw_curve("MACE-MP-0 (medium, dispersion=False / bare PBE)", "mace_mp0_nod3", calc_mp_nod3)

# 3. MACE-MP-0 WITH D3 dispersion (PBE-D3 -> restores vdW attraction)
calc_mp_d3 = mace_mp(model="medium", dispersion=True, device=device, default_dtype="float64")
res_mp_d3 = evaluate_vdw_curve("MACE-MP-0 (medium, dispersion=True / PBE+D3)", "mace_mp0_d3", calc_mp_d3)

print("\n=== Van der Waals Benchmark Summary (Parallel-Displaced Coronene Dimer, 72 atoms) ===")
print(f"MACE-OFF23 (native D3BJ) : z_eq ≈ {res_off[0]:.2f} A | E_int = {res_off[1]:+7.2f} kJ/mol")
print(f"MACE-MP-0 (no dispersion): z_eq ≈ {res_mp_nod3[0]:.2f} A | E_int = {res_mp_nod3[1]:+7.2f} kJ/mol")
print(f"MACE-MP-0 (+ D3)         : z_eq ≈ {res_mp_d3[0]:.2f} A | E_int = {res_mp_d3[1]:+7.2f} kJ/mol")
print("Reference CCSD(T)/FN-DMC : z_eq ≈ 3.38 A | E_int ≈  -82.00 kJ/mol")

