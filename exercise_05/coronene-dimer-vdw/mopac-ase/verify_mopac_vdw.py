from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from ase.units import kJ, mol
from ase.calculators.mopac import MOPAC

print("=================================================================")
print("  MOPAC Van der Waals (Dispersion) Verification: Coronene Dimer")
print("=================================================================")
print("Methods            : PM7 (built-in damped dispersion) & PM6-D3H4 (analytical D3)")
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

def make_mopac_calc(label, method):
    return MOPAC(
        label=label,
        method=method,
        task="1SCF GRADIENTS PRECISE",
        command="mopac PREFIX.mop > /dev/null 2>&1",
    )

def evaluate_mopac_vdw(method, tag, slip_x=1.42):
    print(f"\n--- MOPAC ({method}) ---")
    # 1. Relax isolated coronene monomer (36 atoms)
    monomer = build_coronene()
    calc_mono = make_mopac_calc(f"tmp_{tag}_mono", method)
    monomer.calc = calc_mono
    mono_traj = f"{tag}_monomer_opt.traj"
    opt = BFGS(monomer, trajectory=mono_traj, logfile=None)
    opt.run(fmax=0.005)

    # Ensure results['final_hof'] is populated and read Heat of Formation (in eV)
    _ = monomer.get_potential_energy()
    e_mono = calc_mono.get_final_heat_of_formation()
    write(f"{tag}_monomer_opt.xyz", monomer)
    print(
        f"Relaxed Monomer dHf: {e_mono / (kJ / mol):+.2f} kJ/mol "
        f"({len(monomer)} atoms, BFGS steps={opt.nsteps})"
    )

    # 2. Scan vertical interplanar distance z (slip-stacked by slip_x = 1.42 A)
    distances = [3.0, 3.2, 3.3, 3.4, 3.5, 3.6, 3.8, 4.0, 4.5, 5.0, 6.5]
    print(f"{'z (A)':>7} | {'E_int (eV)':>12} | {'E_int (kJ/mol)':>15}")
    print("-" * 41)

    scan_frames = []
    min_z, min_e_kj, best_dimer = None, 1e9, None

    for z in distances:
        top = monomer.copy()
        top.translate([slip_x, 0.0, z])
        dimer = monomer + top
        dimer.pbc = False

        calc_dimer = make_mopac_calc(f"tmp_{tag}_dimer", method)
        dimer.calc = calc_dimer
        _ = dimer.get_potential_energy()
        e_dimer = calc_dimer.get_final_heat_of_formation()

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
    scan_traj = f"{tag}_dimer_zscan.traj"
    eq_xyz = f"{tag}_dimer_eq_z{min_z:.2f}A.xyz"
    write(scan_traj, scan_frames)
    write(eq_xyz, best_dimer)

    print(f"=> Equilibrium minimum : z ≈ {min_z:.2f} A, E_vdW = {min_e_kj:+.2f} kJ/mol ({len(best_dimer)} atoms)")
    print(f"=> Saved trajectories  : {mono_traj}, {scan_traj}, {eq_xyz}")
    return min_z, min_e_kj

res_pm7 = evaluate_mopac_vdw("PM7", "mopac_pm7")
res_pm6 = evaluate_mopac_vdw("PM6-D3H4", "mopac_pm6_d3h4")

# Clean up temporary MOPAC scratch files
for ext in ("tmp_mopac_*.mop", "tmp_mopac_*.out", "tmp_mopac_*.arc"):
    for f in Path(".").glob(ext):
        f.unlink()

print("\n=== Van der Waals Benchmark Summary (Parallel-Displaced Coronene Dimer, 72 atoms) ===")
print(f"MOPAC (PM7)              : z_eq ≈ {res_pm7[0]:.2f} A | E_int = {res_pm7[1]:+7.2f} kJ/mol")
print(f"MOPAC (PM6-D3H4)         : z_eq ≈ {res_pm6[0]:.2f} A | E_int = {res_pm6[1]:+7.2f} kJ/mol")
print("Reference CCSD(T)/FN-DMC : z_eq ≈ 3.38 A | E_int ≈  -82.00 kJ/mol")

