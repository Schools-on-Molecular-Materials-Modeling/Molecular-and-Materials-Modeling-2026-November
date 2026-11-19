import time
import numpy as np
import torch
from pathlib import Path
from tqdm import tqdm
from ase import Atoms, units
from ase.build import molecule
from ase.constraints import FixCom
from ase.io import read, write
from ase.md.velocitydistribution import thermalize_momenta, Stationary
from ase.md.langevin import Langevin
from ase.geometry.rdf import get_rdf
from mace.calculators import mace_off

print("=================================================================")
print("  ASE Molecular Dynamics Demo: Liquid Water (32 H2O, 300 K)")
print("=================================================================")

n_mol = 32
box_len = 9.865  # experimental liquid water density 0.997 g/cm^3 -> L = 9.865 A
sample_every = 10
dt_fs = 0.5
traj_file = Path("liquid_water_300K.traj")

# Set rerun_md = True to run a fresh MD trajectory with tqdm progress bars;
# set rerun_md = False if liquid_water_300K.traj is already on disk and you only want RDF/MSD analysis.
rerun_md = True

if rerun_md or not traj_file.is_file():
    # 1. Build a 32-molecule water box (4 x 4 x 2 grid with randomized orientations)
    grid = (4, 4, 2)
    water_box = Atoms(cell=[box_len, box_len, box_len], pbc=True)
    base_h2o = molecule("H2O")
    rng = np.random.default_rng(42)

    for ix in range(grid[0]):
        for iy in range(grid[1]):
            for iz in range(grid[2]):
                mol = base_h2o.copy()
                angles = rng.uniform(0, 360, size=3)
                mol.euler_rotate(phi=angles[0], theta=angles[1], psi=angles[2], center="COM")
                pos = np.array([
                    (ix + 0.5) * (box_len / grid[0]),
                    (iy + 0.5) * (box_len / grid[1]),
                    (iz + 0.5) * (box_len / grid[2]),
                ])
                mol.translate(pos - mol.get_center_of_mass())
                water_box += mol

    print(f"System             : {n_mol} H2O molecules ({len(water_box)} atoms)")
    print(f"Cubic Box Length   : {box_len:.3f} A (density = 0.997 g/cm^3)")

    # 2. Attach MACE-OFF23 (small) on GPU
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Calculator         : MACE-OFF23 (small) on {device}")
    water_box.calc = mace_off(model="small", device=device, default_dtype="float32")

    # 3. Initialize velocities at 300 K and fix center-of-mass drift
    T_target = 300.0
    thermalize_momenta(water_box, temperature_K=T_target, rng=rng)
    Stationary(water_box)
    water_box.set_constraint(FixCom())

    # 4. Run NVT Langevin MD with live tqdm progress bars
    dt = dt_fs * units.fs
    equil_steps = 1000   # 0.50 ps equilibration
    prod_steps = 4000    # 2.00 ps production
    pbar_interval = 10

    dyn = Langevin(
        water_box,
        timestep=dt,
        temperature_K=T_target,
        friction=0.01 / units.fs,
        fixcm=False,
    )

    t0 = time.perf_counter()

    print(f"\nEquilibration ({equil_steps * dt_fs / 1000:.2f} ps):")
    pbar_eq = tqdm(total=equil_steps, unit="step", ncols=85)

    def update_eq_pbar():
        pbar_eq.update(pbar_interval)
        pbar_eq.set_postfix(T=f"{water_box.get_temperature():.1f}K")

    dyn.attach(update_eq_pbar, interval=pbar_interval)
    dyn.run(equil_steps)
    pbar_eq.close()
    dyn.observers.clear()

    traj_frames = []
    temps = []

    def record_frame():
        traj_frames.append(water_box.copy())
        temps.append(water_box.get_temperature())

    print(f"\nProduction ({prod_steps * dt_fs / 1000:.2f} ps):")
    pbar_prod = tqdm(total=prod_steps, unit="step", ncols=85)

    def update_prod_pbar():
        pbar_prod.update(pbar_interval)
        pbar_prod.set_postfix(T=f"{water_box.get_temperature():.1f}K")

    dyn.attach(record_frame, interval=sample_every)
    dyn.attach(update_prod_pbar, interval=pbar_interval)
    dyn.run(prod_steps)
    pbar_prod.close()

    elapsed = time.perf_counter() - t0
    write(traj_file, traj_frames)
    print(f"\nCompleted in {elapsed:.1f} s ({(equil_steps + prod_steps) / elapsed:.1f} steps/s)")
    print(f"Average Production Temperature : {np.mean(temps):.1f} +/- {np.std(temps):.1f} K")
else:
    print(f"Loading existing trajectory: {traj_file}")
    traj_frames = read(traj_file, index=":")
    print(f"Loaded {len(traj_frames)} frames ({len(traj_frames[0])} atoms)")

# 5. Compute Oxygen-Oxygen Radial Distribution Function g_OO(r) using ase.geometry.rdf.get_rdf
rmax = 4.80
nbins = 120
rdf_per_frame = []
r_grid = None
for frame in traj_frames:
    g_r, r_grid = get_rdf(frame, rmax=rmax, nbins=nbins, elements=(8, 8))
    rdf_per_frame.append(g_r)

g_oo = np.mean(rdf_per_frame, axis=0)

# Locate first solvation shell peak (2.4 A < r < 3.3 A)
mask_1st = (r_grid > 2.4) & (r_grid < 3.3)
peak_idx = np.argmax(g_oo[mask_1st])
r1_sim = r_grid[mask_1st][peak_idx]
g1_sim = g_oo[mask_1st][peak_idx]

# 6. Compute Self-Diffusion Coefficient D from Oxygen Mean-Squared Displacement (MSD)
o_indices = [i for i, sym in enumerate(traj_frames[0].get_chemical_symbols()) if sym == "O"]
pos_series = np.array([f.get_positions()[o_indices] for f in traj_frames])  # (n_frames, 32, 3)

diffs = np.diff(pos_series, axis=0)
diffs -= box_len * np.round(diffs / box_len)
unwrapped = np.concatenate([pos_series[0:1], pos_series[0:1] + np.cumsum(diffs, axis=0)], axis=0)

msd = np.mean(np.sum((unwrapped - unwrapped[0]) ** 2, axis=-1), axis=1)  # in A^2
times_ps = np.arange(len(traj_frames)) * (sample_every * dt_fs * 1e-3)   # in ps

# Linear fit over ballistic-free regime (t >= 0.5 ps): MSD(t) = 6 * D * t
fit_mask = times_ps >= 0.5
slope_A2_ps, _ = np.polyfit(times_ps[fit_mask], msd[fit_mask], 1)
D_sim = (slope_A2_ps / 6.0) * 10.0  # in 10^-5 cm^2/s

print("\n=== Liquid Water MD Observables vs. Experimental Reference (298-300 K) ===")
print(f"1st O-O Peak Position r1 : {r1_sim:.2f} A   (Expt X-ray/Neutron: 2.80 A)")
print(f"1st O-O Peak Height g(r1): {g1_sim:.2f}     (Expt X-ray/Neutron: 2.57 - 2.80)")
print(f"Self-Diffusion Coeff D   : {D_sim:.2f} x 10^-5 cm^2/s  (Expt NMR: 2.30 x 10^-5 cm^2/s)")
print(f"Saved trajectory         : {traj_file} (view with: ase gui {traj_file})")

