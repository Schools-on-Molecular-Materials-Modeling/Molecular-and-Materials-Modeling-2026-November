import time
import numpy as np
import torch
from ase import Atoms, units
from ase.build import molecule
from ase.io import write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from ase.md.langevin import Langevin
from ase.geometry.analysis import Analysis
from mace.calculators import mace_off

print("=================================================================")
print("  ASE Molecular Dynamics Demo: Liquid Water (32 H2O, 300 K)")
print("=================================================================")

# 1. Build a 32-molecule water box at experimental liquid density (0.997 g/cm^3 -> L = 9.865 A)
n_mol = 32
box_len = 9.865  # (32 * 18.015 g/mol / (0.997 g/cm^3 * N_A))^(1/3) = 9.865 A
grid = [4, 4, 2] # 4 x 4 x 2 = 32 lattice sites

water_box = Atoms(cell=[box_len, box_len, box_len], pbc=True)
base_h2o = molecule("H2O")
rng = np.random.default_rng(42)

idx = 0
for ix in range(grid[0]):
    for iy in range(grid[1]):
        for iz in range(grid[2]):
            mol = base_h2o.copy()
            # Apply random rotation to break crystal symmetry
            angles = rng.uniform(0, 360, size=3)
            mol.euler_rotate(phi=angles[0], theta=angles[1], psi=angles[2], center="COM")
            pos = np.array([
                (ix + 0.5) * (box_len / grid[0]),
                (iy + 0.5) * (box_len / grid[1]),
                (iz + 0.5) * (box_len / grid[2]),
            ])
            mol.translate(pos - mol.get_center_of_mass())
            water_box += mol
            idx += 1

print(f"System             : {n_mol} H2O molecules ({len(water_box)} atoms)")
print(f"Cubic Box Length   : {box_len:.3f} A (density = 0.997 g/cm^3)")

# 2. Attach MACE-OFF23 (small) on GPU (float32 is fast and accurate for MD)
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Calculator         : MACE-OFF23 (small) on {device}")
water_box.calc = mace_off(model="small", device=device, default_dtype="float32")

# 3. Initialize velocities at 300 K and remove net momentum/angular momentum
T_target = 300.0  # K
MaxwellBoltzmannDistribution(water_box, temperature_K=T_target, rng=rng)
Stationary(water_box)
ZeroRotation(water_box)

# 4. Run NVT Langevin MD: 0.5 fs timestep (safe for O-H bond vibrations)
dt = 0.5 * units.fs
equil_steps = 1000   # 0.5 ps equilibration
prod_steps  = 4000   # 2.0 ps production
sample_every = 10    # save every 5 fs (400 frames)

dyn = Langevin(water_box, timestep=dt, temperature_K=T_target, friction=0.01 / units.fs)

t0 = time.perf_counter()
print(f"\nRunning {equil_steps} equilibration steps ({equil_steps * 0.5 / 1000:.2f} ps)...")
dyn.run(equil_steps)

traj_frames = []
temps = []

def record_frame():
    frame = water_box.copy()
    traj_frames.append(frame)
    temps.append(water_box.get_temperature())

dyn.attach(record_frame, interval=sample_every)
print(f"Running {prod_steps} production steps ({prod_steps * 0.5 / 1000:.2f} ps)...")
dyn.run(prod_steps)
elapsed = time.perf_counter() - t0

write("liquid_water_300K.traj", traj_frames)
print(f"Completed in {elapsed:.1f} s ({(equil_steps + prod_steps) / elapsed:.1f} steps/s)")
print(f"Average Production Temperature : {np.mean(temps):.1f} +/- {np.std(temps):.1f} K")

# 5. Compute Oxygen-Oxygen Radial Distribution Function g_OO(r)
# Note: rmax must be < box_len / 2 = 4.93 A
rmax = 4.80
nbins = 120
ana = Analysis(traj_frames)
rdf_list = ana.get_rdf(rmax=rmax, nbins=nbins, elements="O")
g_oo = np.mean(rdf_list, axis=0)
r_grid = np.linspace(rmax / (2 * nbins), rmax - rmax / (2 * nbins), nbins)

# Locate first solvation shell peak (2.4 A < r < 3.3 A)
mask_1st = (r_grid > 2.4) & (r_grid < 3.3)
peak_idx = np.argmax(g_oo[mask_1st])
r1_sim = r_grid[mask_1st][peak_idx]
g1_sim = g_oo[mask_1st][peak_idx]

# 6. Compute Self-Diffusion Coefficient D from Oxygen Mean-Squared Displacement (MSD)
# Unwrap periodic boundary crossings across production frames
o_indices = [i for i, sym in enumerate(water_box.get_chemical_symbols()) if sym == "O"]
pos_series = np.array([f.get_positions()[o_indices] for f in traj_frames])  # (n_frames, 32, 3)

# Unwrap jumps > box_len / 2
diffs = np.diff(pos_series, axis=0)
diffs -= box_len * np.round(diffs / box_len)
unwrapped = np.concatenate([pos_series[0:1], pos_series[0:1] + np.cumsum(diffs, axis=0)], axis=0)

# MSD(t) relative to t=0 of production window
msd = np.mean(np.sum((unwrapped - unwrapped[0]) ** 2, axis=-1), axis=1)  # in A^2
times_ps = np.arange(len(traj_frames)) * (sample_every * 0.5 * 1e-3)     # in ps

# Linear fit over ballistic-free regime (t > 0.5 ps): MSD(t) = 6 * D * t
fit_mask = times_ps >= 0.5
slope_A2_ps, _ = np.polyfit(times_ps[fit_mask], msd[fit_mask], 1)
# Convert A^2/ps -> 10^-5 cm^2/s (1 A^2/ps = 10^-16 cm^2 / 10^-12 s = 10^-4 cm^2/s = 10 * 10^-5 cm^2/s)
D_sim = (slope_A2_ps / 6.0) * 10.0

print("\n=== Liquid Water MD Observables vs. Experimental Reference (298-300 K) ===")
print(f"1st O-O Peak Position r1 : {r1_sim:.2f} A   (Expt X-ray/Neutron: 2.80 A)")
print(f"1st O-O Peak Height g(r1): {g1_sim:.2f}     (Expt X-ray/Neutron: 2.57 - 2.80)")
print(f"Self-Diffusion Coeff D   : {D_sim:.2f} x 10^-5 cm^2/s  (Expt NMR: 2.30 x 10^-5 cm^2/s)")
print("Saved trajectory         : liquid_water_300K.traj (view with: ase gui liquid_water_300K.traj)")

