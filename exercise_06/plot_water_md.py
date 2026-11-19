import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from ase.io import read
from ase.geometry.rdf import get_rdf

# 1. Load saved MD trajectory
traj_file = Path("liquid_water_300K.traj")
if not traj_file.is_file():
    raise FileNotFoundError(f"Trajectory file not found: {traj_file}")

traj_frames = read(traj_file, index=":")
box_len = traj_frames[0].cell.lengths()[0]
dt_frame_ps = 0.005  # 10 steps * 0.5 fs = 5 fs = 0.005 ps per saved frame

print(f"Loaded {len(traj_frames)} frames from {traj_file} (box = {box_len:.3f} A)")

# 2. Compute Oxygen-Oxygen Radial Distribution Function g_OO(r)
rmax = 4.80
nbins = 120
rdf_per_frame = []
r_grid = None

for frame in traj_frames:
    g_r, r_grid = get_rdf(frame, rmax=rmax, nbins=nbins, elements=(8, 8))
    rdf_per_frame.append(g_r)

g_oo = np.mean(rdf_per_frame, axis=0)
g_oo_std = np.std(rdf_per_frame, axis=0) / np.sqrt(len(traj_frames))  # standard error

mask_1st = (r_grid > 2.4) & (r_grid < 3.3)
peak_idx = np.argmax(g_oo[mask_1st])
r1_sim = r_grid[mask_1st][peak_idx]
g1_sim = g_oo[mask_1st][peak_idx]

# 3. Compute Oxygen Mean-Squared Displacement MSD(t) and Self-Diffusion Coefficient D
o_indices = [i for i, sym in enumerate(traj_frames[0].get_chemical_symbols()) if sym == "O"]
pos_series = np.array([f.get_positions()[o_indices] for f in traj_frames])

diffs = np.diff(pos_series, axis=0)
diffs -= box_len * np.round(diffs / box_len)
unwrapped = np.concatenate([pos_series[0:1], pos_series[0:1] + np.cumsum(diffs, axis=0)], axis=0)

msd = np.mean(np.sum((unwrapped - unwrapped[0]) ** 2, axis=-1), axis=1)  # A^2
times_ps = np.arange(len(traj_frames)) * dt_frame_ps

# Linear Einstein fit over diffusive regime (t >= 0.50 ps): MSD(t) = 6 * D * t + c
fit_mask = times_ps >= 0.50
slope, intercept = np.polyfit(times_ps[fit_mask], msd[fit_mask], 1)
D_sim = (slope / 6.0) * 10.0  # in 10^-5 cm^2/s

# Experimental NMR slope reference (D_expt = 2.30 x 10^-5 cm^2/s -> slope = 6 * 0.230 = 1.380 A^2/ps)
D_expt = 2.30
slope_expt = 6.0 * (D_expt / 10.0)
msd_expt_line = intercept + slope_expt * times_ps

# 4. Plot Styling Helper
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 9.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
})

def plot_rdf(ax):
    ax.plot(r_grid, g_oo, color="#1f77b4", lw=2.2, label=f"MACE-OFF23 (300 K): $r_1={r1_sim:.2f}$ Å, $g(r_1)={g1_sim:.2f}$")
    ax.fill_between(r_grid, g_oo - g_oo_std, g_oo + g_oo_std, color="#1f77b4", alpha=0.2)
    ax.axhline(1.0, color="gray", ls=":", lw=1.0)
    ax.axvline(2.80, color="#d62728", ls="--", lw=1.4, label="Expt. X-ray/Neutron ($r_1 = 2.80$ Å)")
    ax.axhspan(2.57, 2.80, color="#d62728", alpha=0.12, label="Expt. $g(r_1)$ range (2.57–2.80)")
    ax.scatter([r1_sim], [g1_sim], color="#1f77b4", zorder=5, s=40)
    ax.set_xlim(2.0, rmax)
    ax.set_ylim(0.0, max(3.3, g1_sim + 0.35))
    ax.set_xlabel(r"O$\cdots$O Distance $r$ (Å)")
    ax.set_ylabel(r"Radial Distribution Function $g_{\mathrm{OO}}(r)$")
    ax.set_title("Liquid Water Structure: $g_{\mathrm{OO}}(r)$")
    ax.legend(loc="upper right", frameon=True)

def plot_msd(ax):
    ax.plot(times_ps, msd, color="#2ca02c", lw=2.2, label="MACE-OFF23 Oxygen MSD(t)")
    ax.plot(
        times_ps[fit_mask],
        slope * times_ps[fit_mask] + intercept,
        color="#111111",
        ls="--",
        lw=1.8,
        label=f"Einstein fit ($t \\geq 0.5$ ps): $D = {D_sim:.2f} \\times 10^{{-5}}$ cm$^2$/s",
    )
    ax.plot(
        times_ps[fit_mask],
        msd_expt_line[fit_mask],
        color="#d62728",
        ls=":",
        lw=1.8,
        label=f"Expt. NMR slope: $D = {D_expt:.2f} \\times 10^{{-5}}$ cm$^2$/s",
    )
    ax.axvline(0.50, color="gray", ls="-.", lw=1.0, alpha=0.7)
    ax.set_xlim(0.0, times_ps[-1])
    ax.set_ylim(0.0, max(msd.max(), msd_expt_line.max()) * 1.12)
    ax.set_xlabel("Simulation Time $t$ (ps)")
    ax.set_ylabel(r"Mean-Squared Displacement $\langle |\Delta \mathbf{r}_{\mathrm{O}}(t)|^2 \rangle$ (Å$^2$)")
    ax.set_title("Liquid Water Dynamics: Self-Diffusion")
    ax.legend(loc="upper left", frameon=True)

# 5. Export combined two-panel figure + individual figures
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.8), dpi=300)
plot_rdf(ax1)
plot_msd(ax2)
fig.tight_layout()
fig.savefig("water_md_analysis.png", dpi=300)
plt.close(fig)

fig1, ax_single1 = plt.subplots(figsize=(6.2, 4.6), dpi=300)
plot_rdf(ax_single1)
fig1.tight_layout()
fig1.savefig("g_oo_rdf.png", dpi=300)
plt.close(fig1)

fig2, ax_single2 = plt.subplots(figsize=(6.2, 4.6), dpi=300)
plot_msd(ax_single2)
fig2.tight_layout()
fig2.savefig("msd_diffusion.png", dpi=300)
plt.close(fig2)

print("Saved plots: water_md_analysis.png, g_oo_rdf.png, msd_diffusion.png")


