# The Coronene Dimer


Coronene Dimer (\(\text{C}_{24}\text{H}_{12}\), 72 Atoms) Is an Ideal Test

• Pure dispersion regime: Two parallel, planar polycyclic aromatic hydrocarbons (PAHs) have zero net dipole moment; without accurate dispersion physics, the interaction energy is purely repulsive.


## Challenge

With the help of AI agents, adapt the van der Waals scan scripts (`MACE-OFF23`, `MACE-MP-0`, `GFN2-xTB`, and `MOPAC PM6-D3H4` / `PM7`) to benchmark the **benzene dimer ($\text{C}_6\text{H}_6 \cdots \text{C}_6\text{H}_6$, 24 atoms)**—the foundational aromatic $\pi\text{–}\pi$ interaction system from the **S22** and **S66** non-covalent benchmark datasets.

### Reference Benchmark Data (CCSD(T)/CBS, S22 / S66 Datasets)

Unlike the larger coronene dimer (where dispersion dominates over electrostatics and favors parallel stacking), the benzene dimer exhibits a delicate balance between London dispersion and quadrupole–quadrupole electrostatics, resulting in two nearly isoenergetic minima:

| Stacking Configuration | Equilibrium Geometry | Reference $E_{\text{int}}$ ($\text{kcal/mol}$) | Reference $E_{\text{int}}$ ($\text{kJ/mol}$) | Physical Role |
| :--- | :--- | :---: | :---: | :--- |
| **Parallel-Displaced (Slip-Stacked, $C_{2h}$)** | $z_{\text{eq}} \approx 3.40\text{–}3.50\text{ \AA}$ *(lateral slip $\approx 1.6\text{ \AA}$)* | **$-2.73\text{ kcal/mol}$** | **$-11.42\text{ kJ/mol}$** | Co-winner global minimum (dispersion + favorable offset quadrupoles) |
| **T-Shaped (Edge-to-Face, $C_{2v}$)** | $R_{\text{COM}} \approx 4.95\text{ \AA}$ *(ring center-to-center)* | **$-2.74\text{ kcal/mol}$** | **$-11.46\text{ kJ/mol}$** | Co-winner global minimum ($\text{C–H}^{\delta+}\cdots\pi^{\delta-}$ quadrupole attraction) |
| **Sandwich (Eclipsed, $D_{6h}$)** | $z_{\text{eq}} \approx 3.80\text{–}3.90\text{ \AA}$ | $-1.70\text{ kcal/mol}$ | $-7.11\text{ kJ/mol}$ | Repulsive quadrupole–quadrupole saddle point ($\sim 4.3\text{ kJ/mol}$ higher) |

### Suggested Tasks

1. **Parallel-Displaced vs. Sandwich Scan:** Build a benzene monomer (`ase.build.molecule("C6H6")`), relax it tightly, and scan the vertical interplanar separation $z \in [3.0, 6.0]\text{ \AA}$ for both the **eclipsed sandwich** ($\text{slip}_x = 0.0\text{ \AA}$) and **parallel-displaced** ($\text{slip}_x = 1.6\text{ \AA}$) configurations.
2. **Cutoff Sensitivity Test (`MACE-OFF23` vs. Coronene):** Notice that for the smaller benzene dimer (ring diameter $\approx 2.8\text{ \AA}$), nearly all inter-ring $\text{C}\cdots\text{C}$ pairs lie inside `MACE-OFF23`'s $r_{\text{cut}} = 5.0\text{ \AA}$ cutoff—verify whether `MACE-OFF23` recovers close to $100\%$ of the **$-11.42\text{ kJ/mol}$** CCSD(T)/CBS benchmark (unlike the $68\%$ recovery observed for the $9.5\text{ \AA}$-wide coronene dimer).
3. **T-Shaped Configuration (Bonus):** Rotate the top benzene ring by $90^\circ$ around the $x$-axis (`top.rotate(90, "x")`), place its center of mass at $z \in [4.2, 6.5]\text{ \AA}$ above the bottom ring, and check whether each method correctly predicts the **T-shaped** and **parallel-displaced** dimers to be nearly degenerate ($\Delta E < 1\text{ kJ/mol}$) and more stable than the **sandwich** dimer.

