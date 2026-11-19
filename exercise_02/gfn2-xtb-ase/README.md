# Methane Combustion Enthalpy Calculation with GFN2-xTB

The `methane_combustion_xtb.py` benchmark demonstrates a classic characteristic of extended tight-binding methods: **GFN2-xTB** reproduces molecular geometries and vibrational frequencies accurately, but severely overestimates the combustion exothermicity ($-1222.49\text{ kJ/mol}$ vs. $-802.30\text{ kJ/mol}$ experimental).

## Physical Breakdown of the GFN2-xTB Results

1. **Near-Experimental Geometries & Frequencies:**
   - **Triplet $\text{O}_2$:** Setting `uhf=2` gives an $\text{O=O}$ bond length of $1.2098\text{ \AA}$ (Expt: $1.208\text{ \AA}$) and a stretching frequency of $1629.4\text{ cm}^{-1}$ (Expt: $1580\text{ cm}^{-1}$).
   - **$\text{CH}_4$, $\text{CO}_2$, and $\text{H}_2\text{O}$:** All bond lengths fall within $0.016\text{ \AA}$ of experiment, and $\Delta\text{ZPE} = +8.74\text{ kJ/mol}$ closely matches the true zero-point shift ($\approx +8$ to $+10\text{ kJ/mol}$).

2. **Why GFN2-xTB Overestimates Combustion Exothermicity by $\sim 420\text{ kJ/mol}$:**
   - **Focus of the GFN family:** As indicated by its name (**G**eometries, **F**requencies, **N**on-covalent interactions), Grimme's GFN2-xTB is parameterized specifically for structures, vibrational modes, and intermolecular interactions—**not** absolute covalent bond breaking/formation thermochemistry involving $\text{O}_2$.
   - **Minimal valence basis & $\text{O}_2$ underbinding:** Using a minimal valence basis without polarization functions on oxygen leaves high-spin triplet $\text{O}_2$ strongly underbound relative to closed-shell $\text{CO}_2$ and $\text{H}_2\text{O}$, artificially driving $\Delta_c H^\circ(298\text{ K})$ down to $-1222.49\text{ kJ/mol}$.


