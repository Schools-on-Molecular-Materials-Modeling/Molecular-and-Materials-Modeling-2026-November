import numpy as np
from ase import Atoms
from ase.units import kJ, mol
from xtb.ase.calculator import XTB


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


# --- Monomer reference (fresh calculator) ---
mono = build_coronene()
mono.calc = XTB(method="GFN2-xTB")
e_mono = mono.get_potential_energy()
print(f"E_mono = {e_mono:.6f} eV")

# --- Dimer at large separation, fresh calculator for every object ---
for z in [6.5, 6.0, 5.5, 5.0, 4.5, 4.0]:
    top = mono.copy()
    top.translate([1.42, 0.0, z])
    top.calc = XTB(method="GFN2-xTB")          # fresh

    dimer = mono + top
    dimer.pbc = False
    dimer.calc = XTB(method="GFN2-xTB")        # fresh, NOT reused

    e_dimer = dimer.get_potential_energy()
    e_int = (e_dimer - 2.0 * e_mono) / (kJ / mol)
    print(f"z = {z:.2f} A   E_int = {e_int:+8.2f} kJ/mol")
