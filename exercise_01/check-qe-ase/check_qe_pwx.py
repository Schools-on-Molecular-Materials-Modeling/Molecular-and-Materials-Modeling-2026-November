import os
import shutil
from pathlib import Path
import numpy as np
from ase.build import molecule
from ase.calculators.espresso import Espresso, EspressoProfile

nprocs = 2
print(f"=== Quantum ESPRESSO (pw.x, np={nprocs}) Smoke Test ===")
pwx_path = shutil.which("pw.x")
assert pwx_path, "pw.x executable not found in PATH"
print(f"Executable         : {pwx_path}")

# Locate H and O pseudopotential (.UPF / .upf) files
search_dirs = [
    os.environ.get("ESPRESSO_PSEUDO", ""),
    "./pseudo",
    str(Path.home() / "pseudo"),
    str(Path.home() / "pseudopotentials"),
]
pseudo_dir, h_upf, o_upf = None, None, None
for d in search_dirs:
    if d and Path(d).is_dir():
        p = Path(d)
        h_matches = sorted(p.glob("H*.[uU][pP][fF]"))
        o_matches = sorted(p.glob("O*.[uU][pP][fF]"))
        if h_matches and o_matches:
            pseudo_dir = str(p.resolve())
            h_upf = h_matches[0].name
            o_upf = o_matches[0].name
            break

if not pseudo_dir:
    raise FileNotFoundError(
        "pw.x is installed, but no H/O .upf files were found in "
        "$ESPRESSO_PSEUDO, ./pseudo, or ~/pseudo."
    )

print(f"Pseudo Directory   : {pseudo_dir} (H={h_upf}, O={o_upf})")
atoms = molecule("H2O")
atoms.center(vacuum=6.0)
atoms.pbc = True

profile = EspressoProfile(
    command=f"mpirun -np {nprocs} pw.x -in PREFIX.pwi > PREFIX.pwo",
    pseudo_dir=pseudo_dir,
)
atoms.calc = Espresso(
    profile=profile,
    input_data={
        "control": {"calculation": "scf", "tprnfor": True, "outdir": "./tmp_qe_check"},
        "system": {"ecutwfc": 30.0, "ecutrho": 240.0},
        "electrons": {"conv_thr": 1e-6},
    },
    pseudopotentials={"H": h_upf, "O": o_upf},
    kpts=(1, 1, 1),
)

try:
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    print("Status             : PASS")
    print(f"Total Energy (H2O) : {energy:.6f} eV")
    print(f"Max Force Component: {np.abs(forces).max():.6f} eV/A")
finally:
    for fname in ("espresso.pwi", "espresso.pwo", "espresso.err", "CRASH"):
        p = Path(fname)
        if p.is_file():
            p.unlink()
    if Path("tmp_qe_check").is_dir():
        shutil.rmtree("tmp_qe_check")

