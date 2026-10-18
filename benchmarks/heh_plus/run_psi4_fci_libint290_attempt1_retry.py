#!/usr/bin/env python3
"""No-overwrite Psi4 FCI retry with an explicitly created scratch directory."""
import json
import platform
from pathlib import Path

import psi4

ROOT = Path(__file__).resolve().parent
R_BOHR = 1.463283
OUT = ROOT / "psi4_fci_output_libint290_attempt1_retry.txt"
RESULT = ROOT / "psi4_fci_reference_libint290_attempt1_retry.json"
SCRATCH = ROOT / ".psi4-scratch-libint290-attempt1-retry"
SCRATCH.mkdir(parents=True, exist_ok=False)
psi4.core.set_output_file(str(OUT), False)
psi4.core.set_num_threads(1)
psi4.set_memory("1 GB")
psi4.core.IOManager.shared_object().set_default_path(str(SCRATCH))
psi4.set_options({
    "basis": "sto-3g",
    "reference": "rhf",
    "scf_type": "pk",
    "e_convergence": 1e-12,
    "d_convergence": 1e-12,
    "maxiter": 100,
    "num_roots": 1,
    "r_convergence": 1e-10,
    "ci_maxiter": 100,
    "calc_s_squared": True,
    "freeze_core": False,
})
mol = psi4.geometry(f"""
1 1
He 0.0 0.0 0.0
H  0.0 0.0 {R_BOHR:.12f}
units au
symmetry c1
no_reorient
no_com
""")
energy, wfn = psi4.energy("fci", molecule=mol, return_wfn=True)
nuclear = float(psi4.variable("NUCLEAR REPULSION ENERGY"))
electronic = float(energy - nuclear)
try:
    s_squared = float(wfn.variable("S^2"))
except Exception:
    try:
        s_squared = float(psi4.variable("S^2"))
    except Exception:
        s_squared = None

if s_squared is None:
    spin_lines = [
        line for line in OUT.read_text(encoding="utf-8").splitlines()
        if "Computed <S^2> vector 0 =" in line
    ]
    if spin_lines:
        s_squared = float(spin_lines[-1].split("=")[-1].strip())

payload = {
    "program": "Psi4",
    "version": psi4.__version__,
    "python": platform.python_version(),
    "method": "FCI (Psi4 DETCI), RHF closed-shell singlet reference",
    "geometry": {"atoms": [["He", 0.0, 0.0, 0.0], ["H", 0.0, 0.0, R_BOHR]], "unit": "bohr", "distance_bohr": R_BOHR},
    "basis": "STO-3G",
    "charge": 1,
    "multiplicity": 1,
    "electron_count": 2,
    "nalpha": 1,
    "nbeta": 1,
    "spatial_orbitals": int(wfn.nmo()),
    "frozen_docc": 0,
    "frozen_uocc": 0,
    "energies_hartree": {
        "electronic_fci": electronic,
        "nuclear_repulsion": nuclear,
        "total_born_oppenheimer": float(energy),
    },
    "spin_diagnostic_s_squared": s_squared,
    "settings": {
        "SCF": "RHF; scf_type=pk; e_convergence=1e-12; d_convergence=1e-12; maxiter=100",
        "FCI": "num_roots=1; r_convergence=1e-10; ci_maxiter=100; CALC_S_SQUARED=true",
        "symmetry": "C1", "threads": 1, "memory": "1 GB",
        "energy_convention": "Psi4 FCI total energy includes nuclear repulsion; electronic energy is total minus nuclear repulsion",
        "scratch_directory": str(SCRATCH),
    },
}
RESULT.write_text(json.dumps(payload, indent=2) + "\n")
print(json.dumps(payload, indent=2))
