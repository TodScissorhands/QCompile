#!/usr/bin/env python3
"""Run a separate Psi4 singlet-FCI calculation for H2 at a fixed distance."""
import argparse
import json
import math
import os
import platform
import tempfile
from pathlib import Path

ANGSTROM_TO_BOHR = 1.0 / 0.52917721092


def write_new_json(path, payload):
    path = path.resolve()
    if not path.parent.is_dir():
        raise FileNotFoundError(f"output directory does not exist: {path.parent}")
    if path.exists():
        raise FileExistsError(f"refusing to overwrite: {path}")
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent,
        prefix=f".{path.name}.", suffix=".tmp", delete=False,
    ) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(payload, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
            os.link(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance-angstrom", type=float, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--generator-provenance", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.distance_angstrom) or args.distance_angstrom <= 0.0:
        parser.error("distance must be finite and greater than zero")
    output = args.output_json.resolve()
    log_path = Path(str(output) + ".psi4.log")
    scratch = args.scratch.resolve()
    if output.exists() or log_path.exists() or scratch.exists():
        parser.error("refusing to overwrite output, log, or scratch path")
    if not output.parent.is_dir():
        parser.error(f"output parent directory does not exist: {output.parent}")

    import psi4

    distance_bohr = args.distance_angstrom * ANGSTROM_TO_BOHR
    scratch.mkdir(parents=False, exist_ok=False)
    psi4.core.set_output_file(str(log_path), False)
    psi4.core.set_num_threads(1)
    psi4.set_memory("1 GB")
    psi4.core.IOManager.shared_object().set_default_path(str(scratch))
    psi4.set_options({
        "basis": "sto-3g",
        "reference": "rhf",
        "scf_type": "pk",
        "e_convergence": 1.0e-12,
        "d_convergence": 1.0e-12,
        "maxiter": 100,
        "num_roots": 1,
        "r_convergence": 1.0e-10,
        "ci_maxiter": 100,
        "calc_s_squared": True,
        "freeze_core": False,
    })
    molecule = psi4.geometry(f"""
0 1
H 0.0 0.0 0.0
H 0.0 0.0 {distance_bohr:.16f}
units au
symmetry c1
no_reorient
no_com
""")
    molecule.set_molecular_charge(0)
    molecule.set_multiplicity(1)
    molecule.update_geometry()
    total, wavefunction = psi4.energy("fci", molecule=molecule, return_wfn=True)
    nuclear = float(psi4.variable("NUCLEAR REPULSION ENERGY"))
    total = float(total)
    electronic = total - nuclear
    try:
        s_squared = float(wavefunction.variable("S^2"))
    except Exception:
        try:
            s_squared = float(psi4.variable("S^2"))
        except Exception:
            s_squared = None
    if s_squared is None:
        spin_lines = [
            line for line in Path(log_path).read_text(encoding="utf-8").splitlines()
            if "Computed <S^2> vector 0 =" in line
        ]
        if spin_lines:
            s_squared = float(spin_lines[-1].split("=")[-1].strip())
    nalpha = int(wavefunction.nalpha())
    nbeta = int(wavefunction.nbeta())
    electron_count = nalpha + nbeta
    if electron_count != 2 or nalpha != 1 or nbeta != 1:
        raise RuntimeError("Psi4 did not construct the required two-electron singlet sector")
    if wavefunction.nmo() != 2:
        raise RuntimeError(f"expected two spatial orbitals, got {wavefunction.nmo()}")
    if s_squared is not None and abs(s_squared) > 1.0e-8:
        raise RuntimeError(f"Psi4 singlet spin diagnostic is not zero: <S^2>={s_squared}")

    generator = json.loads(args.generator_provenance.read_text(encoding="utf-8"))
    if abs(generator["model"]["distance_angstrom"] - args.distance_angstrom) > 1.0e-12:
        raise RuntimeError("generator provenance geometry does not match Psi4 input")
    generator_energies = generator["energies_hartree"]
    generator_electronic = float(generator_energies["electronic_singlet_pyscf_fci"])
    generator_nuclear = float(generator_energies["nuclear_repulsion"])
    generator_total = float(generator_energies["total_born_oppenheimer_pyscf_fci"])
    comparison = {
        "generator_provenance": str(args.generator_provenance.resolve()),
        "electronic_abs_difference_hartree": abs(electronic - generator_electronic),
        "electronic_tolerance_hartree": 1.0e-8,
        "nuclear_abs_difference_hartree": abs(nuclear - generator_nuclear),
        "nuclear_tolerance_hartree": 1.0e-12,
        "total_abs_difference_hartree": abs(total - generator_total),
        "total_tolerance_hartree": 1.0e-8,
    }
    if (comparison["electronic_abs_difference_hartree"] > 1.0e-8
            or comparison["nuclear_abs_difference_hartree"] > 1.0e-12
            or comparison["total_abs_difference_hartree"] > 1.0e-8):
        raise RuntimeError(f"Psi4 and generated H2 reference disagree: {comparison}")
    payload = {
        "schema": "QCOMPILE-H2-Psi4-FCI-1",
        "program": "Psi4 DETCI FCI",
        "version": psi4.__version__,
        "python": platform.python_version(),
        "model": {
            "system": "H2", "distance_angstrom": args.distance_angstrom,
            "distance_bohr": distance_bohr, "geometry_units": "bohr",
            "basis": "STO-3G", "charge": 0, "electrons": electron_count,
            "multiplicity": 1, "spin_s": 0, "nalpha": nalpha, "nbeta": nbeta,
            "spatial_orbitals": int(wavefunction.nmo()), "frozen_orbitals": 0,
        },
        "method": "RHF reference; Psi4 DETCI full configuration interaction",
        "settings": {
            "scf_type": "PK", "scf_energy_convergence_hartree": 1.0e-12,
            "scf_density_convergence": 1.0e-12, "scf_max_iterations": 100,
            "fci_roots": 1, "fci_residual_convergence": 1.0e-10,
            "fci_max_iterations": 100, "calc_s_squared": True,
            "frozen_core": False, "symmetry": "C1", "threads": 1,
            "memory": "1 GB",
            "energy_convention": "Psi4 energy includes nuclear repulsion; electronic energy = total - nuclear",
            "scratch_directory": str(scratch), "log_file": str(log_path),
        },
        "energies_hartree": {
            "electronic_fci": electronic,
            "nuclear_repulsion": nuclear,
            "total_born_oppenheimer": total,
        },
        "spin_diagnostic_s_squared": s_squared,
        "comparison_to_generator": comparison,
    }
    result_path = write_new_json(output, payload)
    print(json.dumps({"result": str(result_path), "log": str(log_path), **payload}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
