#!/usr/bin/env python3
"""Independently validate the fixed LiH (1 core, 2 active orbital) CASCI case."""
import json
import math

from pyscf import gto, mcscf, scf

DISTANCE_ANGSTROM = 1.595
ANGSTROM_TO_BOHR = 0.52917721092
EXPECTED_RHF_TOTAL = -7.86202386012712
EXPECTED_CORE_ENERGY = -6.80297354998635
EXPECTED_ACTIVE_ENERGY = -1.059311983913
EXPECTED_CASCI_TOTAL = -7.86228553389953
REFERENCE_TOLERANCE = 1.0e-7
IDENTITY_TOLERANCE = 1.0e-12


def require_close(name, actual, expected, tolerance):
    difference = abs(actual - expected)
    if not math.isfinite(actual) or difference > tolerance:
        raise RuntimeError(
            f"{name} differs from reference: actual={actual:.15g}, "
            f"expected={expected:.15g}, difference={difference:.3e}, "
            f"tolerance={tolerance:.3e}"
        )
    return difference


def main():
    distance_bohr = DISTANCE_ANGSTROM / ANGSTROM_TO_BOHR
    molecule = gto.M(
        atom=f"Li 0 0 0; H 0 0 {distance_bohr:.16g}",
        unit="Bohr",
        basis="sto-3g",
        charge=0,
        spin=0,
        verbose=0,
    )
    rhf = scf.RHF(molecule)
    rhf.conv_tol = 1.0e-12
    rhf.max_cycle = 100
    rhf.init_guess = "minao"
    rhf.kernel()
    if not rhf.converged:
        raise RuntimeError("PySCF RHF did not converge")

    casci = mcscf.CASCI(rhf, 2, 2)
    casci.ncore = 1
    casci_result = casci.kernel()
    _, core_energy = casci.get_h1eff()
    active_energy = float(casci_result[1])
    total_energy = float(casci_result[0])
    nuclear_repulsion = float(molecule.energy_nuc())
    rhf_total = float(rhf.e_tot)
    identity_residual = abs(total_energy - (core_energy + active_energy))

    if casci.ncore != 1 or casci.ncas != 2 or tuple(casci.nelecas) != (1, 1):
        raise RuntimeError("PySCF CASCI did not retain the requested (1 core, 2 active, 1 alpha/1 beta) model")
    if identity_residual > IDENTITY_TOLERANCE:
        raise RuntimeError(f"CASCI energy identity residual is {identity_residual:.3e} Ha")

    differences = {
        "rhf_total_abs_difference_hartree": require_close(
            "RHF total energy", rhf_total, EXPECTED_RHF_TOTAL, REFERENCE_TOLERANCE
        ),
        "core_energy_abs_difference_hartree": require_close(
            "frozen-core scalar", core_energy, EXPECTED_CORE_ENERGY, REFERENCE_TOLERANCE
        ),
        "active_energy_abs_difference_hartree": require_close(
            "active-space CASCI energy", active_energy, EXPECTED_ACTIVE_ENERGY, REFERENCE_TOLERANCE
        ),
        "casci_total_abs_difference_hartree": require_close(
            "CASCI total energy", total_energy, EXPECTED_CASCI_TOTAL, REFERENCE_TOLERANCE
        ),
    }
    print(json.dumps({
        "program": "PySCF CASCI",
        "pyscf_version": __import__("pyscf").__version__,
        "geometry": {
            "formula": "LiH",
            "distance_angstrom": DISTANCE_ANGSTROM,
            "distance_bohr": distance_bohr,
            "basis": "STO-3G",
            "charge": 0,
            "electrons": 4,
        },
        "method": {
            "reference": "RHF",
            "active_space": {"frozen_core_orbitals": 1, "active_spatial_orbitals": 2, "active_electrons": [1, 1]},
            "rhf_energy_convergence_hartree": 1.0e-12,
            "rhf_max_iterations": 100,
            "initial_guess": "MINAO",
        },
        "energies_hartree": {
            "rhf_total": rhf_total,
            "nuclear_repulsion": nuclear_repulsion,
            "core_energy_including_nuclear": float(core_energy),
            "active_space_casci_energy": active_energy,
            "casci_total": total_energy,
        },
        "energy_identity_abs_residual_hartree": identity_residual,
        "reference_tolerance_hartree": REFERENCE_TOLERANCE,
        "reference_abs_differences_hartree": differences,
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
