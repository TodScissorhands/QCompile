#!/usr/bin/env python3
"""Generate and internally cross-check the fixed HeH+ STO-3G/JW reference."""
import json
import platform
from pathlib import Path

import numpy as np
import pyscf
from pyscf import ao2mo, gto, scf
import qiskit
import qiskit_nature
from qiskit_nature.second_q.drivers import InitialGuess, MethodType, PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.units import DistanceUnit

ROOT = Path(__file__).resolve().parent
R_BOHR = 1.463283
NMODES = 4


def ladder_matrix(terms, dimension=1 << NMODES):
    """Build a Fock matrix directly by applying rightmost ladder operators first."""
    result = np.zeros((dimension, dimension), dtype=np.complex128)
    for label, coefficient in terms:
        operators = [] if not label else label.split()
        for ket in range(dimension):
            state = ket
            amplitude = complex(coefficient)
            for token in reversed(operators):
                action, index_text = token.split("_")
                index = int(index_text)
                occupied = (state >> index) & 1
                if (action == "+" and occupied) or (action == "-" and not occupied):
                    amplitude = 0.0
                    break
                if (state & ((1 << index) - 1)).bit_count() & 1:
                    amplitude = -amplitude
                state ^= 1 << index
            result[state, ket] += amplitude
    return result

def integral_ladder_terms(h1, eri):
    """Expand chemist-ordered MO integrals into spin-orbital ladder terms."""
    norb = h1.shape[0]
    terms = []
    for spin in range(2):
        for p in range(norb):
            for q in range(norb):
                terms.append((f"+_{p + spin*norb} -_{q + spin*norb}", h1[p, q]))
    for spin_p in range(2):
        for spin_r in range(2):
            for p in range(norb):
                for q in range(norb):
                    for r in range(norb):
                        for s in range(norb):
                            coefficient = 0.5 * eri[p, q, r, s]
                            label = (
                                f"+_{p + spin_p*norb} +_{r + spin_r*norb} "
                                f"-_{s + spin_r*norb} -_{q + spin_p*norb}"
                            )
                            terms.append((label, coefficient))
    return terms


def matrix_from_integrals(h1, eri):
    return ladder_matrix(integral_ladder_terms(h1, eri))


def as_json_number(value):
    z = complex(value)
    if abs(z.imag) < 1e-15:
        return float(z.real)
    return {"real": float(z.real), "imag": float(z.imag)}


def main():
    atom = "He 0.0 0.0 0.0; H 0.0 0.0 1.463283"
    driver = PySCFDriver(
        atom=atom,
        unit=DistanceUnit.BOHR,
        charge=1,
        spin=0,
        basis="sto3g",
        method=MethodType.RHF,
        conv_tol=1e-12,
        max_cycle=100,
        init_guess=InitialGuess.MINAO,
        max_memory=2000,
    )
    problem = driver.run()
    fermionic = problem.hamiltonian.second_q_op()
    mapped = JordanWignerMapper().map(fermionic)
    pauli_terms = [(label, as_json_number(coefficient)) for label, coefficient in mapped.to_list()]
    h_pauli = mapped.to_matrix()
    h_fermion = ladder_matrix(list(fermionic.items()))
    matrix_residual = float(np.max(np.abs(h_pauli - h_fermion)))
    hermiticity_residual = float(np.max(np.abs(h_fermion - h_fermion.conj().T)))
    if matrix_residual > 2e-12 or hermiticity_residual > 2e-12:
        raise AssertionError(f"matrix residuals too large: map={matrix_residual}, herm={hermiticity_residual}")

    mol = gto.M(
        atom=atom,
        unit="Bohr",
        basis="sto-3g",
        charge=1,
        spin=0,
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.max_cycle = 100
    mf.init_guess = "minao"
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("independent explicit PySCF RHF did not converge")
    nmo = mf.mo_coeff.shape[1]
    h1_mo = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri_mo = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), nmo)
    qn_h1 = np.asarray(problem.hamiltonian.electronic_integrals.alpha["+-"])
    h1_residual = float(np.max(np.abs(h1_mo - qn_h1)))
    qn_nuclear = float(problem.nuclear_repulsion_energy)
    pyscf_nuclear = float(mol.energy_nuc())
    nuclear_residual = abs(qn_nuclear - pyscf_nuclear)
    integral_terms = integral_ladder_terms(h1_mo, eri_mo)
    h_integrals = ladder_matrix(integral_terms)
    integral_vs_fermion_residual = float(np.max(np.abs(h_integrals - h_fermion)))
    integral_vs_pauli_residual = float(np.max(np.abs(h_integrals - h_pauli)))
    if max(integral_vs_fermion_residual, integral_vs_pauli_residual) > 2e-12:
        raise AssertionError(
            "integral/fermion/Pauli matrix mismatch: "
            f"{integral_vs_fermion_residual}, {integral_vs_pauli_residual}"
        )

    # Restrict to N_alpha=N_beta=1, then independently select S=0 by the
    # eigenvalues of S^2 constructed from the same fermionic ladder rules.
    sector = [state for state in range(1 << NMODES)
              if state.bit_count() == 2
              and (state & 0b0011).bit_count() == 1
              and (state & 0b1100).bit_count() == 1]
    plus_terms = [(f"+_{i} -_{i + 2}", 1.0) for i in range(2)]
    minus_terms = [(f"+_{i + 2} -_{i}", 1.0) for i in range(2)]
    s_plus = ladder_matrix(plus_terms)
    s_minus = ladder_matrix(minus_terms)
    s2 = 0.5 * (s_plus @ s_minus + s_minus @ s_plus)
    s2_sector = s2[np.ix_(sector, sector)]
    spin_values, spin_vectors = np.linalg.eigh(s2_sector)
    singlet = spin_vectors[:, np.abs(spin_values) < 1e-10]
    if singlet.shape[1] != 3:
        raise AssertionError(f"expected 3 singlet CSFs in the Ms=0 determinant sector; got {singlet.shape[1]}")
    h_sector = h_integrals[np.ix_(sector, sector)]
    h_singlet = singlet.conj().T @ h_sector @ singlet
    electronic = float(np.linalg.eigvalsh(h_singlet)[0].real)
    nuclear = pyscf_nuclear
    total = electronic + nuclear

    # Independently verify qubit convention: Qiskit Pauli label's rightmost
    # character acts on q0, whose occupation is bit 0 (least significant).
    # The full matrix comparison above uses basis index=sum(bit[q]*2**q).
    payload = {
        "benchmark": {
            "system": "HeH+", "geometry": atom, "distance_bohr": R_BOHR,
            "basis": "STO-3G", "charge": 1, "electrons": 2,
            "spin": 0, "state": "singlet S=0", "method": "RHF integrals; FCI by determinant diagonalization",
            "orbitals_retained": nmo, "spin_orbitals": NMODES, "qubits": NMODES,
            "mapping": "Jordan-Wigner", "symmetry_reduction": "none",
            "spin_orbital_order": ["0 alpha", "1 alpha", "0 beta", "1 beta"],
            "basis_bit_order": "q0 is least-significant bit; occupation bit q equals spin-orbital q",
            "pauli_label_order": "Qiskit label left-to-right is q3..q0; rightmost character is q0",
            "hamiltonian_energy_convention": "electronic only; nuclear repulsion excluded",
        },
        "settings": {
            "rhf_conv_tol_hartree": 1e-12, "rhf_max_cycle": 100,
            "initial_guess": "MINAO", "frozen_core": False,
            "precision": "binary64 full precision JSON decimal serialization",
            "checks_abs_tolerance": 2e-12,
            "generator_versions": {
                "python": platform.python_version(), "numpy": np.__version__,
                "pyscf": pyscf.__version__, "qiskit": qiskit.__version__,
                "qiskit_nature": qiskit_nature.__version__,
            },
        },
        "energies_hartree": {
            "rhf_total_including_nuclear": float(mf.e_tot),
            "rhf_electronic": float(mf.e_tot - pyscf_nuclear),
            "nuclear_repulsion": nuclear,
            "fci_electronic_singlet": electronic,
            "fci_total_born_oppenheimer": total,
            "singlet_sector_dimension": int(singlet.shape[1]),
            "Ms0_determinant_sector_dimension": len(sector),
            "Ms0_s2_eigenvalues": spin_values.tolist(),
        },
        "integrals_hartree": {
            "one_electron_mo_h_pq": h1_mo.tolist(),
            "two_electron_mo_chemist_eri_pqrs": eri_mo.tolist(),
            "nuclear_repulsion": nuclear,
            "qn_driver_one_electron_max_abs_residual_vs_direct_pyscf": h1_residual,
            "qn_driver_nuclear_max_abs_residual_vs_direct_pyscf": nuclear_residual,
        },
        "fermionic_operator_terms": [
            {"ladder_string": label, "coefficient": as_json_number(coefficient)}
            for label, coefficient in fermionic.items()
        ],
        "pauli_hamiltonian_terms_hartree": [
            {"label_q3_to_q0": label, "coefficient": coefficient}
            for label, coefficient in pauli_terms
        ],
        "matrix_validation": {
            "method": "independent 16x16 determinant matrix from frozen h_pq and chemist-ordered (pq|rs) integrals, compared to Qiskit SparsePauliOp.to_matrix and Qiskit Nature FermionicOp ladder terms",
            "maximum_absolute_pauli_vs_qn_fermionic_matrix_residual_hartree": matrix_residual,
            "maximum_absolute_integral_vs_qn_fermionic_matrix_residual_hartree": integral_vs_fermion_residual,
            "maximum_absolute_integral_vs_pauli_matrix_residual_hartree": integral_vs_pauli_residual,
            "maximum_absolute_hermiticity_residual_hartree": hermiticity_residual,
            "exact_singlet_energy_from_direct_integral_determinant_matrix_hartree": electronic,
        },
    }
    (ROOT / "reference_data.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload["energies_hartree"], indent=2))
    print(json.dumps(payload["matrix_validation"], indent=2))
    print("Pauli terms:")
    for term in payload["pauli_hamiltonian_terms_hartree"]:
        print(term)


if __name__ == "__main__":
    main()
