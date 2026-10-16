#!/usr/bin/env python3
"""Generate a geometry-specific H2 STO-3G two-qubit Hamiltonian.

The output is electronic-only. Nuclear repulsion is recorded separately.
All integral and matrix checks are calculated before any output is written.
"""
import argparse
import json
import math
import os
import platform
import tempfile
from pathlib import Path

import numpy as np
import pyscf
import qiskit
import qiskit_nature
from pyscf import ao2mo, fci, gto, scf
from pyscf.data import nist
from qiskit_nature.second_q.drivers import InitialGuess, MethodType, PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper, ParityMapper
from qiskit_nature.units import DistanceUnit

MATRIX_TOLERANCE_HARTREE = 2.0e-12
INTEGRAL_SYMMETRY_TOLERANCE = 2.0e-12
ENERGY_CROSSCHECK_TOLERANCE_HARTREE = 2.0e-10
NUCLEAR_CROSSCHECK_TOLERANCE_HARTREE = 2.0e-12
FIXED_H2_COEFFICIENTS = {
    "II": -1.05237325,
    "IZ": 0.39793742,
    "ZI": -0.39793742,
    "ZZ": -0.01128010,
    "XX": 0.18093120,
}
FIXED_H2_ELECTRONIC_REFERENCE = -1.857275030202
FIXED_H2_NUCLEAR_REPULSION = 0.719968994449
FIXED_H2_TOTAL_REFERENCE = -1.137306035753
ANGSTROM_TO_BOHR = 1.0 / nist.BOHR
NORB_EXPECTED = 2
NMODES = 2 * NORB_EXPECTED


def ladder_matrix(terms, dimension=1 << NMODES):
    """Apply fermionic ladder strings to occupation-bit determinants."""
    matrix = np.zeros((dimension, dimension), dtype=np.complex128)
    for label, coefficient in terms:
        tokens = [] if not label else label.split()
        for ket in range(dimension):
            state = ket
            amplitude = complex(coefficient)
            for token in reversed(tokens):
                action, index_text = token.split("_")
                index = int(index_text)
                occupied = (state >> index) & 1
                if (action == "+" and occupied) or (action == "-" and not occupied):
                    amplitude = 0.0
                    break
                if (state & ((1 << index) - 1)).bit_count() & 1:
                    amplitude = -amplitude
                state ^= 1 << index
            matrix[state, ket] += amplitude
    return matrix


def integral_ladder_terms(h1, eri):
    """Expand real chemist-ordered MO integrals into spin-orbital terms.

    Spin orbitals are alpha spatial orbitals first, then beta spatial orbitals.
    (pq|rs) is chemist notation; the 1/2 factor is explicit.
    """
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


def parity_index(occupation, nmodes=NMODES):
    """Map occupation bits to cumulative parity bits, with mode 0 as LSB."""
    running = 0
    parity = 0
    for mode in range(nmodes):
        running ^= (occupation >> mode) & 1
        parity |= running << mode
    return parity


def direct_parity_matrix(determinant_matrix):
    """Permute the independently built determinant matrix into parity basis."""
    dimension = determinant_matrix.shape[0]
    result = np.zeros_like(determinant_matrix)
    for ket in range(dimension):
        for bra in range(dimension):
            result[parity_index(bra), parity_index(ket)] = determinant_matrix[bra, ket]
    return result


def reduced_parity_matrix(parity_matrix):
    """Project Nalpha=Nbeta=1 parity bits (p1=1,p3=0) to two qubits.

    The surviving QCompile bits are q0=p0 and q1=p2. Reduced basis index
    is q0 + 2*q1, preserving q0 as the least-significant bit.
    """
    def full_index(reduced_index):
        p0 = reduced_index & 1
        p2 = (reduced_index >> 1) & 1
        return p0 | (1 << 1) | (p2 << 2)

    indices = [full_index(index) for index in range(4)]
    return parity_matrix[np.ix_(indices, indices)]


def pauli_label_matrix(label):
    """Build Qiskit q1..q0 label matrix in q0-LSB computational order."""
    paulis = {
        "I": np.array([[1, 0], [0, 1]], dtype=np.complex128),
        "X": np.array([[0, 1], [1, 0]], dtype=np.complex128),
        "Y": np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
        "Z": np.array([[1, 0], [0, -1]], dtype=np.complex128),
    }
    if len(label) != 2 or any(char not in paulis for char in label):
        raise ValueError(f"unexpected two-qubit Pauli label: {label!r}")
    # Kronecker order is q1 ⊗ q0; the rightmost label character is q0.
    return np.kron(paulis[label[0]], paulis[label[1]])


def pauli_sum_matrix(terms):
    result = np.zeros((4, 4), dtype=np.complex128)
    for label, coefficient in terms:
        result += float(complex(coefficient).real) * pauli_label_matrix(label)
    return result


def maximum_absolute_difference(left, right):
    return float(np.max(np.abs(np.asarray(left) - np.asarray(right))))


def json_coefficient(value):
    z = complex(value)
    if abs(z.imag) > 1.0e-12:
        raise ValueError(f"Hamiltonian coefficient is not real: {value!r}")
    if not math.isfinite(z.real):
        raise ValueError(f"Hamiltonian coefficient is not finite: {value!r}")
    return float(z.real)


def check_integrals(h1, eri):
    checks = {
        "h1_max_abs_asymmetry_hartree": maximum_absolute_difference(h1, h1.T),
        "eri_pqrs_vs_qprs_max_abs_hartree": maximum_absolute_difference(eri, eri.transpose(1, 0, 2, 3)),
        "eri_pqrs_vs_pqsr_max_abs_hartree": maximum_absolute_difference(eri, eri.transpose(0, 1, 3, 2)),
        "eri_pqrs_vs_rspq_max_abs_hartree": maximum_absolute_difference(eri, eri.transpose(2, 3, 0, 1)),
        "h1_all_finite": bool(np.isfinite(h1).all()),
        "eri_all_finite": bool(np.isfinite(eri).all()),
        "symmetry_tolerance_hartree": INTEGRAL_SYMMETRY_TOLERANCE,
    }
    if not checks["h1_all_finite"] or not checks["eri_all_finite"]:
        raise ValueError("one- or two-electron integrals contain non-finite values")
    residual = max(value for key, value in checks.items() if key.endswith("hartree") and "tolerance" not in key)
    if residual > INTEGRAL_SYMMETRY_TOLERANCE:
        raise ValueError(f"integral symmetry residual {residual:.3e} exceeds tolerance")
    return checks


def independent_singlet_energy(h1, eri, determinant_matrix):
    """Project a direct determinant matrix onto its S=0 sector, then diagonalize."""
    sector = [
        state for state in range(1 << NMODES)
        if (state & 0b0011).bit_count() == 1
        and (state & 0b1100).bit_count() == 1
    ]
    plus_terms = [(f"+_{i} -_{i + NORB_EXPECTED}", 1.0) for i in range(NORB_EXPECTED)]
    minus_terms = [(f"+_{i + NORB_EXPECTED} -_{i}", 1.0) for i in range(NORB_EXPECTED)]
    s_plus = ladder_matrix(plus_terms)
    s_minus = ladder_matrix(minus_terms)
    s2 = 0.5 * (s_plus @ s_minus + s_minus @ s_plus)
    s2_sector = s2[np.ix_(sector, sector)]
    spin_values, spin_vectors = np.linalg.eigh(s2_sector)
    singlets = spin_vectors[:, np.abs(spin_values) < 1.0e-10]
    if singlets.shape != (4, 3):
        raise ValueError(f"expected 3 singlet CSFs in the four determinant Ms=0 sector; got {singlets.shape}")
    h_sector = determinant_matrix[np.ix_(sector, sector)]
    h_singlet = singlets.conj().T @ h_sector @ singlets
    energy = float(np.linalg.eigvalsh(h_singlet)[0].real)
    return energy, {
        "ms0_determinant_count": len(sector),
        "singlet_csf_dimension": int(singlets.shape[1]),
        "ms0_s_squared_eigenvalues": [float(value) for value in spin_values],
    }


def format_qham(distance, nuclear, fci_energy, terms):
    rows = [
        "QCOMPILE-H2-PAULI 1",
        f"distance_angstrom {distance:.17g}",
        f"distance_bohr {distance * ANGSTROM_TO_BOHR:.17g}",
        "charge 0",
        "electrons 2",
        "spin 0",
        "multiplicity 1",
        "basis STO-3G",
        "energy_unit Hartree",
        "spatial_orbitals 2",
        "spin_orbitals 4",
        "qubits 2",
        "mapping PARITY",
        "particle_counts_alpha_beta 1 1",
        "spin_orbital_order alpha-block-then-beta-block",
        "parity_reduction fixed_p1=1_fixed_p3=0",
        "pauli_label_order q1q0_rightmost_is_q0",
        "qcompile_operator_order q0_to_q1_q0_is_lsb",
        "electronic_energy_convention nuclear_repulsion_excluded",
        f"nuclear_repulsion_hartree {nuclear:.17g}",
        f"singlet_fci_electronic_reference_hartree {fci_energy:.17g}",
        f"term_count {len(terms)}",
        "TERMS",
    ]
    rows.extend(f"{label} {coefficient:.17g}" for label, coefficient in terms)
    return "\n".join(rows) + "\n"


def write_new_pair_atomic(output, qham_contents, provenance):
    """Install qham and JSON with no-overwrite links; roll back partial pairs."""
    output = output.resolve()
    sidecar = Path(str(output) + ".provenance.json")
    if not output.parent.is_dir():
        raise FileNotFoundError(f"output parent directory does not exist: {output.parent}")
    if output.exists() or sidecar.exists():
        raise FileExistsError(f"refusing to overwrite output or provenance: {output}, {sidecar}")

    temporary_paths = []
    installed = []
    try:
        for target, content in (
            (output, qham_contents),
            (sidecar, json.dumps(provenance, indent=2, allow_nan=False) + "\n"),
        ):
            with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", dir=output.parent,
                prefix=f".{target.name}.", suffix=".tmp", delete=False,
            ) as stream:
                temporary = Path(stream.name)
                temporary_paths.append(temporary)
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
        for temporary, target in zip(temporary_paths, (output, sidecar)):
            os.link(temporary, target)
            installed.append(target)
        for temporary in temporary_paths:
            temporary.unlink()
    except Exception:
        for target in installed:
            target.unlink(missing_ok=True)
        for temporary in temporary_paths:
            temporary.unlink(missing_ok=True)
        raise
    return output, sidecar


def generate(distance_angstrom, output):
    if not math.isfinite(distance_angstrom) or distance_angstrom <= 0.0:
        raise ValueError("bond distance must be finite and greater than zero")
    distance_bohr = distance_angstrom * ANGSTROM_TO_BOHR
    atom = f"H 0.0 0.0 0.0; H 0.0 0.0 {distance_bohr:.17g}"
    driver = PySCFDriver(
        atom=atom,
        unit=DistanceUnit.BOHR,
        charge=0,
        spin=0,
        basis="sto3g",
        method=MethodType.RHF,
        conv_tol=1.0e-12,
        max_cycle=100,
        init_guess=InitialGuess.MINAO,
        max_memory=2000,
    )
    problem = driver.run()
    fermionic = problem.hamiltonian.second_q_op()
    mapped = ParityMapper(num_particles=(1, 1)).map(fermionic)
    mapped_terms = sorted(
        ((label, json_coefficient(coefficient)) for label, coefficient in mapped.to_list()),
        key=lambda item: item[0],
    )
    if mapped.num_qubits != 2 or not mapped_terms:
        raise ValueError(f"expected a non-empty two-qubit mapped Hamiltonian, got {mapped.num_qubits} qubits")
    if len(mapped_terms) != len({label for label, _ in mapped_terms}):
        raise ValueError("mapper emitted duplicate Pauli labels; refusing ambiguous output")

    mol = gto.M(
        atom=atom, unit="Bohr", basis="sto-3g", charge=0, spin=0,
    )
    mf = scf.RHF(mol)
    mf.conv_tol = 1.0e-12
    mf.max_cycle = 100
    mf.init_guess = "minao"
    mf.max_memory = 2000
    hf_total = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError("explicit PySCF RHF did not converge")
    if mf.mo_coeff.shape != (NORB_EXPECTED, NORB_EXPECTED):
        raise ValueError(f"expected two retained spatial orbitals; got {mf.mo_coeff.shape}")

    h1 = mf.mo_coeff.T @ mf.get_hcore() @ mf.mo_coeff
    eri = ao2mo.restore(1, ao2mo.kernel(mol, mf.mo_coeff), NORB_EXPECTED)
    integral_checks = check_integrals(h1, eri)
    nuclear = float(mol.energy_nuc())
    nuclear_qn = float(problem.nuclear_repulsion_energy)
    nuclear_qn_residual = abs(nuclear - nuclear_qn)
    if nuclear_qn_residual > NUCLEAR_CROSSCHECK_TOLERANCE_HARTREE:
        raise ValueError(f"PySCF/Qiskit Nature nuclear energy residual {nuclear_qn_residual:.3e} Ha")
    h1_qn = np.asarray(problem.hamiltonian.electronic_integrals.alpha["+-"])
    h1_qn_residual = maximum_absolute_difference(h1, h1_qn)
    if h1_qn_residual > MATRIX_TOLERANCE_HARTREE:
        raise ValueError(f"PySCF/Qiskit Nature one-electron integral residual {h1_qn_residual:.3e} Ha")

    integral_terms = integral_ladder_terms(h1, eri)
    determinant_matrix = ladder_matrix(integral_terms)
    # FermionicOp has no dense-matrix API in Qiskit Nature 0.8; map its ladder
    # terms through Jordan-Wigner and compare in occupation-basis ordering.
    qn_fermionic_matrix = JordanWignerMapper().map(fermionic).to_matrix()
    fermionic_residual = maximum_absolute_difference(determinant_matrix, qn_fermionic_matrix)
    if fermionic_residual > MATRIX_TOLERANCE_HARTREE:
        raise ValueError(f"integral-built/Qiskit Nature fermionic residual {fermionic_residual:.3e} Ha")

    full_parity = ParityMapper().map(fermionic)
    full_parity_matrix = full_parity.to_matrix()
    direct_parity = direct_parity_matrix(determinant_matrix)
    full_parity_residual = maximum_absolute_difference(direct_parity, full_parity_matrix)
    reduced_direct = reduced_parity_matrix(direct_parity)
    reduced_pauli = pauli_sum_matrix(mapped_terms)
    reduced_qiskit_residual = maximum_absolute_difference(reduced_direct, mapped.to_matrix())
    reduced_pauli_residual = maximum_absolute_difference(reduced_direct, reduced_pauli)
    matrix_residual = max(
        fermionic_residual, full_parity_residual,
        reduced_qiskit_residual, reduced_pauli_residual,
    )
    if matrix_residual > MATRIX_TOLERANCE_HARTREE:
        raise ValueError(f"maximum independent matrix residual {matrix_residual:.3e} Ha")

    singlet_energy, sector_checks = independent_singlet_energy(h1, eri, determinant_matrix)
    fci_solver = fci.direct_spin0.FCI(mol)
    fci_solver.conv_tol = 1.0e-12
    fci_solver.max_cycle = 100
    pyscf_fci_energy, _ = fci_solver.kernel(
        h1, eri, NORB_EXPECTED, (1, 1), ecore=0.0,
    )
    pyscf_fci_energy = float(pyscf_fci_energy)
    fci_residual = abs(singlet_energy - pyscf_fci_energy)
    if fci_residual > ENERGY_CROSSCHECK_TOLERANCE_HARTREE:
        raise ValueError(f"direct singlet projection/PySCF FCI residual {fci_residual:.3e} Ha")

    fixed_comparison = None
    if abs(distance_angstrom - 0.735) <= 1.0e-12:
        observed = dict(mapped_terms)
        labels = set(observed) | set(FIXED_H2_COEFFICIENTS)
        coefficient_residual = max(
            abs(observed.get(label, 0.0) - FIXED_H2_COEFFICIENTS.get(label, 0.0))
            for label in labels
        )
        fixed_matrix = pauli_sum_matrix(sorted(FIXED_H2_COEFFICIENTS.items()))
        fixed_electronic = float(np.linalg.eigvalsh(fixed_matrix)[0].real)
        energy_residual = abs(singlet_energy - FIXED_H2_ELECTRONIC_REFERENCE)
        nuclear_residual = abs(nuclear - FIXED_H2_NUCLEAR_REPULSION)
        total_energy = singlet_energy + nuclear
        total_residual = abs(total_energy - FIXED_H2_TOTAL_REFERENCE)
        fixed_comparison = {
            "coefficient_max_abs_difference_hartree": coefficient_residual,
            "coefficient_tolerance_hartree": 1.0e-8,
            "singlet_energy_minus_documented_electronic_reference_hartree": singlet_energy - FIXED_H2_ELECTRONIC_REFERENCE,
            "singlet_energy_reference_abs_difference_hartree": energy_residual,
            "energy_tolerance_hartree": 1.0e-8,
            "fixed_rounded_hamiltonian_lowest_eigenvalue_hartree": fixed_electronic,
            "nuclear_repulsion_abs_difference_hartree": nuclear_residual,
            "nuclear_repulsion_tolerance_hartree": NUCLEAR_CROSSCHECK_TOLERANCE_HARTREE,
            "fixed_nuclear_repulsion_hartree": FIXED_H2_NUCLEAR_REPULSION,
            "total_energy_hartree": total_energy,
            "fixed_total_reference_hartree": FIXED_H2_TOTAL_REFERENCE,
            "total_energy_abs_difference_hartree": total_residual,
            "total_energy_tolerance_hartree": 1.0e-8,
        }
        if (coefficient_residual > 1.0e-8 or energy_residual > 1.0e-8
                or nuclear_residual > NUCLEAR_CROSSCHECK_TOLERANCE_HARTREE
                or total_residual > 1.0e-8):
            raise ValueError(f"generated H2 differs from fixed reference: {fixed_comparison}")

    provenance = {
        "schema": "QCOMPILE-H2-PROVENANCE-1",
        "model": {
            "system": "H2", "geometry_angstrom": [["H", 0.0, 0.0, 0.0], ["H", 0.0, 0.0, distance_angstrom]],
            "geometry_bohr": [["H", 0.0, 0.0, 0.0], ["H", 0.0, 0.0, distance_bohr]],
            "geometry_units": "angstrom input; bohr passed explicitly to both PySCF paths",
            "distance_angstrom": distance_angstrom, "distance_bohr": distance_bohr,
            "basis": "STO-3G", "integral_method": "RHF", "frozen_orbitals": 0,
            "spatial_orbitals": NORB_EXPECTED, "spin_orbitals": NMODES,
            "spin_orbital_order": ["MO0 alpha", "MO1 alpha", "MO0 beta", "MO1 beta"],
            "mapping": "ParityMapper(num_particles=(1,1))",
            "parity_mapping_definition": "p_q = XOR of occupation bits n_0..n_q",
            "fixed_parities": {"p1_alpha_parity": 1, "p3_total_parity": 0},
            "reduced_qubits": {"q0": "p0", "q1": "p2"},
            "state_index_convention": "occupation/parity bit q is 2^q; q0 least significant",
            "qiskit_pauli_label_convention": "left-to-right q1q0; rightmost character acts on q0",
            "energy_convention": "Pauli Hamiltonian is electronic only; nuclear repulsion is separate",
        },
        "settings": {
            "pyscf_rhf_conv_tol_hartree": 1.0e-12, "pyscf_rhf_max_cycle": 100,
            "pyscf_rhf_initial_guess": "MINAO", "pyscf_max_memory_mb": 2000,
            "pyscf_fci_solver": "pyscf.fci.direct_spin0.FCI", "pyscf_fci_conv_tol": 1.0e-12,
            "pyscf_fci_max_cycle": 100, "pyscf_fci_nelec_alpha_beta": [1, 1],
            "frozen_core": False, "integral_order": "chemist (pq|rs)",
            "precision": "IEEE-754 binary64; JSON round-trip decimal serialization; qham coefficients use 17 significant digits",
            "matrix_abs_tolerance_hartree": MATRIX_TOLERANCE_HARTREE,
            "integral_symmetry_abs_tolerance_hartree": INTEGRAL_SYMMETRY_TOLERANCE,
            "energy_crosscheck_tolerance_hartree": ENERGY_CROSSCHECK_TOLERANCE_HARTREE,
            "nuclear_crosscheck_tolerance_hartree": NUCLEAR_CROSSCHECK_TOLERANCE_HARTREE,
        },
        "versions": {
            "python": platform.python_version(), "numpy": np.__version__,
            "pyscf": pyscf.__version__, "qiskit": qiskit.__version__,
            "qiskit_nature": qiskit_nature.__version__,
        },
        "scf": {
            "converged": bool(mf.converged), "cycles": int(mf.cycles),
            "total_including_nuclear_hartree": hf_total,
            "electronic_hartree": hf_total - nuclear,
            "nuclear_repulsion_analytic_1_over_r_hartree": 1.0 / distance_bohr,
            "mo_coefficient_shape": list(mf.mo_coeff.shape),
        },
        "integrals_hartree": {
            "one_electron_mo_h_pq": h1.tolist(),
            "two_electron_mo_chemist_eri_pqrs": eri.tolist(),
            "checks": integral_checks,
            "qiskit_nature_h1_max_abs_difference_hartree": h1_qn_residual,
            "qiskit_nature_nuclear_max_abs_difference_hartree": nuclear_qn_residual,
        },
        "fermionic_hamiltonian": {
            "terms": [
                {"ladder_string": label, "coefficient": json_coefficient(coefficient)}
                for label, coefficient in fermionic.items()
            ],
            "term_count": len(list(fermionic.items())),
            "matrix_energy_convention": "electronic only",
        },
        "pauli_hamiltonian": {
            "terms": [{"label_q1q0": label, "coefficient_hartree": coefficient,
                       "qcompile_operators_q0_to_q1": list(reversed(label))}
                      for label, coefficient in mapped_terms],
            "term_count": len(mapped_terms),
        },
        "energies_hartree": {
            "electronic_singlet_direct_determinant_projection": singlet_energy,
            "electronic_singlet_pyscf_fci": pyscf_fci_energy,
            "nuclear_repulsion": nuclear,
            "total_born_oppenheimer_direct_singlet": singlet_energy + nuclear,
            "total_born_oppenheimer_pyscf_fci": pyscf_fci_energy + nuclear,
        },
        "physical_sector": sector_checks,
        "validation": {
            "method": "independent determinant matrix from PySCF MO h1 and chemist ERIs; ladder operators applied directly to occupation determinants; parity permutation then explicit p1=1,p3=0 projection; Pauli matrices built independently with q1 tensor q0 order",
            "integral_determinant_vs_qn_jordan_wigner_matrix_max_abs_hartree": fermionic_residual,
            "full_parity_determinant_vs_qiskit_mapper_max_abs_hartree": full_parity_residual,
            "reduced_determinant_vs_qiskit_operator_max_abs_hartree": reduced_qiskit_residual,
            "reduced_determinant_vs_independent_pauli_sum_max_abs_hartree": reduced_pauli_residual,
            "maximum_matrix_residual_hartree": matrix_residual,
            "required_maximum_matrix_residual_hartree": MATRIX_TOLERANCE_HARTREE,
            "pyscf_direct_spin0_vs_determinant_singlet_abs_difference_hartree": fci_residual,
            "pyscf_direct_spin0_tolerance_hartree": ENERGY_CROSSCHECK_TOLERANCE_HARTREE,
            "fixed_h2_comparison_at_0_735_angstrom": fixed_comparison,
            "labels_explicitly_reversed_for_qcompile_q0_lsb": True,
        },
    }
    qham_contents = format_qham(distance_angstrom, nuclear, pyscf_fci_energy, mapped_terms)
    return write_new_pair_atomic(Path(output), qham_contents, provenance), provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance-angstrom", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        (qham_path, provenance_path), provenance = generate(args.distance_angstrom, args.output)
    except Exception as error:
        parser.exit(2, f"generation failed: {error}\n")
    print(json.dumps({
        "qham": str(qham_path),
        "provenance": str(provenance_path),
        "distance_angstrom": args.distance_angstrom,
        "pauli_terms": provenance["pauli_hamiltonian"]["terms"],
        "energies_hartree": provenance["energies_hartree"],
        "nuclear_repulsion_hartree": provenance["energies_hartree"]["nuclear_repulsion"],
        "validation": provenance["validation"],
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
