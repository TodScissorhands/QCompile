# HeH+ Reference Data — Gate A Status

**Status: ACCEPTED.** The fixed STO-3G HeH+ model has an independently reproduced Psi4 DETCI FCI reference. The earlier Psi4 environment failure is retained below as historical evidence; it is superseded by the successful isolated Libint2 environment run.

## Fixed model

- HeH+ at He `(0,0,0)` and H `(0,0,1.463283)` bohr; this is the high-accuracy Born–Oppenheimer equilibrium geometry reported by [Pachucki (2012)](https://arxiv.org/abs/1203.3927), not a claim that it is the STO-3G equilibrium geometry. That paper's high-accuracy total energy (`-2.978708310771 Ha`) is not a STO-3G reference; the model comparisons below use the same STO-3G Hamiltonian in PySCF/Qiskit Nature and Psi4.
- STO-3G, charge `+1`, two electrons, multiplicity `1`, `Nα=Nβ=1`, RHF orbitals, two spatial orbitals, no frozen orbitals.
- Four spin orbitals/qubits, Jordan–Wigner, no symmetry tapering. Mode order is `q0=MO0 α`, `q1=MO1 α`, `q2=MO0 β`, `q3=MO1 β`; q0 is the least-significant occupation bit. Qiskit labels are printed `q3…q0`.
- The qubit Hamiltonian and electronic FCI energy exclude nuclear repulsion. Total Born–Oppenheimer energy is electronic energy plus `1.3667896093920315 Ha` exactly once.

## Independent result and validation

The successful cross-package calculation used Psi4 1.9.1 with Python 3.11.17 from the isolated `benchmarks/heh_plus/.conda-psi4-libint290-attempt1` environment. It ran one-thread RHF/PK and DETCI FCI in STO-3G, with `1e-12` SCF energy/density thresholds, `1e-10` FCI residual threshold, one root, and no frozen orbitals. The Psi4 output log reports `<S²> = 0.0`; the retained JSON now records this value as well.

| Quantity | PySCF/Qiskit Nature reference | Independent Psi4 result | Absolute difference |
|---|---:|---:|---:|
| Electronic FCI energy (Ha) | `-4.218264154490395` | `-4.218264154490385` | `1.0e-14` |
| Nuclear repulsion (Ha) | `+1.3667896093920315` | `+1.3667896093920315` | `0` |
| Total Born–Oppenheimer energy (Ha) | `-2.8514745450983634` | `-2.8514745450983536` | `9.8e-15` |

The PySCF-side independently assembled 16×16 determinant matrix agrees with the Qiskit Nature fermionic and mapped Pauli matrices at maximum residuals `1.11e-16 Ha` and `4.88e-15 Ha`, respectively. Its singlet-sector eigenvalue agrees with Psi4 FCI within `1.0e-14 Ha`. The C++ fixed HeH+ optimizer reports `-4.218264154490334 Ha` electronic energy, `-4.218264154490395 Ha` exact singlet energy, `+1.366789609392032 Ha` nuclear repulsion, and `-2.851474545098302 Ha` total energy; its exact electronic result differs from Psi4 by `1.0e-14 Ha`.

## Reproduction

Run these commands from the repository root. The checked-in retry script writes outputs alongside itself and requires a new scratch directory. Copy it to a fresh temporary directory to avoid overwriting retained evidence:

```bash
repo=$(pwd)
workdir=$(mktemp -d /tmp/qcompile-heh-psi4.XXXXXX)
cp "$repo/benchmarks/heh_plus/run_psi4_fci_libint290_attempt1_retry.py" "$workdir/run_psi4_fci.py"
(cd "$workdir" && "$repo/benchmarks/heh_plus/.conda-psi4-libint290-attempt1/bin/python" -B "$workdir/run_psi4_fci.py")
```

The script emits JSON and a Psi4 log in `workdir`; the JSON contains the electronic, nuclear, and total energies, model settings, and spin diagnostic. To reproduce the PySCF/Qiskit Nature reference without overwriting tracked data, copy the generator into a fresh directory and run it there:

```bash
repo=$(pwd)
workdir=$(mktemp -d /tmp/qcompile-heh-reference.XXXXXX)
cp "$repo/benchmarks/heh_plus/generate_reference.py" "$workdir/generate_reference.py"
(cd "$workdir" && "$repo/benchmarks/heh_plus/.venv-generator/bin/python" -B "$workdir/generate_reference.py")
```

## Historical failed attempt

An earlier Psi4 1.9.1 environment (`benchmarks/heh_plus/.conda-psi4`) exited with signal 11 (139) during basis loading. Switching to the separately installed Libint2 environment resolved that runtime failure. The earlier failure remains relevant provenance, but it is no longer a gate blocker; the successful output and JSON are retained as `psi4_fci_output_libint290_attempt1_retry.txt` and `psi4_fci_reference_libint290_attempt1_retry.json`.

## Limitations

This cross-check validates one fixed geometry and one frozen model. It does not establish an STO-3G optimized bond length, validate other geometries, or establish general-molecule support. The selected geometry is a high-accuracy physical Born–Oppenheimer equilibrium distance, not the minimum of the STO-3G curve. The published high-accuracy energy is a different computational model and is not compared directly with the STO-3G results above.
