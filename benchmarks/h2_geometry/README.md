# Geometry-specific H₂ reference generation

This directory contains an **offline reference-data generator**, not a C++ runtime feature. It generates a fresh STO-3G Hamiltonian for neutral singlet H₂ from the requested bond length. C++ parsing and benchmark integration are out of scope for this step.

## Fixed model and conventions

- H₂, charge 0, two electrons, singlet `S=0` / multiplicity 1.
- Geometry input: `H (0,0,0); H (0,0,R)` with `R` in Ångström. The generator converts Å to bohr using the pinned PySCF `nist.BOHR=0.52917721092 Å/bohr` value, then passes the identical explicit bohr coordinates to Qiskit Nature and direct PySCF. The Psi4 validation script uses the same conversion, avoiding package-specific Å-to-bohr rounding differences.
- STO-3G, RHF molecular orbitals, all two spatial orbitals retained; no frozen orbitals.
- PySCF and Qiskit Nature generate MO integrals and the electronic Hamiltonian.
- Spin-orbital order: MO0 alpha, MO1 alpha, MO0 beta, MO1 beta.
- Qiskit Nature parity mapping with `num_particles=(1,1)`. The parity bits fixed by that particle count are `p1=1` and `p3=0`; the surviving QCompile qubits are `q0=p0` and `q1=p2`.
- Occupation/parity bit `q` has weight `2^q`; q0 is the least-significant bit. Serialized Pauli labels use Qiskit's left-to-right `q1q0` convention, so the rightmost character is q0. Provenance records the reversed QCompile `q0`-to-`q1` operator sequence.
- The Pauli Hamiltonian contains **electronic energy only**. Nuclear repulsion is a separate field and is added once to calculate total Born–Oppenheimer energy.
- One-electron MO integrals are `h[p,q]`; two-electron integrals are real chemist-ordered `(pq|rs)`, with the `1/2` factor explicit in the spin-orbital Hamiltonian.

The generator rejects invalid distances, non-converged RHF, unexpected orbital or qubit counts, non-finite integrals or coefficients, failed integral/matrix checks, and incompatible FCI results. It writes full round-trip binary64 values. It refuses to overwrite either the requested output or its `.provenance.json` sidecar. Both files are staged and installed without replacement only after the calculations and validation succeed.

## Environment

Use the existing generation-only environment and pinned package set from `benchmarks/heh_plus/requirements-generator.lock`. The currently inspected environment ran the generator stack with Python 3.11.9, NumPy 2.4.6, PySCF 2.8.0, Qiskit 2.5.2, and Qiskit Nature 0.8.0. These packages are not C++ build or runtime dependencies.

The environment was smoke-tested by importing the packages, running `PySCFDriver` for H₂ at 0.735 Å, and mapping to two qubits. For clean environment recreation, follow the pinned environment instructions in `benchmarks/heh_plus/GATE_A_STATUS.md`; do not install these packages system-wide.

For a clean recreation, use the recorded `uv 0.13.1` toolchain and synchronize the exact lock:

```bash
uv python install 3.11.9
uv venv --python 3.11.9 benchmarks/heh_plus/.venv-generator
uv pip sync --python benchmarks/heh_plus/.venv-generator/bin/python \
  benchmarks/heh_plus/requirements-generator.lock
```


## Generate data

Run from the repository root. Output parent directories must already exist; the examples use `/tmp` to avoid adding generated artifacts to the repository.

```bash
benchmarks/heh_plus/.venv-generator/bin/python -B \
  benchmarks/h2_geometry/generate_h2.py \
  --distance-angstrom 0.735 \
  --output /tmp/h2-0735.qham

benchmarks/heh_plus/.venv-generator/bin/python -B \
  benchmarks/h2_geometry/generate_h2.py \
  --distance-angstrom 1.400 \
  --output /tmp/h2-1400.qham
```

Each command emits the versioned `QCOMPILE-H2-PAULI 1` file, a sibling `.qham.provenance.json` sidecar, and a JSON summary to stdout. The sidecar includes package versions, RHF settings/convergence, complete one- and two-electron MO integrals, fermionic terms, Pauli terms, singlet-sector information, electronic/nuclear/total energies, and the numerical checks.

The generator calculates each requested geometry; no coefficients or Hamiltonians are selected from a geometry lookup table. The two prescribed distances are validation points, not a claim that every positive distance is chemically or numerically supported.

## Validation performed by the generator

- Finite integrals; one-electron symmetry and the three real two-electron permutation symmetries, each checked against `2e-12 Ha`.
- A custom determinant-space matrix built by applying fermionic ladder operators to all 16 occupation states from the direct PySCF MO integrals. It is compared with the occupation-basis matrix obtained by Jordan–Wigner mapping the Qiskit Nature `FermionicOp`.
- Independent parity-bit permutation of that determinant matrix, projection onto `p1=1,p3=0`, and comparison of the resulting 4×4 matrix with both the reduced Qiskit Pauli operator and a separately constructed Pauli matrix. Pauli matrices are assembled as `q1 ⊗ q0`; label reversal is explicitly tested in the reported residual.
- Every matrix-comparison maximum must be at most `2e-12 Ha`.
- The determinant-space `Nα=Nβ=1` sector has four determinants. The generator constructs `S²`, retains its three `S=0` configuration-state functions, and compares the lowest projected eigenvalue to PySCF `direct_spin0` FCI. This energy comparison uses `2e-10 Ha` as a fixed solver cross-check tolerance.
- At 0.735 Å only, generated coefficients are compared with the existing rounded H₂ Hamiltonian within `1e-8 Ha` per coefficient; electronic and total energies are compared with their documented references within `1e-8 Ha`, and nuclear repulsion within `1e-12 Ha`. The nuclear constant is not part of the Pauli operator.

PySCF FCI and the direct determinant-space diagonalization use the same generated integrals, so they are complementary internal checks, not independent chemistry packages. `validate_psi4_h2.py` performs a separate Psi4 DETCI FCI calculation and compares electronic, nuclear, and total energies against the matching generator provenance.

## Independent Psi4 calculation

Use the compatible, isolated Psi4 environment. The output JSON, log, and scratch directory must all be new. Run from a temporary working directory so Psi4's process-local `timer.dat` cannot touch repository timer artifacts. Example (adjust only the `/tmp` names to unused paths):

```bash
mkdir /tmp/qcompile-h2-psi4-run
cd /tmp/qcompile-h2-psi4-run
/home/tod/QCompile/benchmarks/heh_plus/.conda-psi4-libint290-attempt1/bin/python -B \
  /home/tod/QCompile/benchmarks/h2_geometry/validate_psi4_h2.py \
  --distance-angstrom 0.735 \
  --generator-provenance /tmp/h2-0735-final.qham.provenance.json \
  --output-json /tmp/h2-psi4-0735.json \
  --scratch /tmp/h2-psi4-scratch-0735
```

Repeat for `1.400` with the matching `/tmp/h2-1400-final.qham.provenance.json` and distinct JSON/scratch paths. This script runs RHF/FCI in STO-3G, charge 0, `Nα=Nβ=1`, with no frozen orbitals. Psi4 reports total energy including nuclear repulsion; the script records electronic energy as total minus nuclear repulsion. It obtains `<S²>` from Psi4's FCI output log if the wavefunction API does not expose it, and rejects a non-singlet result.

The Psi4 check is a separate implementation and is the cross-package reference. If the installed environment crashes or fails to converge, record that actual failure; do not treat PySCF's own reference as Psi4 validation.

## Scope limits

At the generator-only Milestone 5A boundary, this directory did not include C++ file parsing, geometry optimization, or visualization. Milestone 5B later added the restricted C++ loader and execution path; Milestone 6 adds a native H₂ laboratory that invokes this generator. Arbitrary molecules, bases, active spaces, a general mapper, geometry optimization, and density visualization remain out of scope.
