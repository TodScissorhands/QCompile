# QCompile — Milestone 1 Project Report

## Purpose and scope

QCompile is an educational C++17 compiler for a small quantum-circuit DSL and a classical state-vector simulator. Milestone 1 adds a fixed, exact-state-vector evaluation of the H₂ electronic ground-state energy using an established VQE-style variational calculation. It does not add general molecular-integral generation, geometry-to-Hamiltonian generation, support for additional molecules, shot-based measurement, or a graphical interface.

This is an engineering and educational implementation of established methods—Pauli Hamiltonians, state-vector expectation values, variational optimization, and exact diagonalization. It makes no claim of a novel quantum algorithm.

## Existing compiler and simulator

The compiler pipeline is lexer → recursive-descent parser → AST → semantic checker → instruction generator → simulator. The DSL supports qubit and variable declarations, arithmetic expressions for numeric parameters, `H`, `X`, `Z`, `CNOT`, `RZ`, and `measure`. Instructions retain execution order and use zero-based qubit indices.

The simulator stores complex amplitudes in a state vector, applies gates to those amplitudes, and reports basis-state probabilities. Measurement samples an outcome using the existing fixed seed (`12345`) and prints the requested bits; it does **not** collapse or otherwise update the state. The 24-qubit safety limit remains in place. The original compiler path and measurement behavior remain available and are covered by regression tests.

## Fixed H₂ model

The benchmark is neutral singlet H₂ at a 0.735 Å internuclear distance, in the STO-3G minimal basis, with two electrons. Qiskit Nature’s parity mapping with particle counts `(1, 1)` reduces the electronic Hamiltonian to two qubits. In QCompile’s convention, q0 is the least-significant state-vector bit.

The fixed electronic Hamiltonian is

```text
-1.05237325 I
+0.39793742 Z0
-0.39793742 Z1
-0.01128010 Z0 Z1
+0.18093120 X0 X1
```

Coefficients and energies are in Hartree. Qiskit Nature reports an electronic reference of approximately `-1.857275030202 Ha` and nuclear repulsion of `+0.719968994449 Ha`, giving a total Born–Oppenheimer reference of `-1.137306035753 Ha`. Electronic energy excludes nuclear repulsion; the total adds the nuclear term once.

The rounded Hamiltonian’s independently diagonalized minimum is `-1.8572750292882292 Ha`. Its difference from the cited electronic reference is approximately `9.14e-10 Ha`, consistent with coefficient rounding. The corresponding total from the rounded Hamiltonian is `-1.1373060348392292 Ha`.

Sources:

- [Qiskit Nature: Electronic structure](https://qiskit-community.github.io/qiskit-nature/tutorials/01_electronic_structure.html)
- [Qiskit Nature: Qubit mappers](https://qiskit-community.github.io/qiskit-nature/tutorials/06_qubit_mappers.html)

## Pauli Hamiltonian and expectation evaluation

`include/chemistry/hamiltonian.h` and `src/chemistry/hamiltonian.cpp` define I, X, Y, and Z Pauli operators, weighted Pauli strings, and a weighted-sum `QubitHamiltonian`. `pauliExpectation` applies each operator to complex amplitudes using bit flips and phases, then evaluates the inner product. It correctly handles X-containing terms; it does not use measurement probabilities as a substitute. The implementation checks state dimension, finite amplitudes, normalization, and the real-valued result expected for a Hermitian Pauli product.

`Simulator::expectation` evaluates against the simulator’s current state without exposing the state vector or changing the existing measurement path. Hamiltonian construction validates qubit count, coefficient finiteness, and Pauli-string width.

## Trial state, optimization, and exact reference

`src/chemistry/h2_energy.cpp` creates a fresh two-qubit simulator and instruction list for each parameter value. The circuit starts in the encoded reference state and rotates within the odd-parity two-state sector. Its parameterized rotation is synthesized from the existing H and RZ gates; no new gate family or compiler-language feature is introduced. The family can represent the ground state of the benchmark Hamiltonian in that sector.

A deterministic golden-section search minimizes the electronic expectation over `[0, π/2]`, stops when the bracket width is at most `1e-12` radians, and has a 128-iteration limit. The final energy is checked against an independently constructed 4×4 Hamiltonian matrix. That matrix is diagonalized using Jacobi rotations, not the VQE optimizer. The optimizer also checks its rounded-Hamiltonian minimum against the documented electronic reference within `1e-8 Ha` and rejects a final variational energy more than `1e-10 Ha` above the exact minimum.

The run reported parameter `0.223536992823` radians after 62 energy evaluations. The optimized electronic energy and independent exact minimum both print as `-1.857275029288 Ha`; nuclear repulsion is `+0.719968994449 Ha`; total energy is `-1.137306034839 Ha`. The parameter is within the test’s `1e-7` tolerance of the analytic ground-state mixing angle. The total differs from the cited rounded total reference by about `9.14e-10 Ha`.

## Verification

Commands executed from the repository root; each completed with exit code 0:

| Command | Observed result |
|---|---|
| `make all` | Built `qcompile` and `qcompile_h2` with `-std=c++17 -Wall -Wextra -pedantic`; no compiler diagnostics. |
| `make test` | `Passed 78 assertions.` |
| `./qcompile_h2` | Printed the H₂ model, optimized and exact energies, nuclear repulsion, total energy, and 62 evaluations as stated above. |
| `./qcompile samples/bell_state.txt` | Compiler pipeline completed; probabilities were `00: 0.500000` and `11: 0.500000`; seeded sample was `11`. |
| `git diff --check` | No output; no whitespace errors reported. |

Before building, `qcompile`, `qcompile_h2`, and `qcompile_tests` were absent. The Makefile targets therefore did not overwrite pre-existing executables. `make all` created the first two binaries; `make test` created `qcompile_tests`. No `make clean` was run. The existing `QCompile Project Report.pdf` was inspected and left unchanged; it is an earlier academic report describing the compiler/simulator, not this H₂ milestone.

The regression suite covers Pauli expectations on known and superposition states, multi-qubit products, direct dense-matrix comparisons, invalid state inputs, parameter-sensitive energy evaluation, deterministic optimization, the rounded-Hamiltonian reference, nuclear-repulsion separation, and existing compiler, simulator, measurement, and capacity behavior.

## Limitations and future work

The molecule-specific Hamiltonian, ansatz, optimizer, and 4×4 exact reference routine are fixed to this H₂ benchmark. The simulator remains an ideal state-vector simulator; measurement is terminal-style sampling without collapse. The result does not demonstrate hardware execution or a general molecular chemistry pipeline.

Future work, not implemented in this milestone:

1. Define and independently validate a separate fixed HeH⁺ benchmark before adding it; specify geometry, basis, electron count, charge, spin, mapping, qubit order, Hamiltonian coefficients, and energy conventions.
2. Build a small library of documented fixed Hamiltonians with independent exact references.
3. Separately scope geometry-derived Hamiltonian generation after the fixed-library path is validated.
4. Add educational visualizations only when tied to verified calculations.

No HeH⁺ support, arbitrary molecular support, geometry generator, or GUI is included here.
