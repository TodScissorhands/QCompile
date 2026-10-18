# QCompile

QCompile is a domain-specific language (DSL) compiler and classical statevector simulator for quantum circuits, developed in C++17. It translates high-level quantum circuit definitions into simulator instructions, validates semantics and arithmetic expressions, performs adjacent-gate cancellation optimizations, and simulates state evolution and Pauli expectations on a statevector. It also provides fixed H₂ and HeH⁺ VQE benchmarks, geometry-generated H₂ and a fixed LiH active-space benchmark, an interactive Qt6 desktop workbench, and a local loopback dashboard.

---

## Architecture Overview

```text
Source File (.txt)
       │
       ▼
 ┌───────────┐
 │   Lexer   │  Tokenizes keywords, gates, identifiers, literals, operators
 └─────┬─────┘
       │
       ▼
 ┌───────────┐
 │  Parser   │  Constructs AST with recursive-descent arithmetic expressions
 └─────┬─────┘
       │
       ▼
 ┌───────────┐
 │ Semantic  │  Validates qubit declarations, gate arities, unique CNOT pairs,
 │  Checker  │  and statically evaluates parameter expressions
 └─────┬─────┘
       │
       ▼
 ┌───────────┐
 │   Code    │  Emits flat simulator instructions with integer qubit indices
 │ Generator │
 └─────┬─────┘
       │
       ▼
 ┌───────────┐
 │ Optimizer │  Cancels adjacent identical self-inverse gates (H, X, Z, CNOT)
 └─────┬─────┘
       │
       ▼
 ┌───────────┐
 │ Simulator │  Evolves statevector ($2^n$ complex amplitudes), computes Pauli
 └───────────┘  expectations, and samples measurement outcomes
```

### Subsystems and Directory Structure

| Subsystem | Headers | Implementation | Purpose |
| :--- | :--- | :--- | :--- |
| **Compiler** | `include/compiler/` | `src/compiler/` | Lexer, AST parser, semantic checker, code generator, and DAG-based circuit optimizer |
| **Simulator** | `include/simulator/` | `src/simulator/` | Statevector simulator, single/two-qubit gates, measurement, and Pauli expectations |
| **Chemistry** | `include/chemistry/` | `src/chemistry/` | Qubit Hamiltonians, Jacobi eigensolver, H₂ VQE, HeH⁺ VQE, and benchmark catalog |
| **GUI** | — | `src/gui/h2_lab.cpp` | Interactive Qt6 desktop workbench for H₂ potential energy surface inspection |
| **Apps** | — | `src/apps/` | Entry points: `qcompile`, `qcompile_h2`, `qcompile_heh`, `qcompile_benchmarks` |

---

## Supported Language Constructs

The QCompile DSL uses line-oriented syntax:

- **Qubit Declarations**: `qubit q0, q1, q2`
- **Variable Declarations**: `let theta = 0.5`, `let phi = 2 * 3.14159 / 4`
- **Unitary Gates**:
  - Fixed: `H <qubit>`, `X <qubit>`, `Z <qubit>`, `CNOT <control>, <target>`
  - Parameterized: `RZ(<expression>) <qubit>`
- **Arithmetic Expressions**: Supports literals, identifiers, unary `-`, binary `+`, `-`, `*`, `/`, and parentheses with standard operator precedence.
- **Measurements**: `measure <qubit> [, <qubit>...]`

Example Bell-state circuit (`samples/bell_state.txt`):
```text
qubit q0, q1
H q0
CNOT q0, q1
measure q0, q1
```

---

## Implementation

The implementation provides command-line compiler/simulator tools alongside visual workbench interfaces:

| Terminal Simulation Output | Interactive Qt6 H₂ Workbench |
| :---: | :---: |
| ![Bell-state compiler and simulator output](docs/images/bell_state_output.png) | ![Qt6 H2 Workbench displaying potential energy curve and circuit inspection](docs/images/h2_workbench.png) |
| *Compilation trace, optimization pass, and statevector measurement output for a 2-qubit Bell state.* | *Native Qt6 H₂ workbench displaying potential energy surface scan and statevector inspection.* |

---

## Prerequisites and Build Instructions

### Requirements
- **C++17 Compiler**: `g++` (GCC $\ge$ 9) or `clang++`
- **Build System**: GNU `make`
- **Qt6 Widgets** *(optional, for `qcompile_lab`)*: `qt6-base-dev` and `pkg-config`
- **Node.js** *(optional, for loopback dashboard)*: Node.js $\ge$ 18

### Building

From the repository root:

```bash
# Build compiler, simulator, and benchmark CLIs (qcompile, qcompile_h2, qcompile_heh, qcompile_benchmarks)
make all

# Build the Qt6 workbench (requires Qt6 development libraries)
make qcompile_lab

# Build the regression test suite
make qcompile_tests

# Remove build artifacts
make clean
```

---

## Execution Guide

### 1. Quantum Circuit DSL Compiler
Compile and simulate a DSL circuit:
```bash
./qcompile samples/bell_state.txt
```
Run other provided sample programs:
```bash
./qcompile samples/single_gate.txt
./qcompile samples/sequential_gates.txt
./qcompile samples/three_qubit.txt
./qcompile samples/parameterized_gate.txt
./qcompile samples/precedence_check.txt
```

### 2. Regression Test Suite
Run the complete regression suite:
```bash
make test
# or directly:
./qcompile_tests
```
The suite covers the compiler, optimizer, simulator, workbench integration, molecular integrals, and VQE.

Test modules under `tests/` include `test_lexer.cpp`, `test_parser.cpp`, `test_semantics.cpp`, `test_optimizer.cpp`, `test_simulator.cpp`, `test_chemistry.cpp`, `test_compiler_facade.cpp`, `test_simulation_session.cpp`, `test_workbench_integration.cpp`, `test_experiment_bundle.cpp`, `test_molecular_integral_engine.cpp`, `test_circuit_dag.cpp`, `test_coupling_map.cpp`, and `test_vqe_optimizer.cpp`. Shared fixtures and assertions are in `test_support.h` / `test_support.cpp`.

### 3. Molecular Chemistry Benchmarks
Run the fixed molecular benchmarks and catalog CLI:
```bash
# Fixed H2 VQE benchmark (2 qubits, parity mapping, STO-3G at 0.735 Å)
./qcompile_h2

# Fixed HeH+ VQE benchmark (4 qubits, Jordan-Wigner, STO-3G at 1.463283 bohr)
./qcompile_heh

# Native STO-3G H2 integral-engine checks (distances in Angstrom)
./qcompile_benchmarks --engine-h2 0.735
./qcompile_benchmarks --engine-h2 1.400

# Fixed LiH (1 frozen core orbital, 2 active orbitals; exact active-space benchmark)
./qcompile_benchmarks --engine-lih

# Independent PySCF LiH active-space reference
benchmarks/heh_plus/.venv-generator/bin/python -B \
  benchmarks/lih_active_space/validate_pyscf_lih.py

# Benchmark catalog CLI
./qcompile_benchmarks --list
./qcompile_benchmarks h2
./qcompile_benchmarks heh_plus
./qcompile_benchmarks --json-list
./qcompile_benchmarks --json h2
```

### 4. Interactive Qt6 H₂ Workbench
Launch the native graphical workbench:
```bash
./qcompile_lab
```
Run the automated offscreen smoke verification:
```bash
QT_QPA_PLATFORM=offscreen ./qcompile_lab --smoke-test
```

### 5. Loopback Scientific Dashboard
Start the local HTTP dashboard:
```bash
node dashboard/server.js
```
Open `http://127.0.0.1:4173` in a web browser. The server executes allowlisted benchmarks via `qcompile_benchmarks --json <id>` over loopback.

---

## Technical Limitations

- **Statevector Scaling**: The simulator allocates $2^n$ complex amplitudes (`std::complex<double>`) and supports at most 24 qubits; one 24-qubit statevector uses about 256 MiB. Interactive `SimulationSession`s additionally cap the combined statevector payload of the active simulator and retained step snapshots at 512 MiB, excluding allocator and other process overhead. At 24 qubits the initial session uses that full cap, so its first execution step is rejected before applying the instruction.
- **Measurement Model**: Measurements calculate probabilities, sample basis outcomes using a seeded PRNG (`std::mt19937`), project the statevector onto the observed subspace, and renormalize, preserving conditional states and entanglement on unmeasured qubits.
- **Gate Set**: Supported unitary operations are restricted to `H`, `X`, `Z`, `CNOT`, and `RZ`. General arbitrary single-qubit rotations or multi-qubit controlled gates require decomposition before input.
- **Optimizer Scope**: The DAG optimizer cancels recognized matching H/X/Z and ordered CNOT pairs and merges finite RZ rotations. It can bring eligible operations together across explicitly commuting gates, including disjoint supports and specified Z/RZ–CNOT cases; measurements are barriers. It does not perform general gate synthesis or hardware scheduling.
- **Chemistry Benchmarks**: The H₂, HeH⁺, and LiH commands validate fixed STO-3G models and specified symmetry/active-space sectors. They are not an arbitrary quantum chemistry platform; geometry-generated support is restricted to the documented H₂ path.

---

## References and Documentation

For complete technical specifications, language grammar, error catalogs, eigensolver derivations, and chemistry provenance, refer to the [Technical Reference Manual (`DOCUMENTATION.md`)](DOCUMENTATION.md).

The current scientific benchmark evidence and reproduction notes are consolidated in [the MS10 validation report](docs/MS10_FINAL_REPORT.md).
