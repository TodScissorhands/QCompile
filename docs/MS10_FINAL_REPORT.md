# MS10 — Scientific Documentation and Benchmark Validation

## Outcome

**MS10 complete.** All scope-derived exit criteria below passed; no material in-scope blocker remains. No standalone MS10 specification exists in the repository: before this report, a tracked-source search found no MS10/Milestone 10 document, and `DOCUMENTATION.md`'s hardening phase register ends at Phase 9. This work uses the explicit MS10 requirements supplied in the task as its controlling specification; it does not redefine that scope.

| Exit criterion | Result | Evidence |
|---|---|---|
| Identify and complete the explicit MS10 deliverables | Pass | Scientific claims, benchmark procedures, report, documentation, and validation evidence are consolidated here and in the linked benchmark notes. |
| Support scientific and benchmark claims with reproducible evidence | Pass | Fresh QCompile, PySCF/Qiskit Nature, Psi4, and PySCF CASCI runs; model and tolerances recorded below. |
| Validate required results against appropriate references | Pass | H₂ and HeH⁺ cross-package FCI checks; LiH independent CASCI comparison. Comparisons hold basis, geometry, and active-space model fixed. |
| Correct material discrepancies | Pass | Stale HeH⁺ Gate A “blocked” status and missing Psi4 spin diagnostic corrected; stale README test count and incomplete module inventory corrected. |
| Align documentation with implementation, assumptions, and limitations | Pass | Fixed-model, electronic-vs-total energy, LiH core-scalar accounting, reference provenance, and supported-scope limitations documented. |
| Consolidate a final report without duplicating obsolete MS10 reports | Pass | This is the consolidated MS10 report; detailed model-specific reproduction notes remain with their benchmarks. |
| Run relevant tests, builds, benchmarks, and repository checks | Pass | The full suite passed 1,519 assertions after the final label correction. The earlier ASan/UBSan suite passed 1,519 before that output-label-only correction; chronology is stated below. Required benchmark paths and `git diff --check` passed. |
| Keep the changeset limited to intended MS10 work | Pass | The initial worktree was clean; the reviewed delta contains only MS10 work. The report's Git snapshot is explicitly pre-publication. |
| Leave no material in-scope blocker | Pass | No unresolved blocker. Fixed-model and general-chemistry limitations are explicitly documented. |

## Independent findings and changes

1. **HeH⁺ reference status contradicted the retained evidence.** `benchmarks/heh_plus/GATE_A_STATUS.md` still described Gate A as blocked after an earlier Psi4 process crashed, while a later Libint2 retry output and JSON already recorded a successful DETCI FCI run. Re-executing the retry in a temporary output directory with Psi4 1.9.1 reproduced the result. The report now marks the gate accepted, preserves the earlier crash as provenance, records the matching model and energy convention, and provides safe temporary-directory reproduction commands.
2. **The retained Psi4 spin diagnostic was incomplete.** The Python interface did not expose the value in the way the script expected, but the Psi4 log reported `<S²> = 0.0`. `run_psi4_fci_libint290_attempt1_retry.py` now reads the last `Computed <S^2> vector 0 =` entry as a fallback; a fresh run produced JSON `spin_diagnostic_s_squared: 0.0`. The retained retry JSON was corrected accordingly.
3. **The HeH⁺ geometry source could be confused with an STO-3G energy reference.** Pachucki's verified high-accuracy Born–Oppenheimer result gives `R = 1.463283 a₀` and `E = -2.978708310771 Ha`. That is not the STO-3G model used by QCompile. The Gate A report now states the geometry provenance and explicitly excludes that high-accuracy energy from STO-3G comparisons. PySCF/Qiskit Nature and Psi4 comparisons below use the same STO-3G model.
4. **README test information was stale.** Its current assertion count said 391 and its test-module list was incomplete. The README no longer pins a transient assertion total and lists the current modules and benchmark reproduction paths. `DOCUMENTATION.md`'s Phase 9 table retains its explicitly labelled historical “Reported Verification Summary”; that separate CLI-hardening entry is not the current suite result and was not rewritten as MS10 scope.
5. **LiH needed an independently runnable active-space cross-check.** Added `benchmarks/lih_active_space/validate_pyscf_lih.py` and its README. It computes a fixed PySCF RHF/CASCI reference, checks the active-space model and `E_total = E_core + E_active`, and checks component/total references to `1e-7 Ha`. This is a calculation-specific numerical reference, not an experimental energy.
6. **LiH frozen-core accounting is consistent.** The C++ path stores the frozen-core scalar separately (`src/chemistry/molecular_integral_engine.cpp` around the frozen-core construction and result assembly); the active Pauli operator contains only active-space electronic terms. Tests assert the scalar-plus-active energy identity. An independent PySCF CASCI run and a separate two-qubit LiH VQE execution confirmed the documented accounting.
7. **LiH CLI terminology overstated the active-space energy as correlation energy.** `fciExactElectronicEnergy` is the eigenvalue of the complete frozen-core-adjusted active-space electronic Hamiltonian, not the correlation-energy difference from RHF. Changed the human-readable output to `E_active` and corrected the corresponding regression-test descriptions; the calculated values and JSON field are unchanged.

## Scientific and benchmark validation

All energies are in Hartree. Nuclear repulsion is excluded from the electronic/qubit Hamiltonians and included once in the reported total unless the LiH value is identified as `E_core`, which already includes the nuclear term.

| Model and method | Independent/reference energy | QCompile result | Agreement / interpretation |
|---|---:|---:|---|
| H₂, fixed rounded STO-3G Hamiltonian at 0.735 Å; parity-reduced two-qubit singlet | Qiskit Nature electronic reference `-1.857275030202000` | Fixed-model exact `-1.857275029288229`; VQE electronic `-1.857275029288228`; total `-1.137306034839228` | Rounded coefficients shift the electronic minimum by `9.14e-10 Ha`, below the documented `1e-9 Ha`; this is not a disagreement in the generated full-precision model. |
| H₂, freshly generated STO-3G at 0.735 Å; independent Psi4 DETCI FCI | PySCF/Qiskit Nature FCI `-1.857275030202380`; Psi4 `-1.857275030202397` | Generated FCI `-1.857275030202399`; optimized `-1.857275030202398`; total `-1.137306035753399` | Psi4-vs-generator electronic and total differences `1.69e-14 Ha` (script limits `1e-8 Ha`); nuclear difference `0` (limit `1e-12 Ha`). Generator matrix residual `8.88e-16 Ha` (limit `2e-12 Ha`). |
| H₂, freshly generated STO-3G at 1.400 Å; independent Psi4 DETCI FCI | PySCF/direct singlet FCI `-1.393451971373960`; Psi4 `-1.393451971373954` | Generated FCI `-1.393451971373955`; optimized `-1.393451971373954`; total `-1.015468249288244` | Psi4-vs-generator electronic and total differences `5.77e-15 Ha`; nuclear difference `0`. Generator matrix residual `5.55e-16 Ha` (limit `2e-12 Ha`). |
| HeH⁺, STO-3G at 1.463283 bohr, 2 electrons, singlet, 2 spatial orbitals, Jordan–Wigner | PySCF/Qiskit Nature `-4.218264154490395`; Psi4 DETCI `-4.218264154490385`; Psi4 `<S²> = 0.0` | Exact singlet `-4.218264154490395`; variational `-4.218264154490334`; total `-2.851474545098302` | Independent electronic difference `1.0e-14 Ha`; total-reference difference `9.8e-15 Ha`. Direct determinant/Pauli matrix maximum residual `4.88e-15 Ha`. VQE variational excess over exact is `6.1e-14 Ha`; 6,803 objective evaluations. |
| LiH, STO-3G at 1.595 Å, neutral 4-electron singlet; freeze 1 core orbital, retain 2 active orbitals and 2 active electrons | PySCF 2.8.0 CASCI: `E_core = -6.802973549986354`, `E_active = -1.059311983913182`, `E_total = -7.862285533899536` | C++ active-space exact: `E_core = -6.80297354995545`, `E_active = -1.05931198390915`, `E_total = -7.8622855338646`; two-qubit VQE total `-7.8622855338646` | C++/PySCF total difference `3.49e-11 Ha`; PySCF identity residual `0`; independent script limits component and total differences to `1e-7 Ha`. VQE used 323 objective evaluations and matched the C++ active-space exact value at printed precision. |

### Model interpretation and limitations

- The fixed H₂ CLI uses a rounded, fixed Hamiltonian; the generated H₂ path recalculates STO-3G integrals for the requested geometry. Their small 0.735 Å difference is expected from coefficient rounding. The fresh generated cases were checked by direct determinant-space matrix assembly, PySCF FCI, and an independent Psi4 DETCI run.
- The HeH⁺ bond length comes from the high-accuracy equilibrium geometry in [Pachucki (2012)](https://arxiv.org/abs/1203.3927). That paper's high-accuracy energy is not the same basis/model and is not used as a STO-3G reference. The model-matched comparison is PySCF/Qiskit Nature versus Psi4.
- The LiH validation is for one fixed RHF/STO-3G geometry and active space. `E_core` includes nuclear repulsion and the frozen-core electronic contribution; it must not be placed in the qubit Hamiltonian and then added again. The VQE expectation is active electronic energy, and total energy is `E_core + E_active` exactly once.
- These benchmarks do not demonstrate arbitrary-molecule support, arbitrary basis/active-space generation, experimental accuracy, or quantum advantage. Geometry-generated support remains restricted to the documented H₂ route; LiH and HeH⁺ use fixed models.

## Reproduction procedures and execution evidence

Commands below were executed during this MS10 work. Full command output was observed; the listed exit status is the process result.

### Project build and tests

```bash
make -j2 all && make test
```

Exit `0`; output: `Passed 1519 assertions.` The current suite includes optimizer, chemistry, molecular integral, VQE, compiler/simulator, and workbench-related tests.

After the final terminology correction, `make -j2 all && make test && ./qcompile_benchmarks --engine-lih` exited `0`; the full suite again passed 1,519 assertions. The CLI now reports `Active-space FCI electronic energy (E_active)` and `Active-space exact total energy (E_core + E_active)`.

### Fixed/native chemistry benchmark commands

```bash
./qcompile_h2
./qcompile_heh
./qcompile_benchmarks h2 && ./qcompile_benchmarks heh_plus
./qcompile_benchmarks --engine-lih
./qcompile_benchmarks --engine-h2 0.735 && ./qcompile_benchmarks --engine-h2 1.400
```

All exited `0`. Outputs are reported in the table above. The fixed H₂ and HeH⁺ catalog runs each reported `Converged: yes`; the generated H₂ optimizers reported 62 evaluations at both distances.

### Fresh H₂ integral generation, C++ execution, and Psi4 comparison

The generated files were written only under `/tmp/qcompile-ms10-h2-validation-eicq0d89/`.

```bash
/home/tod/QCompile/benchmarks/heh_plus/.venv-generator/bin/python -B /home/tod/QCompile/benchmarks/h2_geometry/generate_h2.py --distance-angstrom 0.735 --output /tmp/qcompile-ms10-h2-validation-eicq0d89/h2-0735.qham
/home/tod/QCompile/qcompile_benchmarks --generated-h2 /tmp/qcompile-ms10-h2-validation-eicq0d89/h2-0735.qham
/home/tod/QCompile/benchmarks/heh_plus/.conda-psi4-libint290-attempt1/bin/python -B /home/tod/QCompile/benchmarks/h2_geometry/validate_psi4_h2.py --distance-angstrom 0.735 --generator-provenance /tmp/qcompile-ms10-h2-validation-eicq0d89/h2-0735.qham.provenance.json --output-json /tmp/qcompile-ms10-h2-validation-eicq0d89/psi4-0735.json --scratch /tmp/qcompile-ms10-h2-validation-eicq0d89/psi4-scratch-0735

/home/tod/QCompile/benchmarks/heh_plus/.venv-generator/bin/python -B /home/tod/QCompile/benchmarks/h2_geometry/generate_h2.py --distance-angstrom 1.400 --output /tmp/qcompile-ms10-h2-validation-eicq0d89/h2-1400.qham
/home/tod/QCompile/qcompile_benchmarks --generated-h2 /tmp/qcompile-ms10-h2-validation-eicq0d89/h2-1400.qham
/home/tod/QCompile/benchmarks/heh_plus/.conda-psi4-libint290-attempt1/bin/python -B /home/tod/QCompile/benchmarks/h2_geometry/validate_psi4_h2.py --distance-angstrom 1.400 --generator-provenance /tmp/qcompile-ms10-h2-validation-eicq0d89/h2-1400.qham.provenance.json --output-json /tmp/qcompile-ms10-h2-validation-eicq0d89/psi4-1400.json --scratch /tmp/qcompile-ms10-h2-validation-eicq0d89/psi4-scratch-1400
```

All six commands exited `0`. The generator checked integral/matrix consistency, singlet-sector energies, and parity/Pauli ordering; Psi4 reported `S² = 0.0` at both geometries and passed its documented electronic/total (`1e-8 Ha`) and nuclear (`1e-12 Ha`) tolerances. A first Psi4 attempt during exploratory validation used a relative script path while the process working directory was `/tmp` and exited `2` (file not found); the absolute-path invocations above were the successful corrected runs. No source or repository output was produced by the failed invocation.

### HeH⁺ independent reference reproduction

The two generators were copied with Python `shutil.copy2` into `/tmp/qcompile-ms10-heh-validation-ebeuzmsn/` so their output files and Psi4 scratch data could not overwrite repository evidence. The executed commands were:

```bash
/home/tod/QCompile/benchmarks/heh_plus/.venv-generator/bin/python -B /tmp/qcompile-ms10-heh-validation-ebeuzmsn/generate_reference.py
/home/tod/QCompile/benchmarks/heh_plus/.conda-psi4-libint290-attempt1/bin/python -B /tmp/qcompile-ms10-heh-validation-ebeuzmsn/run_psi4_fci_libint290_attempt1_retry.py
```

Both exited `0`. PySCF/Qiskit Nature FCI returned `-4.218264154490395 Ha`; Psi4 returned `-4.218264154490385 Ha`, nuclear repulsion `1.3667896093920315 Ha`, total `-2.8514745450983536 Ha`, and `S² = 0.0`. The independent determinant/PySCF/Qiskit matrix residuals were at most `4.88e-15 Ha`.

### LiH independent CASCI and VQE reproduction

```bash
benchmarks/heh_plus/.venv-generator/bin/python -B benchmarks/lih_active_space/validate_pyscf_lih.py
```

Exit `0`. PySCF 2.8.0 returned `E_core = -6.802973549986354 Ha`, `E_active = -1.0593119839131822 Ha`, `E_total = -7.862285533899536 Ha`; energy-identity residual `0`. Differences from checked-in C++ regression references were at most `1.82e-13 Ha`, below `1e-7 Ha`.

A throwaway C++ harness in `/tmp/qcompile-ms10-lih-vqe-tl8pgzpo/` compiled the same optimizer/ansatz/molecular-engine sources and copied the options from `tests/test_vqe_optimizer.cpp`. Compile and run both exited `0`; output: `core=-6.80297354995545 active_vqe=-1.05931198390915 total=-7.8622855338646 exact_total=-7.8622855338646 evaluations=323`. The harness was not added to the repository.

### Sanitizers and whitespace check

The prior full-suite sanitizer validation used ASan/UBSan, leak detection, and halt-on-error. Its exact command was:
```bash
g++ -std=c++17 -Wall -Wextra -pedantic -g -O1 -fsanitize=address,undefined -fno-omit-frame-pointer -Iinclude tests/qcompile_tests.cpp tests/test_support.cpp tests/test_lexer.cpp tests/test_parser.cpp tests/test_semantics.cpp tests/test_optimizer.cpp tests/test_simulator.cpp tests/test_chemistry.cpp tests/test_compiler_facade.cpp tests/test_simulation_session.cpp tests/test_workbench_integration.cpp tests/test_experiment_bundle.cpp tests/test_molecular_integral_engine.cpp tests/test_circuit_dag.cpp tests/test_coupling_map.cpp tests/test_vqe_optimizer.cpp src/compiler/lexer.cpp src/compiler/parser.cpp src/compiler/semantic_checker.cpp src/compiler/code_generator.cpp src/compiler/circuit_dag.cpp src/compiler/coupling_map.cpp src/compiler/router.cpp src/compiler/compiler_facade.cpp src/simulator/simulator.cpp src/simulator/simulation_session.cpp src/chemistry/hamiltonian.cpp src/chemistry/h2_energy.cpp src/chemistry/h2_geometry.cpp src/chemistry/heh_energy.cpp src/chemistry/benchmark_catalog.cpp src/experiment/experiment_bundle.cpp src/chemistry/molecular_integral_engine.cpp src/chemistry/nelder_mead.cpp src/chemistry/vqe_ansatz.cpp -o /tmp/qcompile-ms10-qcompile-tests-asan && ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 /tmp/qcompile-ms10-qcompile-tests-asan
```
Exit `0`; output `Passed 1519 assertions.` No ASan, UBSan, or leak diagnostics were observed. This run preceded the final correction to the LiH output label and test text. Those final changes did not alter calculation logic; the regular full suite was rerun after them.

```bash
git diff --check
```

Exit `0`; no output.

## Changed files

- `README.md` — removed stale current assertion count, completed the test-module list, documented fixed/generated chemistry commands, and linked this report.
- `DOCUMENTATION.md` — documented LiH active-space energy accounting and the independently reproduced HeH⁺ result with exact evidence paths.
- `benchmarks/heh_plus/GATE_A_STATUS.md` — replaced stale blocked status with current accepted model, cross-check, safe reproduction instructions, and model-comparison caveat.
- `benchmarks/heh_plus/run_psi4_fci_libint290_attempt1_retry.py` — added a Psi4-log fallback for the `<S²>` diagnostic.
- `benchmarks/heh_plus/psi4_fci_reference_libint290_attempt1_retry.json` — recorded the verified `S² = 0.0` value.
- `benchmarks/lih_active_space/README.md` — documented the fixed model, energy convention, and reproducible independent reference.
- `benchmarks/lih_active_space/validate_pyscf_lih.py` — added a reproducible PySCF CASCI reference and energy-identity/tolerance checks.
- `src/apps/benchmark_main.cpp` — labels the LiH active-space eigenvalue as electronic energy, not correlation energy.
- `tests/test_molecular_integral_engine.cpp` — aligns LiH reference and exactly-once identity assertions with `E_active` terminology.
- `docs/MS10_FINAL_REPORT.md` — this consolidated milestone report.

No unrelated implementation changes or tracked reports were included. Temporary data and throwaway harnesses remained outside the repository.

## Pre-publication review snapshot

This snapshot was captured after final-diff review and before staging or publication.

- Branch: `main`
- HEAD: `8d9d32189e19bd21d819a39fd27b5facc76441ad`
- Live `origin/main`: `8d9d32189e19bd21d819a39fd27b5facc76441ad` (`git ls-remote origin refs/heads/main`); synchronized at review time.
- `git status --short --untracked-files=all`:
  ```text
   M DOCUMENTATION.md
   M README.md
   M benchmarks/heh_plus/GATE_A_STATUS.md
   M benchmarks/heh_plus/psi4_fci_reference_libint290_attempt1_retry.json
   M benchmarks/heh_plus/run_psi4_fci_libint290_attempt1_retry.py
   M src/apps/benchmark_main.cpp
   M tests/test_molecular_integral_engine.cpp
  ?? benchmarks/lih_active_space/README.md
  ?? benchmarks/lih_active_space/validate_pyscf_lih.py
  ?? docs/MS10_FINAL_REPORT.md
  ```
- `git diff --stat`: 7 tracked files changed, 98 insertions, 91 deletions. This excludes the three untracked files listed above.
- `git diff --check`: exit `0`, no output.
- `git diff --cached --stat`: empty before publication.

The commit and push authorized after this snapshot are recorded in the final Git history.

## References

- [Qiskit Nature 0.8 electronic-structure tutorial](https://qiskit-community.github.io/qiskit-nature/tutorials/01_electronic_structure.html) — H₂/STO-3G setup and the distinction between the electronic Hamiltonian and nuclear-repulsion offset.
- [PySCF CASCI documentation](https://pyscf.org/user/mcscf.html#casci) — complete-active-space CI method used for the LiH reference.
- [K. Pachucki (2012), *Accurate Born-Oppenheimer potential for HeH⁺*, arXiv:1203.3927](https://arxiv.org/abs/1203.3927) — high-accuracy HeH⁺ equilibrium geometry and a distinct, non-STO-3G energy calculation.
