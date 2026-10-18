# LiH frozen-core active-space reference

This validation compares the fixed QCompile LiH model against an independent PySCF CASCI calculation. It is a calculation-specific cross-check, not an experimental energy or a general molecular benchmark.

## Model and energy convention

- Neutral singlet LiH at `R = 1.595 Å`; coordinates are converted with `1 bohr = 0.52917721092 Å`.
- STO-3G basis, RHF molecular orbitals, four electrons.
- Freeze the lowest spatial orbital (Li 1s); retain two active spatial orbitals and two active electrons (`Nα=Nβ=1`). Remaining orbitals are inactive virtuals.
- The active-space qubit Hamiltonian contains only the active-space electronic operator. PySCF's CASCI core scalar includes the nuclear-repulsion energy and frozen-core electronic contribution. Total CASCI energy is `Ecore + Eactive`; the scalar must not be added to the qubit Hamiltonian and then counted again.
- The cross-check uses PySCF's RHF/CASCI implementation and is independent of QCompile's integral and determinant-matrix implementation, but shares the same basis, geometry, and active-space definition.

## Reproduction

From the repository root, using the pinned reference environment described in `benchmarks/h2_geometry/README.md`:

```bash
benchmarks/heh_plus/.venv-generator/bin/python -B \
  benchmarks/lih_active_space/validate_pyscf_lih.py

make all
./qcompile_benchmarks --engine-lih
```

The reference script uses RHF energy convergence `1e-12 Ha`, at most 100 iterations, and a MINAO initial guess. It checks the reference values in `tests/test_molecular_integral_engine.cpp` within `1e-7 Ha` and verifies the `Etotal = Ecore + Eactive` identity within `1e-12 Ha`. It prints the measured values, exact model settings, and residuals as JSON. The C++ executable reports its values separately; compare the reported model and energy conventions before comparing totals.

The C++ result and independent PySCF reference are expected to agree within `1e-7 Ha`; small binary64 residuals reflect independently implemented integral transformations and solvers. The two-electron active-space energy is electronic-only. `Ecore` includes the nuclear term and frozen-core contribution and is added exactly once for a total Born–Oppenheimer energy.

PySCF's [CASCI documentation](https://pyscf.org/user/mcscf.html#casci) describes the complete-active-space CI method used for the independent reference. The numerical values are model-specific PySCF results, not experimental molecular energies.
