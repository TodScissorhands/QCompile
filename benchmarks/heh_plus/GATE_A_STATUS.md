# HeH+ Reference Data — Gate A Status

**Verdict: BLOCKED.** The fixed Hamiltonian and an independent determinant-space matrix check were generated successfully. Psi4 1.9.1 imports, but its SCF/FCI execution terminates with signal 11 (exit 139) during basis loading. No independent Psi4 FCI energy was obtained; Gate A is not passed.

## Fixed model

- HeH+, nuclei He (0,0,0), H (0,0,1.463283) bohr; 1.463283 bohr is the approved high-accuracy physical BO equilibrium geometry, not an STO-3G optimized geometry.
- STO-3G, charge +1, 2 electrons, singlet S=0 (multiplicity 1), RHF orbitals, all 2 spatial orbitals retained, no frozen core.
- Four modes/qubits; Jordan–Wigner, no symmetry tapering. Mode/occupation order is q0=MO0 alpha, q1=MO1 alpha, q2=MO0 beta, q3=MO1 beta. State-vector integer index uses q0 as least-significant bit. Qiskit Pauli strings print q3…q0, so the rightmost label character acts on q0.
- Electronic Hamiltonian excludes nuclear repulsion. Nuclear repulsion is 2/R = 1.3667896093920315 Ha and is added once to form total Born–Oppenheimer energy.

## Generated data and internal validation

`reference_data.json` contains full binary64-decimal serialization of the one-electron MO integrals, full chemist-order two-electron MO tensor, 72 Qiskit Nature fermionic terms, all 27 Jordan–Wigner Pauli terms and coefficients, ordering conventions, settings, and matrix/energy results. The generator ran using Qiskit Nature 0.8.0 and PySCF 2.8.0 at RHF convergence tolerance 1e-12.

The determinant matrix was assembled from the frozen h[p,q] and chemist-ordered (pq|rs) integrals by applying spin-orbital creation/annihilation operators to all 16 occupation-bit basis states. Its entries were compared to both the Qiskit Nature fermionic operator and `SparsePauliOp.to_matrix()`. Maximum absolute residuals:

- Integral-built vs Qiskit Nature fermionic matrix: **1.1102230246251565e-16 Ha**.
- Integral-built vs mapped Pauli matrix: **4.884981308350689e-15 Ha**.
- Mapped Pauli vs fermionic matrix: **4.884981308350689e-15 Ha**.
- Hermiticity residual: **1.1102230246251565e-16 Ha**.
- Direct PySCF/Qiskit Nature one-electron integral and nuclear-repulsion residuals: **0** at binary64 precision.

The fixed electron sector is Nα=Nβ=1 (4 determinants). Independently assembled S² has eigenvalues 0,0,0,2 in that sector; only the three S=0 vectors are retained for the reference minimum. Direct integral-matrix singlet minimum: **−4.218264154490395 Ha electronic**, nuclear repulsion **+1.3667896093920315 Ha**, total **−2.8514745450983634 Ha**. These are generator-side results, not yet cross-validated against Psi4.

The internal comparison threshold was fixed at 2e-12 Ha before comparing matrices. Observed residuals are below it. No tolerance was selected from or relaxed to fit an external result.

## Pauli coefficient table

All coefficients are in Hartree; labels are Qiskit strings q3…q0. The complete definitive table is in `reference_data.json`.

| Pauli string | Coefficient |
|---|---:|
| IIII | -2.9120172734747536 |
| IIIZ | +0.7583883292075204 |
| IIYY | +0.05256451032097764 |
| IIXX | +0.05256451032097764 |
| IIZI | +0.19189724869368352 |
| IZII | +0.7583883292075202 |
| YYII | +0.05256451032097764 |
| XXII | +0.05256451032097764 |
| ZIII | +0.19189724869368363 |
| IIZZ | +0.12871109487820281 |
| IZIZ | +0.23577440756967472 |
| YYIZ | +0.04324233541007819 |
| XXIZ | +0.04324233541007819 |
| ZIIZ | +0.16505925176240369 |
| IZYY | +0.04324233541007819 |
| IZXX | +0.04324233541007819 |
| YYYY | +0.03634815688420088 |
| XXYY | +0.03634815688420088 |
| YYXX | +0.03634815688420088 |
| XXXX | +0.03634815688420088 |
| ZIYY | -0.009322174907132445 |
| ZIXX | -0.009322174907132445 |
| IZZI | +0.16505925176240369 |
| YYZI | -0.009322174907132445 |
| XXZI | -0.009322174907132445 |
| ZIZI | +0.18813101682146008 |
| ZZII | +0.12871109487820281 |

## Environment and reproduction

- New reference tools/data and isolated environments are under this previously absent directory; no existing project source/configuration/test/report was changed. **Exception:** Psi4 emitted a repository-root `timer.dat` during its failed runs. It was absent from the initial `git status` summary and appeared afterward; it is preserved untouched rather than removed or moved. No C++ build dependency was added.

- Generator: Python **3.11.9**, NumPy **2.4.6**, PySCF **2.8.0**, Qiskit **2.5.2**, Qiskit Nature **0.8.0**, uv **0.13.1**. Exact generator transitive pins: `requirements-generator.lock`.
- Independent package: Psi4 **1.9.1**, Python **3.11.17**, conda-forge build `py311h2ae924c_7`, libxc-c **6.2.2** build `cpu_h20a523f_6`. Exact Linux-64 package URLs: `psi4-linux-64.explicit.txt`.
- Executed generator setup: `uv python install 3.11.9`; `uv venv --python 3.11.9 benchmarks/heh_plus/.venv-generator`; `uv pip install --python benchmarks/heh_plus/.venv-generator/bin/python 'qiskit-nature==0.8.0' 'pyscf==2.8.0'`. `uv pip freeze` recorded all installed packages in the lock file. Psi4 was installed into the local prefix from conda-forge with `psi4=1.9.1 python=3.11`, then libxc-c was pinned to 6.2.2 after the initial import failed. Reproduction environment creation was local; no C++ build dependencies were installed.
- Re-run generated Hamiltonian/matrix data: `benchmarks/heh_plus/.venv-generator/bin/python benchmarks/heh_plus/generate_reference.py` — **exit 0**.
- Attempt independent FCI: `benchmarks/heh_plus/.conda-psi4/bin/python benchmarks/heh_plus/run_psi4_fci.py` — **exit 139**. Output ends while loading the STO-3G basis. A first Psi4 import emitted `Functional 'XC_GGA_XC_TH_FL' not found`; pinning libxc-c 6.2.2 removed that import-time error, but RHF/FCI still segfaulted during basis loading. `psi4_fci_reference.json` was therefore not produced. `psi4_fci_output.txt` preserves the observed partial log.

Useful primary documentation: [Qiskit Nature PySCFDriver](https://qiskit-community.github.io/qiskit-nature/stubs/qiskit_nature.second_q.drivers.PySCFDriver.html), [Qiskit Nature mapper tutorial](https://qiskit-community.github.io/qiskit-nature/tutorials/06_qubit_mappers.html), [Psi4 DETCI/FCI manual](https://psicode.org/psi4manual/1.9.x/detci.html), and Pachucki, *Accurate Born-Oppenheimer potential for HeH+* ([arXiv:1203.3927](https://arxiv.org/abs/1203.3927)). The published high-accuracy BO equilibrium geometry is not itself an STO-3G FCI reference energy.

## Remaining blocker

Obtain a Psi4 runtime/build that executes on this host, or identify and pin a compatible conda-forge dependency set without changing the approved benchmark. Then rerun the exact Psi4 FCI calculation, record its spin diagnostic and separate electronic/nuclear/total energies, and compare it against the direct integral-sector minimum. Until that independent comparison succeeds within a justified tolerance, Gate A remains **BLOCKED**. No C++ HeH+ implementation is authorized or included.
