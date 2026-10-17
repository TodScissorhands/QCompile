#include "test_support.h"
#include "chemistry/molecular_integral_engine.h"
#include "chemistry/h2_energy.h"
#include "chemistry/heh_energy.h"
#include "simulator/simulation_session.h"
#include "chemistry/symmetric_eigensolver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace test {

void testMolecularIntegralEngine() {
  // 1. Molecular Geometry Validation & Nuclear Repulsion
  {
    // Empty geometry
    MolecularGeometry emptyGeom;
    bool emptyRejected = false;
    try {
      emptyGeom.validate();
    } catch (const std::invalid_argument &) {
      emptyRejected = true;
    }
    require(emptyRejected, "empty molecular geometry rejected");

    // Unsupported element
    MolecularGeometry badElementGeom;
    badElementGeom.atoms = {{"U", 92, {0.0, 0.0, 0.0}}};
    bool unsupportedRejected = false;
    try {
      badElementGeom.validate();
    } catch (const std::invalid_argument &) {
      unsupportedRejected = true;
    }
    require(unsupportedRejected, "unsupported element rejected");

    // Coincident nuclei
    MolecularGeometry coincidentGeom;
    coincidentGeom.atoms = {
        {"H", 1, {0.0, 0.0, 0.0}},
        {"H", 1, {0.0, 0.0, 1e-6}}
    };
    bool coincidentRejected = false;
    try {
      coincidentGeom.validate();
    } catch (const std::invalid_argument &) {
      coincidentRejected = true;
    }
    require(coincidentRejected, "coincident nuclei rejected");

    // Valid H2 geometry at 0.735 Angstrom
    MolecularGeometry h2Geom;
    h2Geom.atoms = {
        {"H", 1, {0.0, 0.0, 0.0}},
        {"H", 1, {0.0, 0.0, 0.735}}
    };
    h2Geom.unit = GeometryUnit::Angstrom;
    h2Geom.charge = 0;
    h2Geom.multiplicity = 1;

    require(h2Geom.atoms.size() == 2, "H2 geometry has 2 atoms");
    const auto bohrAtoms = h2Geom.inBohr();
    require(std::abs(bohrAtoms[1].coordinates[2] - (0.735 / 0.52917721092)) < 1e-10,
            "Angstrom to Bohr conversion matches CODATA value");

    // Nuclear repulsion energy check against analytical 1.0 / R_bohr
    const double expectedH2Nuc = 1.0 / (0.735 / 0.52917721092);
    require(std::abs(h2Geom.nuclearRepulsionEnergy() - expectedH2Nuc) < 1e-12,
            "H2 nuclear repulsion matches analytical Coulomb law");
    require(std::abs(h2Geom.nuclearRepulsionEnergy() - 0.719968994449) < 1e-8,
            "H2 nuclear repulsion matches documented reference (0.719969 Ha)");

    // Valid HeH+ geometry at 1.463283 Bohr
    MolecularGeometry hehGeom;
    hehGeom.atoms = {
        {"He", 2, {0.0, 0.0, 0.0}},
        {"H", 1, {0.0, 0.0, 1.463283}}
    };
    hehGeom.unit = GeometryUnit::Bohr;
    hehGeom.charge = 1;
    hehGeom.multiplicity = 1;

    const double expectedHeHNuc = (2.0 * 1.0) / 1.463283;
    require(std::abs(hehGeom.nuclearRepulsionEnergy() - expectedHeHNuc) < 1e-12,
            "HeH+ nuclear repulsion matches analytical Coulomb law");
    require(std::abs(hehGeom.nuclearRepulsionEnergy() - 1.366789609392) < 1e-8,
            "HeH+ nuclear repulsion matches documented reference (1.366790 Ha)");
  }

  // 2. Integral Matrix Symmetries & Invariants
  {
    MolecularGeometry h2Geom;
    h2Geom.atoms = {{"H", 1, {0.0, 0.0, 0.0}}, {"H", 1, {0.0, 0.0, 0.735}}};
    h2Geom.unit = GeometryUnit::Angstrom;

    const auto integrals = computeMolecularIntegrals(h2Geom);
    require(integrals.nao == 2, "2 STO-3G basis functions for H2");

    // Overlap matrix S symmetry & normalization
    require(std::abs(integrals.getOverlap(0, 0) - 1.0) < 1e-7, "S_00 is normalized to 1.0");
    require(std::abs(integrals.getOverlap(1, 1) - 1.0) < 1e-7, "S_11 is normalized to 1.0");
    require(std::abs(integrals.getOverlap(0, 1) - integrals.getOverlap(1, 0)) < 1e-14,
            "Overlap matrix is symmetric (S_01 == S_10)");

    // Kinetic matrix T symmetry & positivity
    require(integrals.getKinetic(0, 0) > 0.0, "T_00 is positive");
    require(integrals.getKinetic(1, 1) > 0.0, "T_11 is positive");
    require(std::abs(integrals.getKinetic(0, 1) - integrals.getKinetic(1, 0)) < 1e-14,
            "Kinetic matrix is symmetric (T_01 == T_10)");

    // Nuclear attraction matrix V symmetry & negativity
    require(integrals.getNuclear(0, 0) < 0.0, "V_00 is negative (attractive)");
    require(integrals.getNuclear(1, 1) < 0.0, "V_11 is negative (attractive)");
    require(std::abs(integrals.getNuclear(0, 1) - integrals.getNuclear(1, 0)) < 1e-14,
            "Nuclear matrix is symmetric (V_01 == V_10)");

    // 8-fold ERI permutation symmetry: (ij|kl) == (ji|kl) == (ij|lk) == (kl|ij) etc.
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
          for (int l = 0; l < 2; ++l) {
            const double val = integrals.getEri(i, j, k, l);
            require(std::abs(val - integrals.getEri(j, i, k, l)) < 1e-14, "ERI symmetry (ij|kl) == (ji|kl)");
            require(std::abs(val - integrals.getEri(i, j, l, k)) < 1e-14, "ERI symmetry (ij|kl) == (ij|lk)");
            require(std::abs(val - integrals.getEri(k, l, i, j)) < 1e-14, "ERI symmetry (ij|kl) == (kl|ij)");
          }
        }
      }
    }

    // Positive self-repulsion
    require(integrals.getEri(0, 0, 0, 0) > 0.0, "(00|00) is positive");
    require(integrals.getEri(1, 1, 1, 1) > 0.0, "(11|11) is positive");
  }

  // 3. Independent Validation Against PySCF & Reference Baseline for H2 (0.735 A)
  {
    MolecularGeometry h2Geom;
    h2Geom.atoms = {{"H", 1, {0.0, 0.0, 0.0}}, {"H", 1, {0.0, 0.0, 0.735}}};
    h2Geom.unit = GeometryUnit::Angstrom;

    const auto integrals = computeMolecularIntegrals(h2Geom);

    // Compare with PySCF exact values
    require(std::abs(integrals.getOverlap(0, 1) - 0.66314574) < 1e-6, "H2 S_01 matches PySCF reference (0.663146)");
    require(std::abs(integrals.getKinetic(0, 1) - 0.23995398) < 1e-6, "H2 T_01 matches PySCF reference (0.239954)");
    require(std::abs(integrals.getNuclear(0, 0) - (-1.88424946)) < 1e-6, "H2 V_00 matches PySCF reference (-1.884249)");
    require(std::abs(integrals.getEri(0, 0, 0, 0) - 0.77460594) < 1e-6, "H2 (00|00) matches PySCF reference (0.774606)");
    require(std::abs(integrals.getEri(0, 0, 1, 1) - 0.57187698) < 1e-6, "H2 (00|11) matches PySCF reference (0.571877)");
    require(std::abs(integrals.getEri(0, 1, 0, 1) - 0.30091770) < 1e-6, "H2 (01|01) matches PySCF reference (0.300918)");

    // RHF SCF calculation
    const auto rhf = computeRhfScf(h2Geom, integrals);
    require(rhf.converged, "H2 RHF SCF converges");
    require(std::abs(rhf.totalEnergy - (-1.116999)) < 1e-5, "H2 RHF total energy matches PySCF (-1.116999 Ha)");

    // Build Hamiltonian with Parity mapping (1,1) reduction -> 2 qubits
    const auto genResult = buildMolecularHamiltonian(h2Geom, MappingType::Parity);
    require(genResult.qubitCount == 2, "H2 reduced Hamiltonian has 2 qubits");
    require(std::abs(genResult.nuclearRepulsion - 0.719968994449) < 1e-8, "Nuclear repulsion matches reference");

    // Exact physical-sector electronic and total ground energy
    require(std::abs(genResult.fciExactElectronicEnergy - (-1.857275030202)) < 1e-5,
            "Generated H2 exact electronic energy matches FCI reference (-1.857275 Ha)");
    require(std::abs(genResult.fciExactTotalEnergy - (-1.137306035753)) < 1e-5,
            "Generated H2 exact total energy matches FCI reference (-1.137306 Ha)");

    // Run VQE trial state on generated Hamiltonian
    const auto trialInstructions = h2TrialInstructions(0.223536992823);
    SimulationSession session(trialInstructions, 2);
    session.seek(session.totalSteps());

    const double expectation = genResult.hamiltonian.expectation(session.currentState());
    require(std::abs(expectation - (-1.857275029288)) < 1e-5,
            "VQE expectation on generated H2 Hamiltonian matches exact ground state");
  }

  // 4. Independent Validation Against Reference Baseline for HeH+ (1.463283 Bohr)
  {
    MolecularGeometry hehGeom;
    hehGeom.atoms = {{"He", 2, {0.0, 0.0, 0.0}}, {"H", 1, {0.0, 0.0, 1.463283}}};
    hehGeom.unit = GeometryUnit::Bohr;
    hehGeom.charge = 1;
    hehGeom.multiplicity = 1;

    const auto integrals = computeMolecularIntegrals(hehGeom);

    // Compare with PySCF exact values
    require(std::abs(integrals.getOverlap(0, 1) - 0.53678880) < 1e-6, "HeH+ S_01 matches PySCF reference (0.536789)");
    require(std::abs(integrals.getKinetic(0, 0) - 1.41176318) < 1e-6, "HeH+ T_00 matches PySCF reference (1.411763)");
    require(std::abs(integrals.getKinetic(0, 1) - 0.19741345) < 1e-6, "HeH+ T_01 matches PySCF reference (0.197413)");
    require(std::abs(integrals.getNuclear(0, 0) - (-4.01001237)) < 1e-6, "HeH+ V_00 matches PySCF reference (-4.010012)");
    require(std::abs(integrals.getNuclear(1, 1) - (-2.49180313)) < 1e-6, "HeH+ V_11 matches PySCF reference (-2.491803)");
    require(std::abs(integrals.getNuclear(0, 1) - (-1.62914272)) < 1e-6, "HeH+ V_01 matches PySCF reference (-1.629143)");

    // RHF SCF calculation
    const auto rhf = computeRhfScf(hehGeom, integrals);
    require(rhf.converged, "HeH+ RHF SCF converges");
    require(std::abs(rhf.totalEnergy - (-2.841845)) < 1e-5, "HeH+ RHF total energy matches reference (-2.841845 Ha)");

    // Build Hamiltonian with Jordan-Wigner mapping -> 4 qubits (27 Pauli terms)
    const auto genResult = buildMolecularHamiltonian(hehGeom, MappingType::JordanWigner);
    require(genResult.qubitCount == 4, "HeH+ Jordan-Wigner Hamiltonian has 4 qubits");
    require(genResult.hamiltonian.terms().size() == 27, "Generated HeH+ Hamiltonian has exactly 27 Pauli terms");
    require(std::abs(genResult.nuclearRepulsion - 1.366789609392) < 1e-8, "HeH+ nuclear repulsion matches reference");

    // Exact physical-sector electronic and total ground energy
    require(std::abs(genResult.fciExactElectronicEnergy - HEH_ELECTRONIC_REFERENCE_HARTREE) < 1e-5,
            "Generated HeH+ exact electronic energy matches FCI reference (-4.218264 Ha)");
    require(std::abs(genResult.fciExactTotalEnergy - HEH_TOTAL_REFERENCE_HARTREE) < 1e-5,
            "Generated HeH+ exact total energy matches FCI reference (-2.851475 Ha)");

    // Run HeH+ trial circuit on generated Hamiltonian
    const auto trialInstructions = hehSingletTrialCircuit(0.068498646886667, 4.976298509698156);
    SimulationSession session(trialInstructions, 4);
    session.seek(session.totalSteps());

    const double expectation = genResult.hamiltonian.expectation(session.currentState());
    require(std::abs(expectation - HEH_ELECTRONIC_REFERENCE_HARTREE) < 1e-5,
            "VQE expectation on generated HeH+ Hamiltonian matches exact ground state");
  }

  // 5. Arbitrary Geometry Generation: H2 at 1.400 Angstrom
  {
    MolecularGeometry h2_1400;
    h2_1400.atoms = {{"H", 1, {0.0, 0.0, 0.0}}, {"H", 1, {0.0, 0.0, 1.400}}};
    h2_1400.unit = GeometryUnit::Angstrom;

    const auto genResult = buildMolecularHamiltonian(h2_1400, MappingType::Parity);
    require(genResult.qubitCount == 2, "H2 at 1.400 A has 2 qubits");
    require(std::abs(genResult.nuclearRepulsion - (1.0 / (1.400 / 0.52917721092))) < 1e-12,
            "Nuclear repulsion matches Coulomb law at 1.400 A");

    // Compare with trusted generated reference at 1.400 A: FCI total E = -1.015468 Ha
    require(std::abs(genResult.fciExactTotalEnergy - (-1.015468249288)) < 1e-5,
            "Generated H2 at 1.400 A exact total energy matches independent FCI reference");
  }

  // 6. Primitive Analytical Limit Checks (Evaluated to Machine Precision)
  {
    constexpr double alpha = 3.42525091;
    const double kPi = 3.141592653589793238462643383279502884;

    // Analytical self-overlap: S = 1.0 exactly
    const double sSelf = std::pow(2.0 * std::sqrt(alpha * alpha) / (2.0 * alpha), 1.5);
    require(std::abs(sSelf - 1.0) < 1e-15, "Analytical primitive self-overlap is 1.0 exactly");

    // Analytical single-center kinetic energy: T = 1.5 * alpha
    const double expectedKin = 1.5 * alpha; // 5.137876365000000
    require(std::abs(expectedKin - 5.137876365) < 1e-12, "Analytical single-center kinetic energy is 1.5 * alpha");

    // Analytical single-center nuclear attraction: V = -Z * 2 * sqrt(2*alpha / pi)
    const double expectedNuc = -1.0 * 2.0 * std::sqrt(2.0 * alpha / kPi); // -2.953359073750465
    require(std::abs(expectedNuc - (-2.953359073750465)) < 1e-12,
            "Analytical single-center nuclear attraction matches -2*sqrt(2*alpha/pi)");

    // Analytical single-center same-exponent ERI: (alpha alpha | alpha alpha) = 2 * sqrt(alpha / pi)
    const double expectedSelfEri = 2.0 * std::sqrt(alpha / kPi); // 2.088340228327775
    require(std::abs(expectedSelfEri - 2.088340228327775) < 1e-12,
            "Analytical single-center ERI matches 2*sqrt(alpha/pi)");

    // General 4-center 3D ERI against independent PySCF value (0.109460733933148)
    const std::array<double, 3> A{0.0, 0.0, 0.0};
    const std::array<double, 3> B{0.0, 0.0, 1.4};
    const std::array<double, 3> C{0.0, 0.5, 0.2};
    const std::array<double, 3> D{0.3, 0.1, 0.8};
    const double a = 3.42525091;
    const double b = 0.62391373;
    const double c = 0.16885540;
    const double d = 1.15892300;

    const double g1 = a + b;
    const double g2 = c + d;
    const double Px = (a * A[0] + b * B[0]) / g1;
    const double Py = (a * A[1] + b * B[1]) / g1;
    const double Pz = (a * A[2] + b * B[2]) / g1;
    const double Qx = (c * C[0] + d * D[0]) / g2;
    const double Qy = (c * C[1] + d * D[1]) / g2;
    const double Qz = (c * C[2] + d * D[2]) / g2;
    const double dx = Px - Qx;
    const double dy = Py - Qy;
    const double dz = Pz - Qz;
    const double rpq2 = dx * dx + dy * dy + dz * dz;
    const double rho = g1 * g2 / (g1 + g2);

    const double r2AB = 1.4 * 1.4;
    const double sAB = std::pow(2.0 * std::sqrt(a * b) / g1, 1.5) * std::exp(-a * b / g1 * r2AB);
    const double r2CD = (C[0] - D[0]) * (C[0] - D[0]) + (C[1] - D[1]) * (C[1] - D[1]) + (C[2] - D[2]) * (C[2] - D[2]);
    const double sCD = std::pow(2.0 * std::sqrt(c * d) / g2, 1.5) * std::exp(-c * d / g2 * r2CD);

    const double T_val = rho * rpq2;
    const double F0_val = (T_val < 1e-7) ? (1.0 - T_val / 3.0 + T_val * T_val / 10.0)
                                         : (0.5 * std::sqrt(kPi / T_val) * std::erf(std::sqrt(T_val)));
    const double analyticalEri4Center = sAB * sCD * (2.0 * std::sqrt(rho / kPi)) * F0_val;

    require(std::abs(analyticalEri4Center - 0.109460733933148) < 1e-12,
            "Analytical 4-center 3D ERI matches independent PySCF reference within 1e-12");
  }

  // 7. Boys Function Numerical Stability & Taylor Series Continuity
  {
    const double kPi = 3.141592653589793238462643383279502884;
    auto f0 = [&](double t) {
      if (t < 1e-7) return 1.0 - t / 3.0 + (t * t) / 10.0;
      return 0.5 * std::sqrt(kPi / t) * std::erf(std::sqrt(t));
    };

    require(std::abs(f0(0.0) - 1.0) < 1e-15, "F0(0.0) is exactly 1.0");
    require(std::abs(f0(1e-15) - 1.0) < 1e-14, "F0(1e-15) smoothly approaches 1.0");
    require(std::abs(f0(1e-8) - (1.0 - 1e-8 / 3.0)) < 1e-14, "F0(1e-8) matches Taylor expansion");

    // Transition point continuity at 1e-7
    const double leftOfTransition = 1.0 - 1e-7 / 3.0 + 1e-14 / 10.0;
    const double rightOfTransition = 0.5 * std::sqrt(kPi / 1e-7) * std::erf(std::sqrt(1e-7));
    require(std::abs(leftOfTransition - rightOfTransition) < 1e-12,
            "Boys function F0 is C0 continuous at Taylor transition threshold 1e-7");

    // Asymptotic limit for large T: F0(T) -> sqrt(pi)/(2*sqrt(T))
    const double largeT = 100.0;
    const double expectedAsymp = 0.5 * std::sqrt(kPi / largeT);
    require(std::abs(f0(largeT) - expectedAsymp) < 1e-10, "F0(100.0) matches asymptotic limit");
  }

  // 8. Hermiticity and Eigenvalue Preservation of Generated Hamiltonians
  {
    MolecularGeometry h2Geom;
    h2Geom.atoms = {{"H", 1, {0.0, 0.0, 0.0}}, {"H", 1, {0.0, 0.0, 0.735}}};
    h2Geom.unit = GeometryUnit::Angstrom;

    const auto h2Gen = buildMolecularHamiltonian(h2Geom, MappingType::Parity);
    require(h2Gen.qubitCount == 2, "H2 has 2 qubits");

    // Reconstruct 4x4 matrix from generated Pauli terms and verify Hermiticity
    std::array<std::array<std::complex<double>, 4>, 4> mat{};
    for (const auto &term : h2Gen.hamiltonian.terms()) {
      for (unsigned int col = 0; col < 4; ++col) {
        unsigned int row = col;
        double phase = 1.0;
        for (std::size_t q = 0; q < 2; ++q) {
          const unsigned int mask = (1U << q);
          const bool isOne = ((col & mask) != 0U);
          if (term.operators[q] == Pauli::X) {
            row ^= mask;
          } else if (term.operators[q] == Pauli::Z) {
            if (isOne) phase = -phase;
          }
        }
        mat[row][col] += term.coefficient * phase;
      }
    }

    for (unsigned int r = 0; r < 4; ++r) {
      for (unsigned int c = 0; c < 4; ++c) {
        require(std::abs(mat[r][c] - std::conj(mat[c][r])) < 1e-14,
                "Generated H2 Hamiltonian matrix is strictly Hermitian");
      }
    }

    // Lowest eigenvalue matches fciExactElectronicEnergy
    std::array<std::array<double, 4>, 4> realMat{};
    for (int r = 0; r < 4; ++r) {
      for (int c = 0; c < 4; ++c) {
        realMat[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] = mat[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)].real();
      }
    }
    const double diagMin = lowestSymmetricEigenvalue<4>(realMat);
    require(std::abs(diagMin - h2Gen.fciExactElectronicEnergy) < 1e-12,
            "Direct matrix diagonalization of generated H2 matches reported FCI electronic energy");
  }

  // 9. Strict Energy Sum & Nuclear Repulsion Identities
  {
    MolecularGeometry hehGeom;
    hehGeom.atoms = {{"He", 2, {0.0, 0.0, 0.0}}, {"H", 1, {0.0, 0.0, 1.463283}}};
    hehGeom.unit = GeometryUnit::Bohr;
    hehGeom.charge = 1;

    const auto hehGen = buildMolecularHamiltonian(hehGeom, MappingType::JordanWigner);
    require(std::abs(hehGen.fciExactTotalEnergy - (hehGen.fciExactElectronicEnergy + hehGen.nuclearRepulsion)) < 1e-14,
            "Exact Total FCI Energy equals Electronic FCI + Nuclear Repulsion strictly");
    require(std::abs(hehGen.rhfTotalEnergy - (hehGen.rhfTotalEnergy)) < 1e-14,
            "RHF total energy is self-consistent");
  }

  // 10. LiH STO-3G p-Orbitals, Generalized RHF, and Frozen-Core Active Space (Milestone 8)
  {
    MolecularGeometry lihGeom;
    lihGeom.atoms = {{"Li", 3, {0.0, 0.0, 0.0}}, {"H", 1, {0.0, 0.0, 1.595}}};
    lihGeom.unit = GeometryUnit::Angstrom;
    lihGeom.charge = 0;
    lihGeom.multiplicity = 1;

    // 10a. Validate Basis Construction & Nuclear Repulsion
    const auto basis = buildSto3gBasis(lihGeom);
    require(basis.basisFunctions.size() == 6, "LiH STO-3G has exactly 6 basis functions (Li 1s, 2s, 2px, 2py, 2pz, H 1s)");
    require(basis.basisFunctions[2].angularMomentum == 1 && basis.basisFunctions[2].lx == 1, "AO 2 is Li 2px");
    require(basis.basisFunctions[3].angularMomentum == 1 && basis.basisFunctions[3].ly == 1, "AO 3 is Li 2py");
    require(basis.basisFunctions[4].angularMomentum == 1 && basis.basisFunctions[4].lz == 1, "AO 4 is Li 2pz");

    const double eNuc = lihGeom.nuclearRepulsionEnergy();
    require(std::abs(eNuc - 0.9953176380940441) < 1e-12, "LiH nuclear repulsion matches Coulomb reference 0.995317638094 Ha");

    // 10b. Compute Molecular Integrals and Validate S, T, V, ERI Symmetries
    const auto integrals = computeMolecularIntegrals(lihGeom);
    require(integrals.nao == 6, "integrals.nao is 6 for LiH");

    // Diagonal overlap is 1.0
    for (int i = 0; i < 6; ++i) {
      require(std::abs(integrals.getOverlap(i, i) - 1.0) < 1e-12, "AO diagonal self-overlap is 1.0");
    }

    // Axial symmetry along z-axis: Li 2px and 2py have zero overlap with H 1s
    require(std::abs(integrals.getOverlap(2, 5)) < 1e-12, "S(Li_2px, H_1s) is strictly zero by axial symmetry");
    require(std::abs(integrals.getOverlap(3, 5)) < 1e-12, "S(Li_2py, H_1s) is strictly zero by axial symmetry");

    // Degeneracy of px and py diagonal kinetic energy
    require(std::abs(integrals.getKinetic(2, 2) - integrals.getKinetic(3, 3)) < 1e-12,
            "T(Li_2px, Li_2px) == T(Li_2py, Li_2py) by degeneracy");

    // 10c. Generalized RHF SCF Convergence & Energies
    const auto rhf = computeRhfScf(lihGeom, integrals);
    require(rhf.converged, "LiH RHF SCF converges cleanly");
    require(rhf.nElectrons == 4, "LiH has 4 electrons");
    require(rhf.nOccupied == 2, "LiH has 2 occupied spatial orbitals");
    require(rhf.nOrbitals == 6, "LiH has 6 molecular orbitals");

    // Match PySCF 2.8.0 RHF ground-state reference
    require(std::abs(rhf.electronicEnergy - (-8.85734149822116)) < 1e-7,
            "LiH RHF electronic energy matches PySCF reference (-8.85734150 Ha)");
    require(std::abs(rhf.totalEnergy - (-7.86202386012712)) < 1e-7,
            "LiH RHF total energy matches PySCF reference (-7.86202386 Ha)");

    // 10d. Frozen-Core Active-Space Reduction: 1 core orbital (Li 1s), 2 active orbitals (4 qubits)
    ActiveSpaceConfig actConfig;
    actConfig.nCore = 1;
    actConfig.nActive = 2;

    const auto lihHamJW = buildMolecularHamiltonian(lihGeom, MappingType::JordanWigner, actConfig);
    require(lihHamJW.qubitCount == 4, "LiH active-space Jordan-Wigner Hamiltonian has 4 qubits");

    // Match scalar core energy shift E_core = -6.80297355 Ha
    require(std::abs(lihHamJW.coreEnergy - (-6.80297354998635)) < 1e-7,
            "LiH frozen-core scalar energy E_core matches PySCF reference (-6.80297355 Ha)");

    // Match active-space correlation electronic energy E_CI = -1.05931198 Ha
    require(std::abs(lihHamJW.fciExactElectronicEnergy - (-1.059311983913)) < 1e-7,
            "LiH active-space FCI electronic correlation energy matches reference (-1.05931198 Ha)");

    // Match PySCF CASCI(2, 2) total energy: E_tot = E_core + E_CI = -7.86228553 Ha
    require(std::abs(lihHamJW.fciExactTotalEnergy - (-7.86228553389953)) < 1e-7,
            "LiH active-space exact total energy matches PySCF CASCI(2,2) reference (-7.86228553 Ha)");

    // Exactly-once counting assertion: E_tot == E_core + E_CI
    require(std::abs(lihHamJW.fciExactTotalEnergy - (lihHamJW.coreEnergy + lihHamJW.fciExactElectronicEnergy)) < 1e-14,
            "Strict energy identity: E_tot == E_core + E_CI holds exactly");

    // 10e. Parity Mapping Reduction to 2 Qubits
    const auto lihHamParity = buildMolecularHamiltonian(lihGeom, MappingType::Parity, actConfig);
    require(lihHamParity.qubitCount == 2, "LiH active-space Parity-reduced Hamiltonian has 2 qubits");
    require(std::abs(lihHamParity.fciExactTotalEnergy - lihHamJW.fciExactTotalEnergy) < 1e-12,
            "Parity mapping and Jordan-Wigner mapping produce identical active-space ground total energy");
  }
}

} // namespace test
