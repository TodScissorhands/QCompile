#include "test_support.h"
#include "chemistry/molecular_integral_engine.h"
#include "chemistry/h2_energy.h"
#include "chemistry/heh_energy.h"
#include "simulator/simulation_session.h"

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
}

} // namespace test
