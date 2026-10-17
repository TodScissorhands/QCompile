#pragma once

#include "chemistry/hamiltonian.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

enum class GeometryUnit {
  Angstrom,
  Bohr
};

struct MolecularAtom {
  std::string element; // "H", "He", "Li", "Be"
  int atomicNumber = 0;
  std::array<double, 3> coordinates{0.0, 0.0, 0.0};
};

struct MolecularGeometry {
  std::vector<MolecularAtom> atoms;
  GeometryUnit unit = GeometryUnit::Angstrom;
  int charge = 0;
  int multiplicity = 1; // 2S + 1 (1 for singlet)

  void validate() const;
  double nuclearRepulsionEnergy() const; // in Hartree
  std::vector<MolecularAtom> inBohr() const;
};

struct GaussianPrimitive {
  double exponent = 0.0;
  double contractionCoefficient = 0.0;
};

struct ContractedGaussian {
  int atomIndex = 0;
  int angularMomentum = 0; // 0 for s
  int lx = 0, ly = 0, lz = 0;
  std::array<double, 3> centerBohr{0.0, 0.0, 0.0};
  std::vector<GaussianPrimitive> primitives;
};

struct Sto3gBasis {
  std::vector<ContractedGaussian> basisFunctions;
};

Sto3gBasis buildSto3gBasis(const MolecularGeometry &geometry);

struct MolecularIntegrals {
  int nao = 0;
  std::vector<double> overlap; // nao x nao (row-major)
  std::vector<double> kinetic; // nao x nao
  std::vector<double> nuclear; // nao x nao
  std::vector<double> hCore;   // nao x nao (T + V)
  std::vector<double> eri;     // nao x nao x nao x nao (chemist: (uv|ls))
  double nuclearRepulsion = 0.0;

  double getOverlap(int i, int j) const { return overlap[static_cast<std::size_t>(i * nao + j)]; }
  double getKinetic(int i, int j) const { return kinetic[static_cast<std::size_t>(i * nao + j)]; }
  double getNuclear(int i, int j) const { return nuclear[static_cast<std::size_t>(i * nao + j)]; }
  double getHCore(int i, int j) const { return hCore[static_cast<std::size_t>(i * nao + j)]; }
  double getEri(int i, int j, int k, int l) const {
    return eri[static_cast<std::size_t>(((i * nao + j) * nao + k) * nao + l)];
  }
};

MolecularIntegrals computeMolecularIntegrals(const MolecularGeometry &geometry);

struct RhfResult {
  bool converged = false;
  int iterations = 0;
  double electronicEnergy = 0.0;
  double totalEnergy = 0.0;
  int nElectrons = 0;
  int nOccupied = 0;
  int nOrbitals = 0;
  std::vector<double> moCoefficients; // nao x nmo
  std::vector<double> moEnergies;     // nmo
  std::vector<double> h1Mo;           // nmo x nmo (MO 1-electron integrals)
  std::vector<double> eriMo;          // nmo x nmo x nmo x nmo (MO 2-electron integrals)
};

RhfResult computeRhfScf(const MolecularGeometry &geometry, const MolecularIntegrals &integrals);

enum class MappingType {
  Parity,
  JordanWigner
};

struct ActiveSpaceConfig {
  int nCore = 0;   // Number of frozen core spatial orbitals (0 for H2/HeH+; 1 for LiH)
  int nActive = 2; // Number of active spatial orbitals (default 2)
};

struct GeneratedHamiltonianResult {
  QubitHamiltonian hamiltonian = QubitHamiltonian(0, {});
  int qubitCount = 0;
  double nuclearRepulsion = 0.0;
  double coreEnergy = 0.0;              // E_core = E_nuc + core electronic shift
  double rhfTotalEnergy = 0.0;
  double fciExactElectronicEnergy = 0.0; // E_CI in active space
  double fciExactTotalEnergy = 0.0;      // E_tot = E_core + E_CI
  int activeElectrons = 2;
  int activeSpatialOrbitals = 2;
  std::string provenance;
};

GeneratedHamiltonianResult buildMolecularHamiltonian(
    const MolecularGeometry &geometry,
    MappingType mapping = MappingType::Parity,
    const ActiveSpaceConfig &activeSpace = {}
);
