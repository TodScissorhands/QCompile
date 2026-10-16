#pragma once

#include <array>
#include <cstddef>
#include <string_view>

struct BenchmarkMetadata {
  std::string_view id;
  std::string_view displayName;
  std::string_view formula;
  std::string_view geometry;
  std::string_view geometryUnits;
  int charge;
  int electronCount;
  int totalSpin;
  int spinMultiplicity;
  std::string_view basis;
  int spatialOrbitals;
  int spinOrbitals;
  int qubits;
  std::string_view mapping;
  std::string_view physicalSector;
  std::string_view qubitOrdering;
  std::string_view referenceProvenance;
  double referenceElectronicEnergy;
  double nuclearRepulsionEnergy;
  std::string_view energyUnit;
  bool nuclearRepulsionSeparate;
};

struct BenchmarkResult {
  std::array<double, 2> parameters{};
  std::size_t parameterCount = 0;
  double electronicEnergy = 0.0;
  double nuclearRepulsionEnergy = 0.0;
  double totalEnergy = 0.0;
  double exactElectronicEnergy = 0.0;
  double referenceEnergy = 0.0;
  double variationalError = 0.0;
  int energyEvaluations = 0;
  bool converged = false;
};

using BenchmarkRunner = BenchmarkResult (*)();

struct BenchmarkDescriptor {
  BenchmarkMetadata metadata;
  BenchmarkRunner run;
};

const std::array<BenchmarkDescriptor, 2> &benchmarkCatalog();
const BenchmarkDescriptor *findBenchmark(std::string_view id);
BenchmarkResult runBenchmark(std::string_view id);
