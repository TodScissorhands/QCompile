#include "chemistry/benchmark_catalog.h"

#include "chemistry/h2_energy.h"
#include "chemistry/heh_energy.h"

#include <stdexcept>

namespace {

BenchmarkResult runH2() {
  const H2Result result = optimizeH2();
  return {
      {result.parameter, 0.0},
      1,
      result.electronicEnergy,
      result.nuclearRepulsionEnergy,
      result.totalEnergy,
      result.exactElectronicEnergy,
      H2_ELECTRONIC_REFERENCE_HARTREE,
      result.electronicEnergy - result.exactElectronicEnergy,
      result.energyEvaluations,
      true
  };
}

BenchmarkResult runHeHPlus() {
  const HeHResult result = optimizeHeH();
  return {
      {result.theta, result.phi},
      2,
      result.electronicEnergy,
      result.nuclearRepulsionEnergy,
      result.totalEnergy,
      result.exactElectronicEnergy,
      HEH_ELECTRONIC_REFERENCE_HARTREE,
      result.electronicEnergy - result.exactElectronicEnergy,
      result.energyEvaluations,
      true
  };
}

const std::array<BenchmarkDescriptor, 2> kBenchmarks{{
    {{
         "h2",
         "H₂ (fixed STO-3G benchmark)",
         "H2",
         "H 0 0 0; H 0 0 0.735",
         "angstrom",
         0,
         2,
         0,
         1,
         "STO-3G",
         2,
         4,
         2,
         "two-qubit parity mapping; particle-number reduced (1,1)",
         "exact minimum of the reduced 4x4 Hamiltonian; trial remains in odd-parity block",
         "QCompile qubit 0 is least-significant state-vector bit",
         "Qiskit Nature electronic-structure and qubit-mapping tutorials",
         H2_ELECTRONIC_REFERENCE_HARTREE,
         H2_NUCLEAR_REPULSION_HARTREE,
         "Hartree",
         true
     }, runH2},
    {{
         "heh_plus",
         "HeH+ (fixed STO-3G benchmark)",
         "HeH+",
         "He 0 0 0; H 0 0 1.463283",
         "bohr",
         1,
         2,
         0,
         1,
         "STO-3G",
         2,
         4,
         4,
         "Jordan-Wigner; no symmetry reduction",
         "three singlet CSFs within the four-determinant Nalpha=Nbeta=1, Ms=0 space",
         "q0=MO0 alpha, q1=MO1 alpha, q2=MO0 beta, q3=MO1 beta; q0 is LSB",
         "benchmarks/heh_plus/reference_data.json and independent Psi4 FCI",
         HEH_ELECTRONIC_REFERENCE_HARTREE,
         HEH_NUCLEAR_REPULSION_HARTREE,
         "Hartree",
         true
     }, runHeHPlus}
}};

} // namespace

const std::array<BenchmarkDescriptor, 2> &benchmarkCatalog() {
  return kBenchmarks;
}

const BenchmarkDescriptor *findBenchmark(std::string_view id) {
  for (const BenchmarkDescriptor &benchmark : benchmarkCatalog()) {
    if (benchmark.metadata.id == id) {
      return &benchmark;
    }
  }
  return nullptr;
}

BenchmarkResult runBenchmark(std::string_view id) {
  const BenchmarkDescriptor *benchmark = findBenchmark(id);
  if (benchmark == nullptr) {
    throw std::invalid_argument("Unknown fixed benchmark ID");
  }
  return benchmark->run();
}
