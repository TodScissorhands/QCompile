#include "chemistry/benchmark_catalog.h"
#include "chemistry/h2_energy.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {

void printUsage() {
  std::cerr << "Usage: ./qcompile_benchmarks --list | --json-list | "
               "[--json] <benchmark-id> | --generated-h2 <file.qham> | "
               "--generated-h2-json <file.qham>\n";
}

void printJsonString(std::string_view value) {
  std::cout << '"';
  for (char character : value) {
    switch (character) {
    case '"':
      std::cout << "\\\"";
      break;
    case '\\':
      std::cout << "\\\\";
      break;
    case '\n':
      std::cout << "\\n";
      break;
    case '\r':
      std::cout << "\\r";
      break;
    case '\t':
      std::cout << "\\t";
      break;
    default:
      std::cout << character;
    }
  }
  std::cout << '"';
}

void printMetadata(const BenchmarkMetadata &metadata) {
  std::cout << "{\"id\":";
  printJsonString(metadata.id);
  std::cout << ",\"display_name\":";
  printJsonString(metadata.displayName);
  std::cout << ",\"formula\":";
  printJsonString(metadata.formula);
  std::cout << ",\"geometry\":";
  printJsonString(metadata.geometry);
  std::cout << ",\"geometry_units\":";
  printJsonString(metadata.geometryUnits);
  std::cout << ",\"charge\":" << metadata.charge
            << ",\"electron_count\":" << metadata.electronCount
            << ",\"total_spin\":" << metadata.totalSpin
            << ",\"spin_multiplicity\":" << metadata.spinMultiplicity
            << ",\"basis\":";
  printJsonString(metadata.basis);
  std::cout << ",\"spatial_orbitals\":" << metadata.spatialOrbitals
            << ",\"spin_orbitals\":" << metadata.spinOrbitals
            << ",\"qubits\":" << metadata.qubits
            << ",\"mapping\":";
  printJsonString(metadata.mapping);
  std::cout << ",\"physical_sector\":";
  printJsonString(metadata.physicalSector);
  std::cout << ",\"qubit_ordering\":";
  printJsonString(metadata.qubitOrdering);
  std::cout << ",\"reference_provenance\":";
  printJsonString(metadata.referenceProvenance);
  std::cout << ",\"reference_electronic_energy_hartree\":"
            << metadata.referenceElectronicEnergy
            << ",\"nuclear_repulsion_hartree\":"
            << metadata.nuclearRepulsionEnergy
            << ",\"energy_unit\":";
  printJsonString(metadata.energyUnit);
  std::cout << ",\"nuclear_repulsion_separate\":"
            << (metadata.nuclearRepulsionSeparate ? "true" : "false") << '}';
}

void printJsonResult(const BenchmarkDescriptor &benchmark,
                     const BenchmarkResult &result) {
  std::cout << std::setprecision(17) << "{\"benchmark\":";
  printMetadata(benchmark.metadata);
  std::cout << ",\"electronic_energy_hartree\":" << result.electronicEnergy
            << ",\"exact_electronic_energy_hartree\":"
            << result.exactElectronicEnergy
            << ",\"reference_electronic_energy_hartree\":"
            << result.referenceEnergy
            << ",\"variational_error_hartree\":" << result.variationalError
            << ",\"nuclear_repulsion_hartree\":"
            << result.nuclearRepulsionEnergy << ",\"total_energy_hartree\":"
            << result.totalEnergy << ",\"energy_evaluations\":"
            << result.energyEvaluations << ",\"converged\":"
            << (result.converged ? "true" : "false") << ",\"parameters\":[";
  for (std::size_t index = 0; index < result.parameterCount; ++index) {
    if (index != 0) {
      std::cout << ',';
    }
    std::cout << result.parameters[index];
  }
  std::cout << "]}\n";
}


void printGeneratedH2Result(const GeneratedH2 &model,
                            const GeneratedH2EnergyResult &result) {
  std::cout << std::fixed << std::setprecision(15)
            << "Generated H2 geometry-driven benchmark\n"
            << "Geometry: H2, R = " << model.distanceAngstrom << " Angstrom ("
            << model.distanceBohr << " bohr)\n"
            << "Model: neutral singlet, 2 electrons, STO-3G, RHF, parity (1,1), "
               "2 qubits\n"
            << "Optimized circuit parameter (radians): " << result.parameter
            << " rad\n"
            << "Calculated electronic energy: " << result.electronicEnergy
            << " Hartree\n"
            << "Physical-sector exact electronic energy: "
            << result.exactElectronicEnergy << " Hartree\n"
            << "Generated FCI reference electronic energy: "
            << result.referenceElectronicEnergy << " Hartree\n"
            << "Variational error (calculated - exact): "
            << result.variationalError << " Hartree\n"
            << "Nuclear repulsion (separate): "
            << result.nuclearRepulsionEnergy << " Hartree\n"
            << "Calculated total Born-Oppenheimer energy: "
            << result.totalEnergy << " Hartree\n"
            << "Exact total Born-Oppenheimer energy: "
            << result.exactTotalEnergy << " Hartree\n"
            << "Optimizer evaluations: " << result.energyEvaluations << '\n'
            << "Converged: " << (result.converged ? "yes" : "no") << '\n';
}
void printGeneratedH2Json(const GeneratedH2 &model,
                          const GeneratedH2EnergyResult &result) {
  std::cout << std::setprecision(17)
            << "{\"schema\":\"QCOMPILE-GENERATED-H2-RESULT-1\""
            << ",\"benchmark\":\"generated_h2\",\"formula\":\"H2\""
            << ",\"distance_angstrom\":" << model.distanceAngstrom
            << ",\"distance_bohr\":" << model.distanceBohr
            << ",\"geometry_units\":\"Angstrom\""
            << ",\"charge\":0,\"electron_count\":2,\"total_spin\":0"
            << ",\"multiplicity\":1,\"basis\":\"STO-3G\""
            << ",\"spatial_orbitals\":2,\"qubits\":2"
            << ",\"mapping\":\"parity\",\"particle_number_reduction\":[1,1]"
            << ",\"qubit_ordering\":\"q0 is least-significant state-vector bit\""
            << ",\"energy_unit\":\"Hartree\""
            << ",\"electronic_energy_excludes_nuclear_repulsion\":true"
            << ",\"electronic_energy_hartree\":" << result.electronicEnergy
            << ",\"exact_electronic_energy_hartree\":" << result.exactElectronicEnergy
            << ",\"reference_electronic_energy_hartree\":" << result.referenceElectronicEnergy
            << ",\"nuclear_repulsion_hartree\":" << result.nuclearRepulsionEnergy
            << ",\"total_energy_hartree\":" << result.totalEnergy
            << ",\"exact_total_energy_hartree\":" << result.exactTotalEnergy
            << ",\"variational_error_hartree\":" << result.variationalError
            << ",\"optimizer_parameter_radians\":" << result.parameter
            << ",\"evaluation_count\":" << result.energyEvaluations
            << ",\"converged\":" << (result.converged ? "true" : "false")
            << "}\n";
}

int runGeneratedH2(const char *path, bool jsonOutput) {
  try {
    const GeneratedH2 model = loadGeneratedH2(path);
    const GeneratedH2EnergyResult result = optimizeGeneratedH2(model);
    if (jsonOutput) {
      printGeneratedH2Json(model, result);
    } else {
      printGeneratedH2Result(model, result);
    }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}


} // namespace

int main(int argc, char *argv[]) {
  if (argc == 3 && std::string_view(argv[1]) == "--generated-h2") {
    return runGeneratedH2(argv[2], false);
  }
  if (argc == 3 && std::string_view(argv[1]) == "--generated-h2-json") {
    return runGeneratedH2(argv[2], true);
  }
  if ((argc == 2 && std::string_view(argv[1]) == "--generated-h2") ||
      (argc == 2 && std::string_view(argv[1]) == "--generated-h2-json")) {
    printUsage();
    return 2;
  }
  if (argc == 2 && std::string_view(argv[1]) == "--list") {
    for (const BenchmarkDescriptor &benchmark : benchmarkCatalog()) {
      std::cout << benchmark.metadata.id << "\t"
                << benchmark.metadata.displayName << '\n';
    }
    return 0;
  }
  if (argc == 2 && std::string_view(argv[1]) == "--json-list") {
    std::cout << std::setprecision(17) << "{\"benchmarks\":[";
    bool first = true;
    for (const BenchmarkDescriptor &benchmark : benchmarkCatalog()) {
      if (!first) {
        std::cout << ',';
      }
      printMetadata(benchmark.metadata);
      first = false;
    }
    std::cout << "]}\n";
    return 0;
  }

  const bool jsonOutput =
      argc == 3 && std::string_view(argv[1]) == "--json";
  if ((!jsonOutput && argc != 2) || (jsonOutput && argc != 3)) {
    printUsage();
    return 2;
  }
  const std::string_view benchmarkId =
      jsonOutput ? std::string_view(argv[2]) : std::string_view(argv[1]);
  const BenchmarkDescriptor *benchmark = findBenchmark(benchmarkId);
  if (benchmark == nullptr) {
    std::cerr << "Unknown benchmark ID: " << benchmarkId << '\n';
    printUsage();
    return 2;
  }

  try {
    const BenchmarkResult result = runBenchmark(benchmark->metadata.id);
    if (jsonOutput) {
      printJsonResult(*benchmark, result);
      return 0;
    }
    std::cout << std::fixed << std::setprecision(15)
              << benchmark->metadata.displayName << '\n'
              << "Electronic energy: " << result.electronicEnergy << ' '
              << benchmark->metadata.energyUnit << '\n'
              << "Exact electronic energy: " << result.exactElectronicEnergy
              << ' ' << benchmark->metadata.energyUnit << '\n'
              << "Reference electronic energy: " << result.referenceEnergy
              << ' ' << benchmark->metadata.energyUnit << '\n'
              << "Variational error: " << result.variationalError << ' '
              << benchmark->metadata.energyUnit << '\n'
              << "Nuclear repulsion: " << result.nuclearRepulsionEnergy << ' '
              << benchmark->metadata.energyUnit << '\n'
              << "Total Born-Oppenheimer energy: " << result.totalEnergy << ' '
              << benchmark->metadata.energyUnit << '\n'
              << "Energy evaluations: " << result.energyEvaluations << '\n'
              << "Converged: " << (result.converged ? "yes" : "no") << '\n'
              << "Parameters:";
    for (std::size_t index = 0; index < result.parameterCount; ++index) {
      std::cout << ' ' << result.parameters[index];
    }
    std::cout << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}

