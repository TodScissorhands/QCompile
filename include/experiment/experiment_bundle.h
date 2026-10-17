#pragma once

#include "compiler/code_generator.h"
#include "chemistry/vqe_convergence.h"

#include <complex>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

struct ExperimentBundle {
  // Schema and metadata
  std::string schema = "QCOMPILE-EXPERIMENT-1";
  std::string id;
  std::string title;
  std::string timestamp;
  std::string description;

  // DSL source and compiler configuration
  std::string dslSource;
  bool optimize = true;
  std::vector<std::string> diagnostics;

  // Circuit model
  int qubitCount = 0;
  std::vector<Instruction> rawInstructions;
  std::vector<Instruction> instructions;

  // Simulation state and measurement snapshots
  std::size_t currentStep = 0;
  std::size_t totalSteps = 0;
  bool isComplete = false;
  double stateNorm = 1.0;
  std::vector<std::complex<double>> statevector;
  std::vector<double> basisProbabilities;
  std::vector<std::pair<std::size_t, std::string>> recordedMeasurements; // (step, outcomeString)

  // Chemistry model and results (optional)
  bool hasChemistry = false;
  std::string moleculeName;
  std::string formula;
  std::string basis;
  std::string mapping;
  double distance = 0.0;
  std::string distanceUnits;
  double electronicEnergy = 0.0;
  double nuclearRepulsion = 0.0;
  double totalEnergy = 0.0;
  double exactElectronicEnergy = 0.0;
  double exactTotalEnergy = 0.0;
  double variationalError = 0.0;
  int energyEvaluations = 0;
  bool converged = false;

  // Real optimizer convergence trajectory (optional)
  bool hasConvergence = false;
  std::vector<VqeConvergencePoint> convergenceTrajectory;
};

// Serialization & Deserialization (pure C++17)
std::string serializeExperimentBundle(const ExperimentBundle &bundle);
ExperimentBundle deserializeExperimentBundle(const std::string &jsonText);

// Atomic file operations (writes to .tmp and renames)
void saveExperimentBundle(const std::string &filePath, const ExperimentBundle &bundle);
ExperimentBundle loadExperimentBundle(const std::string &filePath);

// Tabular Convergence Export (CSV with 17-digit numerical precision and escaping)
std::string exportConvergenceCsv(const ExperimentBundle &bundle);
void saveConvergenceCsv(const std::string &filePath, const ExperimentBundle &bundle);

// Circuit Vector Export (Scalable Vector Graphics SVG)
std::string exportCircuitSvg(const std::vector<Instruction> &instructions, int qubitCount);
void saveCircuitSvg(const std::string &filePath, const std::vector<Instruction> &instructions, int qubitCount);

// Cross-Run Comparison Inspector
struct ExperimentComparisonItem {
  std::string propertyName;
  std::string valueA;
  std::string valueB;
  std::string delta;
  bool isCompatible = true;
  std::string note;
};

struct ExperimentComparisonResult {
  std::string titleA;
  std::string titleB;
  std::string formulaA;
  std::string formulaB;
  bool moleculesMatch = true;
  bool geometriesMatch = true;
  bool isDirectlyComparable = true;
  std::string summary;
  std::vector<ExperimentComparisonItem> items;
};

ExperimentComparisonResult compareExperiments(const ExperimentBundle &a, const ExperimentBundle &b);
