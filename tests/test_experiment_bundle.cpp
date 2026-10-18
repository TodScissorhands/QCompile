#include "test_support.h"
#include "experiment/experiment_bundle.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace test {

void testExperimentBundle() {
  // 1. Full Round-trip Serialization & Deserialization
  {
    ExperimentBundle original;
    original.schema = "QCOMPILE-EXPERIMENT-1";
    original.id = "exp-bell-001";
    original.title = "Bell State Entangled Pair";
    original.timestamp = "2026-10-17T12:00:00Z";
    original.description = "Full simulation and chemistry test experiment";

    original.dslSource = "qubit q0, q1\nH q0\nCNOT q0, q1\nmeasure q0, q1\n";
    original.optimize = true;
    original.diagnostics = {"Info: compilation clean"};

    original.qubitCount = 2;
    original.rawInstructions = {
        {InstructionType::ApplyGate, "H", {0}, 0.0},
        {InstructionType::ApplyGate, "CNOT", {0, 1}, 0.0},
        {InstructionType::Measure, "measure", {0, 1}, 0.0}
    };
    original.instructions = original.rawInstructions;

    original.currentStep = 3;
    original.totalSteps = 3;
    original.isComplete = true;
    original.stateNorm = 1.0;
    original.statevector = {
        std::complex<double>(0.7071067811865475, 0.0),
        std::complex<double>(0.0, 0.0),
        std::complex<double>(0.0, 0.0),
        std::complex<double>(0.7071067811865475, 0.0)
    };
    original.basisProbabilities = {0.5, 0.0, 0.0, 0.5};
    original.recordedMeasurements = {{3, "11"}};

    // Chemistry metadata
    original.hasChemistry = true;
    original.moleculeName = "HeH+";
    original.formula = "HeH+";
    original.basis = "STO-3G";
    original.mapping = "jordan_wigner";
    original.distance = 0.774331;
    original.distanceUnits = "Angstrom";
    original.electronicEnergy = -4.218264154490334;
    original.nuclearRepulsion = 1.366789609392032;
    original.totalEnergy = -2.851474545098302;
    original.exactElectronicEnergy = -4.218264154490395;
    original.exactTotalEnergy = -2.851474545098363;
    original.variationalError = 0.000000000000061;
    original.energyEvaluations = 6803;
    original.converged = true;
    original.optimizerParameterRadians = 0.068498646886667;
    original.referenceElectronicEnergy = -4.218264154490395;

    // Real convergence trajectory sample
    original.hasConvergence = true;
    original.convergenceTrajectory = {
        {1, 0.0, 0.0, -4.208634703392453, -2.841845094000421, true},
        {2, 0.05, 4.5, -4.215000000000000, -2.848210390607968, true},
        {3, 0.068498646886667, 4.976298509698156, -4.218264154490334, -2.851474545098302, true}
    };

    // Serialize to JSON string
    const std::string jsonText = serializeExperimentBundle(original);
    require(!jsonText.empty(), "serialized JSON string is non-empty");

    // Deserialize back
    const ExperimentBundle loaded = deserializeExperimentBundle(jsonText);

    // Verify metadata
    require(loaded.schema == original.schema, "schema preserved");
    require(loaded.id == original.id, "id preserved");
    require(loaded.title == original.title, "title preserved");
    require(loaded.timestamp == original.timestamp, "timestamp preserved");
    require(loaded.description == original.description, "description preserved");

    // Verify source
    require(loaded.dslSource == original.dslSource, "DSL source preserved");
    require(loaded.optimize == original.optimize, "optimize flag preserved");
    require(loaded.diagnostics == original.diagnostics, "diagnostics preserved");

    // Verify circuit
    require(loaded.qubitCount == original.qubitCount, "qubit count preserved");
    require(loaded.rawInstructions.size() == original.rawInstructions.size(), "raw instructions size preserved");
    require(loaded.instructions.size() == original.instructions.size(), "optimized instructions size preserved");
    require(loaded.instructions[1].gateName == "CNOT", "instruction gate preserved");
    require(loaded.instructions[1].qubitIndices == std::vector<int>({0, 1}), "CNOT operands preserved");

    // Verify simulation
    require(loaded.currentStep == original.currentStep, "current step preserved");
    require(loaded.totalSteps == original.totalSteps, "total steps preserved");
    require(loaded.isComplete == original.isComplete, "completion flag preserved");
    require(std::abs(loaded.stateNorm - original.stateNorm) < 1e-12, "state norm preserved");
    require(loaded.statevector.size() == original.statevector.size(), "statevector size preserved");
    require(std::abs(loaded.statevector[0].real() - original.statevector[0].real()) < 1e-15, "state amplitude real preserved");
    require(loaded.basisProbabilities == original.basisProbabilities, "basis probabilities preserved");
    require(loaded.recordedMeasurements == original.recordedMeasurements, "recorded measurements preserved");

    // Verify chemistry numerical fidelity (17 digits)
    require(loaded.hasChemistry, "hasChemistry preserved");
    require(loaded.formula == original.formula, "formula preserved");
    require(std::abs(loaded.distance - original.distance) < 1e-12, "distance preserved");
    require(std::abs(loaded.electronicEnergy - original.electronicEnergy) < 1e-15, "electronic energy preserved at full precision");
    require(std::abs(loaded.nuclearRepulsion - original.nuclearRepulsion) < 1e-15, "nuclear repulsion preserved");
    require(std::abs(loaded.totalEnergy - original.totalEnergy) < 1e-15, "total energy preserved");
    require(std::abs(loaded.exactTotalEnergy - original.exactTotalEnergy) < 1e-15, "exact total energy preserved");
    require(loaded.energyEvaluations == original.energyEvaluations, "evaluations preserved");
    require(loaded.optimizerParameterRadians.has_value(), "optimizer parameter is present");
    require(std::abs(*loaded.optimizerParameterRadians - *original.optimizerParameterRadians) < 1e-15,
            "optimizer parameter preserved");
    require(loaded.referenceElectronicEnergy.has_value(), "reference electronic energy is present");
    require(std::abs(*loaded.referenceElectronicEnergy - *original.referenceElectronicEnergy) < 1e-15,
            "reference electronic energy preserved");

    require(loaded.converged == original.converged, "converged preserved");

    // Verify convergence trajectory
    require(loaded.hasConvergence, "hasConvergence preserved");
    require(loaded.convergenceTrajectory.size() == 3, "trajectory count preserved");
    require(loaded.convergenceTrajectory[2].evaluation == 3, "eval index preserved");
    require(std::abs(loaded.convergenceTrajectory[2].totalEnergy - original.convergenceTrajectory[2].totalEnergy) < 1e-15,
            "convergence point energy preserved");
    require(loaded.convergenceTrajectory[2].isBest, "isBest flag preserved");
  }
  // Version-1 bundles written before the optional result fields remain loadable.
  {
    const std::string legacyJson =
        R"({"schema":"QCOMPILE-EXPERIMENT-1","source":{"dsl":"","optimize":false},"circuit":{"qubit_count":0},"chemistry":{"has_chemistry":true,"formula":"H2"}})";
    const ExperimentBundle legacy = deserializeExperimentBundle(legacyJson);
    require(legacy.hasChemistry && legacy.formula == "H2", "legacy chemistry metadata remains loadable");
    require(!legacy.optimizerParameterRadians.has_value(), "legacy bundle does not invent an optimizer parameter");
    require(!legacy.referenceElectronicEnergy.has_value(), "legacy bundle does not invent a reference energy");
  }

  // 2. Schema Validation, Malformed Handling & Boundary Rejection
  {
    // Empty JSON
    bool emptyRejected = false;
    try {
      (void)deserializeExperimentBundle("");
    } catch (const std::exception &) {
      emptyRejected = true;
    }
    require(emptyRejected, "empty JSON rejected");

    // Malformed JSON (syntax error)
    bool malformedRejected = false;
    try {
      (void)deserializeExperimentBundle("{ \"schema\": \"QCOMPILE-EXPERIMENT-1\", ");
    } catch (const std::exception &) {
      malformedRejected = true;
    }
    require(malformedRejected, "truncated/malformed JSON rejected");

    // Unsupported schema version
    bool wrongSchemaRejected = false;
    try {
      (void)deserializeExperimentBundle("{\"schema\": \"QCOMPILE-EXPERIMENT-99\", \"source\": {\"dsl\": \"\"}, \"circuit\": {\"qubit_count\": 0}}");
    } catch (const std::runtime_error &err) {
      wrongSchemaRejected = (std::string(err.what()).find("Unsupported .qcomp bundle schema") != std::string::npos);
    }
    require(wrongSchemaRejected, "unsupported schema version rejected with clear error");

    // Missing required fields
    bool missingDslRejected = false;
    try {
      (void)deserializeExperimentBundle("{\"schema\": \"QCOMPILE-EXPERIMENT-1\", \"circuit\": {\"qubit_count\": 1}}");
    } catch (const std::exception &) {
      missingDslRejected = true;
    }
    require(missingDslRejected, "missing source.dsl field rejected");

    // Out of range qubit count
    bool excessiveQubitsRejected = false;
    try {
      (void)deserializeExperimentBundle("{\"schema\": \"QCOMPILE-EXPERIMENT-1\", \"source\": {\"dsl\": \"\"}, \"circuit\": {\"qubit_count\": 25}}");
    } catch (const std::exception &) {
      excessiveQubitsRejected = true;
    }
    require(excessiveQubitsRejected, "qubit_count > 24 in .qcomp bundle rejected");
  }

  // 3. Atomic File Save and Load
  {
    const std::string tempFile = std::filesystem::temp_directory_path() / "qcompile_test_bundle_atomic.qcomp";

    ExperimentBundle bundle;
    bundle.id = "atomic-test";
    bundle.title = "Atomic Test File";
    bundle.dslSource = "qubit q0\nH q0\n";
    bundle.qubitCount = 1;
    bundle.instructions = {{InstructionType::ApplyGate, "H", {0}, 0.0}};

    // Save atomically
    saveExperimentBundle(tempFile, bundle);
    require(std::filesystem::exists(tempFile), "saved .qcomp file exists on disk");

    // Load back
    const ExperimentBundle loaded = loadExperimentBundle(tempFile);
    require(loaded.id == "atomic-test", "loaded bundle id matches");
    require(loaded.qubitCount == 1, "loaded qubit count is 1");
    require(loaded.instructions.size() == 1, "loaded instruction count is 1");

    // Re-save (atomic replace of existing file)
    bundle.title = "Updated Atomic Title";
    saveExperimentBundle(tempFile, bundle);
    const ExperimentBundle updated = loadExperimentBundle(tempFile);
    require(updated.title == "Updated Atomic Title", "atomic replace succeeded");

    // Invalid in-memory data must be rejected before the existing destination is replaced.
    std::ifstream savedFile(tempFile, std::ios::binary);
    const std::string savedBytes((std::istreambuf_iterator<char>(savedFile)), std::istreambuf_iterator<char>());
    const std::string savedTitle = updated.title;
    auto rejectedSavePreservesDestination = [&](ExperimentBundle invalid, const std::string &description) {
      bool rejected = false;
      try {
        saveExperimentBundle(tempFile, invalid);
      } catch (const std::exception &) {
        rejected = true;
      }
      std::ifstream afterFile(tempFile, std::ios::binary);
      const std::string afterBytes((std::istreambuf_iterator<char>(afterFile)), std::istreambuf_iterator<char>());
      require(rejected, description + " save rejected");
      require(afterBytes == savedBytes, description + " save preserved destination bytes");
      require(loadExperimentBundle(tempFile).title == savedTitle, description + " destination remains loadable");
    };

    ExperimentBundle invalid = bundle;
    invalid.qubitCount = 25;
    rejectedSavePreservesDestination(invalid, "unsupported qubit count");
    invalid = bundle;
    invalid.schema = "QCOMPILE-EXPERIMENT-99";
    rejectedSavePreservesDestination(invalid, "unsupported schema");
    invalid = bundle;
    invalid.stateNorm = std::numeric_limits<double>::quiet_NaN();
    rejectedSavePreservesDestination(invalid, "non-finite simulation state");
    invalid = bundle;
    invalid.instructions[0].parameterValue = std::numeric_limits<double>::infinity();
    rejectedSavePreservesDestination(invalid, "non-finite circuit parameter");
    invalid = bundle;
    invalid.hasChemistry = true;
    invalid.electronicEnergy = std::numeric_limits<double>::infinity();
    rejectedSavePreservesDestination(invalid, "non-finite chemistry result");
    invalid = bundle;
    invalid.hasConvergence = true;
    invalid.convergenceTrajectory = {{1, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN(), -1.0, true}};
    rejectedSavePreservesDestination(invalid, "non-finite convergence point");
    invalid = bundle;
    invalid.instructions[0].qubitIndices = {1};
    rejectedSavePreservesDestination(invalid, "out-of-range circuit qubit");
    invalid = bundle;
    invalid.statevector = {std::complex<double>(1.0, 0.0)};
    rejectedSavePreservesDestination(invalid, "statevector dimension mismatch");
    invalid = bundle;
    invalid.currentStep = 1;
    invalid.totalSteps = 1;
    invalid.isComplete = false;
    invalid.statevector = {std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)};
    invalid.basisProbabilities = {1.0, 0.0};
    rejectedSavePreservesDestination(invalid, "inconsistent session completion state");
    ExperimentBundle snapshot = bundle;
    snapshot.totalSteps = 1;
    snapshot.statevector = {std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)};
    snapshot.basisProbabilities = {1.0, 0.0};
    invalid = snapshot;
    invalid.stateNorm = 0.5;
    rejectedSavePreservesDestination(invalid, "state norm inconsistent with statevector");
    invalid = snapshot;
    invalid.basisProbabilities = {0.5, 0.5};
    rejectedSavePreservesDestination(invalid, "basis probabilities inconsistent with statevector");
    invalid = snapshot;
    invalid.statevector = {std::complex<double>(2.0, 0.0), std::complex<double>(0.0, 0.0)};
    invalid.stateNorm = 2.0;
    invalid.basisProbabilities = {4.0, 0.0};
    rejectedSavePreservesDestination(invalid, "unnormalized statevector");
    invalid = snapshot;
    invalid.basisProbabilities = {-1.0, 2.0};
    rejectedSavePreservesDestination(invalid, "negative basis probability");
    invalid = snapshot;
    invalid.instructions = {{InstructionType::Measure, "measure", {0}, 0.0}};
    invalid.currentStep = 1;
    invalid.isComplete = true;
    invalid.recordedMeasurements = {{1, "0"}, {1, "1"}};
    rejectedSavePreservesDestination(invalid, "duplicate recorded measurement step");
    // Clean up temp file
    std::filesystem::remove(tempFile);
  }

  // 4. Cross-Run Comparison Inspector
  {
    // Compatible runs: H2 at two different bond lengths
    ExperimentBundle h2_0735;
    h2_0735.title = "H2 at 0.735 A";
    h2_0735.hasChemistry = true;
    h2_0735.formula = "H2";
    h2_0735.basis = "STO-3G";
    h2_0735.distance = 0.735;
    h2_0735.distanceUnits = "Angstrom";
    h2_0735.totalEnergy = -1.137306035753;
    h2_0735.energyEvaluations = 62;
    h2_0735.qubitCount = 2;
    h2_0735.instructions = {{InstructionType::ApplyGate, "H", {0}, 0.0}};
    h2_0735.stateNorm = 1.0;

    ExperimentBundle h2_1400;
    h2_1400.title = "H2 at 1.400 A";
    h2_1400.hasChemistry = true;
    h2_1400.formula = "H2";
    h2_1400.basis = "STO-3G";
    h2_1400.distance = 1.400;
    h2_1400.distanceUnits = "Angstrom";
    h2_1400.totalEnergy = -1.015468249288;
    h2_1400.energyEvaluations = 62;
    h2_1400.qubitCount = 2;
    h2_1400.instructions = {{InstructionType::ApplyGate, "H", {0}, 0.0}};
    h2_1400.stateNorm = 1.0;

    const auto compSame = compareExperiments(h2_0735, h2_1400);
    require(compSame.moleculesMatch, "molecules match for H2 comparison");
    require(!compSame.geometriesMatch, "geometries differ as expected");
    require(compSame.isDirectlyComparable, "PES points of same molecule are directly comparable");
    require(!compSame.items.empty(), "comparison items generated");

    // Incompatible runs: H2 vs HeH+
    ExperimentBundle heh;
    heh.title = "HeH+ at 0.774 A";
    heh.hasChemistry = true;
    heh.formula = "HeH+";
    heh.basis = "STO-3G";
    heh.distance = 0.774331;
    heh.distanceUnits = "Angstrom";
    heh.totalEnergy = -2.851474545098;
    heh.qubitCount = 4;
    heh.stateNorm = 1.0;

    const auto compDiff = compareExperiments(h2_0735, heh);
    require(!compDiff.moleculesMatch, "molecules differ between H2 and HeH+");
    require(!compDiff.isDirectlyComparable, "different molecules are not directly comparable");
    require(compDiff.summary.find("Mismatched molecules") != std::string::npos, "summary explains mismatch");
  }

  // 5. Circuit SVG Vector Export
  {
    std::vector<Instruction> circuit = {
        {InstructionType::ApplyGate, "H", {0}, 0.0},
        {InstructionType::ApplyGate, "CNOT", {0, 1}, 0.0},
        {InstructionType::Measure, "measure", {0, 1}, 0.0}
    };

    const std::string svg = exportCircuitSvg(circuit, 2);
    require(svg.find("<?xml") != std::string::npos, "SVG header present");
    require(svg.find("<svg") != std::string::npos, "SVG root tag present");
    require(svg.find("</svg>") != std::string::npos, "SVG closing tag present");
    require(svg.find("q0") != std::string::npos && svg.find("q1") != std::string::npos, "qubit wire labels present in SVG");
    require(svg.find("H") != std::string::npos, "H gate present in SVG");
    require(svg.find("<circle") != std::string::npos, "CNOT circle present in SVG");
    require(svg.find("viewBox=") != std::string::npos, "viewBox present in SVG");

    // Test save to disk
    const std::string tempSvg = std::filesystem::temp_directory_path() / "qcompile_circuit_test.svg";
    saveCircuitSvg(tempSvg, circuit, 2);
    require(std::filesystem::exists(tempSvg), "saved SVG file exists");
    std::filesystem::remove(tempSvg);
  }

  // 6. Convergence CSV Tabular Export
  {
    ExperimentBundle bundle;
    bundle.hasChemistry = true;
    bundle.nuclearRepulsion = 0.719968994449;
    bundle.hasConvergence = true;
    bundle.convergenceTrajectory = {
        {1, 0.2, 0.0, -1.85, -1.130031005551, true},
        {2, 0.223536992823, 0.0, -1.857275029288, -1.137306034839, true}
    };

    const std::string csv = exportConvergenceCsv(bundle);
    require(csv.find("evaluation,parameter0,parameter1,electronic_energy_hartree") != std::string::npos,
            "CSV header contains expected columns");
    require(csv.find("-1.857275029288") != std::string::npos, "high precision electronic energy present in CSV");
    require(csv.find("-1.137306034839") != std::string::npos, "high precision total energy present in CSV");
    require(csv.find("0.7199689944") != std::string::npos, "nuclear repulsion present in CSV");

    // Test save to disk
    const std::string tempCsv = std::filesystem::temp_directory_path() / "qcompile_conv_test.csv";
    saveConvergenceCsv(tempCsv, bundle);
    require(std::filesystem::exists(tempCsv), "saved CSV file exists");

    // Parse back CSV to verify rows and format
    std::ifstream in(tempCsv);
    std::string headerLine, line1, line2;
    std::getline(in, headerLine);
    std::getline(in, line1);
    std::getline(in, line2);
    require(!headerLine.empty(), "CSV header line read");
    require(!line1.empty(), "CSV row 1 read");
    require(!line2.empty(), "CSV row 2 read");
    require(line2.rfind("true") != std::string::npos, "is_best_so_far column preserved");

    in.close();
    std::filesystem::remove(tempCsv);
  }
}

} // namespace test
