#include "compiler/code_generator.h"
#include "chemistry/benchmark_catalog.h"
#include "chemistry/h2_energy.h"
#include "chemistry/symmetric_eigensolver.h"
#include "chemistry/heh_energy.h"
#include "chemistry/hamiltonian.h"
#include "compiler/lexer.h"
#include "compiler/parser.h"
#include "compiler/semantic_checker.h"
#include "simulator/simulator.h"
#include <algorithm>
#include <utility>
#include <atomic>
#include <filesystem>
#include <fstream>

#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <system_error>

namespace {

int assertions = 0;

void require(bool condition, const std::string &message) {
  ++assertions;
  if (!condition) {
    throw std::runtime_error(message);
  }
}

Program parse(const std::string &source) {
  std::vector<Token> tokens = Lexer(source).tokenize();
  return Parser(tokens).parse();
}

SemanticResult check(const std::string &source, Program &program) {
  program = parse(source);
  return checkSemantics(program);
}

std::string run(const std::string &source) {
  Program program = parse(source);
  SemanticResult semantic = checkSemantics(program);
  require(semantic.errors.empty(), "simulator input must pass semantic checking");
  std::vector<Instruction> instructions = generateInstructions(program, semantic);
  std::ostringstream output;
  std::streambuf *original = std::cout.rdbuf(output.rdbuf());
  try {
    Simulator simulator(static_cast<int>(semantic.qubitNames.size()));
    simulator.run(instructions);
  } catch (...) {
    std::cout.rdbuf(original);
    throw;
  }
  std::cout.rdbuf(original);
  return output.str();
}

using Complex = std::complex<double>;
using Matrix2 = std::array<std::array<Complex, 2>, 2>;
using Matrix4 = std::array<std::array<Complex, 4>, 4>;

std::string generatedH2Fixture(double angstrom) {
  if (angstrom < 1.0) {
    return
        "QCOMPILE-H2-PAULI 1\n"
        "distance_angstrom 0.735\n"
        "distance_bohr 1.3889487015553204\n"
        "charge 0\n"
        "electrons 2\n"
        "spin 0\n"
        "multiplicity 1\n"
        "basis STO-3G\n"
        "energy_unit Hartree\n"
        "spatial_orbitals 2\n"
        "spin_orbitals 4\n"
        "qubits 2\n"
        "mapping PARITY\n"
        "particle_counts_alpha_beta 1 1\n"
        "spin_orbital_order alpha-block-then-beta-block\n"
        "parity_reduction fixed_p1=1_fixed_p3=0\n"
        "pauli_label_order q1q0_rightmost_is_q0\n"
        "qcompile_operator_order q0_to_q1_q0_is_lsb\n"
        "electronic_energy_convention nuclear_repulsion_excluded\n"
        "nuclear_repulsion_hartree 0.71996899444897966\n"
        "singlet_fci_electronic_reference_hartree -1.8572750302023802\n"
        "term_count 5\n"
        "TERMS\n"
        "II -1.0523732457728596\n"
        "IZ 0.39793742484317934\n"
        "XX 0.18093119978423144\n"
        "ZI -0.39793742484317934\n"
        "ZZ -0.011280104256235268\n";
  }
  return
      "QCOMPILE-H2-PAULI 1\n"
      "distance_angstrom 1.4\n"
      "distance_bohr 2.645616574391086\n"
      "charge 0\n"
      "electrons 2\n"
      "spin 0\n"
      "multiplicity 1\n"
      "basis STO-3G\n"
      "energy_unit Hartree\n"
      "spatial_orbitals 2\n"
      "spin_orbitals 4\n"
      "qubits 2\n"
      "mapping PARITY\n"
      "particle_counts_alpha_beta 1 1\n"
      "spin_orbital_order alpha-block-then-beta-block\n"
      "parity_reduction fixed_p1=1_fixed_p3=0\n"
      "pauli_label_order q1q0_rightmost_is_q0\n"
      "qcompile_operator_order q0_to_q1_q0_is_lsb\n"
      "electronic_energy_convention nuclear_repulsion_excluded\n"
      "nuclear_repulsion_hartree 0.37798372208571435\n"
      "singlet_fci_electronic_reference_hartree -1.3934519713739597\n"
      "term_count 5\n"
      "TERMS\n"
      "II -1.0253590356223656\n"
      "IZ 0.14956793876884486\n"
      "XX 0.22302208907471485\n"
      "ZI -0.1495679387688449\n"
      "ZZ -0.005030536366542548\n";
}

std::filesystem::path generatedH2TempPath() {
  static std::atomic<unsigned int> nextId{0};
  for (;;) {
    const auto candidate = std::filesystem::temp_directory_path() /
        ("qcompile-generated-h2-test-" +
         std::to_string(nextId.fetch_add(1)) + ".qham");
    if (!std::filesystem::exists(candidate)) {
      return candidate;
    }
  }
}

GeneratedH2 loadGeneratedH2Text(const std::string &text) {
  const auto path = generatedH2TempPath();
  {
    std::ofstream file(path);
    if (!file) {
      throw std::runtime_error("could not create temporary generated-H2 fixture");
    }
    file << text;
  }
  try {
    GeneratedH2 result = loadGeneratedH2(path.string());
    std::filesystem::remove(path);
    return result;
  } catch (...) {
    std::filesystem::remove(path);
    throw;
  }
}

bool generatedH2TextRejected(const std::string &text) {
  try {
    (void)loadGeneratedH2Text(text);
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

std::string replaceOnce(std::string text, const std::string &from,
                        const std::string &to) {
  const std::size_t position = text.find(from);
  if (position == std::string::npos) {
    throw std::runtime_error("generated-H2 test fixture field not found");
  }
  text.replace(position, from.size(), to);
  return text;
}

Matrix2 pauliMatrix(Pauli pauli) {
  Matrix2 matrix{};
  if (pauli == Pauli::I) {
    matrix[0][0] = matrix[1][1] = 1.0;
  } else if (pauli == Pauli::X) {
    matrix[0][1] = matrix[1][0] = 1.0;
  } else if (pauli == Pauli::Z) {
    matrix[0][0] = 1.0;
    matrix[1][1] = -1.0;
  } else if (pauli == Pauli::Y) {
    matrix[0][1] = Complex(0.0, -1.0);
    matrix[1][0] = Complex(0.0, 1.0);
  }
  return matrix;
}

Matrix4 kron(const Matrix2 &left, const Matrix2 &right) {
  Matrix4 result{};
  for (int row = 0; row < 4; ++row) {
    for (int column = 0; column < 4; ++column) {
      result[row][column] =
          left[row / 2][column / 2] * right[row % 2][column % 2];
    }
  }
  return result;
}

double denseExpectation(
    const Matrix4 &matrix,
    const std::array<Complex, 4> &state
) {
  Complex result(0.0, 0.0);
  for (int row = 0; row < 4; ++row) {
    for (int column = 0; column < 4; ++column) {
      result += std::conj(state[row]) * matrix[row][column] * state[column];
    }
  }
  require(std::abs(result.imag()) < 1e-12, "dense Hermitian expectation is real");
  return result.real();
}

void testPauliExpectations() {
  const std::vector<Complex> zero{1.0, 0.0};
  const std::vector<Complex> one{0.0, 1.0};
  require(std::abs(pauliExpectation(zero, {Pauli::I}) - 1.0) < 1e-12,
          "identity expectation on |0> is one");
  require(std::abs(pauliExpectation(zero, {Pauli::X})) < 1e-12,
          "X expectation on |0> is zero");
  require(std::abs(pauliExpectation(zero, {Pauli::Z}) - 1.0) < 1e-12,
          "Z expectation on |0> is one");
  require(std::abs(pauliExpectation(one, {Pauli::Z}) + 1.0) < 1e-12,
          "Z expectation on |1> is negative one");

  const double inverseSqrtTwo = 1.0 / std::sqrt(2.0);
  const std::vector<Complex> plus{inverseSqrtTwo, inverseSqrtTwo};
  require(std::abs(pauliExpectation(plus, {Pauli::X}) - 1.0) < 1e-12,
          "X expectation on |+> is one");
  require(std::abs(pauliExpectation(plus, {Pauli::Z})) < 1e-12,
          "Z expectation on |+> is zero");
  const std::vector<Complex> plusI{
      inverseSqrtTwo, Complex(0.0, inverseSqrtTwo)
  };
  require(std::abs(pauliExpectation(plusI, {Pauli::Y}) - 1.0) < 1e-12,
          "Y expectation on |+i> is one");

  const std::vector<Complex> bell{inverseSqrtTwo, 0.0, 0.0, inverseSqrtTwo};
  require(std::abs(pauliExpectation(bell, {Pauli::X, Pauli::X}) - 1.0) <
              1e-12,
          "XX expectation on Bell state is one");
  require(std::abs(pauliExpectation(bell, {Pauli::Z, Pauli::Z}) - 1.0) <
              1e-12,
          "ZZ expectation on Bell state is one");
  require(std::abs(pauliExpectation(bell, {Pauli::I, Pauli::I}) - 1.0) <
              1e-12,
          "two-qubit identity expectation is one");

  bool rejected = false;
  try {
    pauliExpectation({1.0, 0.0}, {Pauli::I, Pauli::I});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "state dimension must match Pauli width");
  rejected = false;
  try {
    pauliExpectation({1.0, 1.0}, {Pauli::I});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "unnormalized state is rejected");
  rejected = false;
  try {
    pauliExpectation({Complex(NAN, 0.0), 0.0}, {Pauli::I});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "non-finite state amplitude is rejected");
}

void testHamiltonianExpectations() {
  const std::vector<Complex> raw{
      Complex(1.0, 2.0), Complex(-2.0, 0.5),
      Complex(0.3, -0.8), Complex(1.1, 0.4)
  };
  double normSquared = 0.0;
  for (const Complex &amplitude : raw) {
    normSquared += std::norm(amplitude);
  }
  std::array<Complex, 4> state{};
  std::vector<Complex> stateVector;
  for (int index = 0; index < 4; ++index) {
    state[index] = raw[index] / std::sqrt(normSquared);
    stateVector.push_back(state[index]);
  }

  const QubitHamiltonian hamiltonian = h2Hamiltonian();
  double termwise = 0.0;
  Matrix4 dense{};
  for (const PauliTerm &term : hamiltonian.terms()) {
    termwise += term.coefficient * pauliExpectation(stateVector, term.operators);
    const Matrix4 product = kron(
        pauliMatrix(term.operators[1]), pauliMatrix(term.operators[0])
    );
    for (int row = 0; row < 4; ++row) {
      for (int column = 0; column < 4; ++column) {
        dense[row][column] += term.coefficient * product[row][column];
      }
    }
  }
  require(std::abs(termwise - denseExpectation(dense, state)) < 1e-12,
          "termwise Hamiltonian expectation matches dense matrix evaluation");

  Simulator simulator(2);
  simulator.run({
      {InstructionType::ApplyGate, "H", {0}, 0.0},
      {InstructionType::ApplyGate, "CNOT", {0, 1}, 0.0}
  });
  require(std::abs(simulator.expectation({Pauli::X, Pauli::X}) - 1.0) <
              1e-12,
          "simulator expectation reads the current Bell state");
  require(std::abs(hamiltonian.expectation(simulator) - (-0.88272215)) <
              1e-8,
          "Hamiltonian expectation evaluates weighted terms on simulator state");

  bool rejected = false;
  try {
    simulator.expectation({Pauli::I});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "simulator rejects Pauli width mismatch");
  rejected = false;
  try {
    QubitHamiltonian invalid(2, {{1.0, {Pauli::Z}}});
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "Hamiltonian rejects wrong-width Pauli term");
}

void testSymmetricEigensolverValidation() {
  std::array<std::array<double, 2>, 2> matrix{{
      {{2.0, 1.0}},
      {{1.0, 2.0}}
  }};
  require(std::abs(lowestSymmetricEigenvalue(matrix) - 1.0) < 1e-12,
          "Jacobi solver returns the lowest eigenvalue of a symmetric matrix");

  auto rejects = [](const std::array<std::array<double, 2>, 2> &input) {
    try {
      lowestSymmetricEigenvalue(input);
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  auto invalid = matrix;
  invalid[0][0] = std::numeric_limits<double>::quiet_NaN();
  require(rejects(invalid), "Jacobi solver rejects NaN matrix entries");
  invalid = matrix;
  invalid[1][1] = std::numeric_limits<double>::infinity();
  require(rejects(invalid), "Jacobi solver rejects infinite matrix entries");
  invalid = matrix;
  invalid[0][1] += 1e-8;
  require(rejects(invalid), "Jacobi solver rejects asymmetric matrices");
}

void testBenchmarkCatalog() {
  const auto &catalog = benchmarkCatalog();
  require(catalog.size() == 2, "catalog contains only the two fixed benchmarks");
  const BenchmarkDescriptor *h2 = findBenchmark("h2");
  const BenchmarkDescriptor *heh = findBenchmark("heh_plus");
  require(h2 != nullptr && heh != nullptr,
          "catalog resolves both stable benchmark IDs");
  require(findBenchmark("lih") == nullptr,
          "catalog does not imply support for unregistered molecules");
  require(h2->metadata.formula == "H2" &&
              h2->metadata.geometryUnits == "angstrom" &&
              h2->metadata.electronCount == 2 &&
              h2->metadata.qubits == 2 &&
              h2->metadata.physicalSector.find("odd-parity") !=
                  std::string_view::npos &&
              h2->metadata.referenceElectronicEnergy ==
                  H2_ELECTRONIC_REFERENCE_HARTREE,
          "H2 metadata describes its fixed reduced parity model and sector");
  require(heh->metadata.formula == "HeH+" &&
              heh->metadata.geometryUnits == "bohr" &&
              heh->metadata.charge == 1 &&
              heh->metadata.electronCount == 2 &&
              heh->metadata.spinMultiplicity == 1 &&
              heh->metadata.qubits == 4 &&
              heh->metadata.physicalSector.find("three singlet CSFs") !=
                  std::string_view::npos &&
              heh->metadata.referenceElectronicEnergy ==
                  HEH_ELECTRONIC_REFERENCE_HARTREE,
          "HeH+ metadata describes its fixed singlet Jordan-Wigner model");
  require(h2->metadata.nuclearRepulsionSeparate &&
              heh->metadata.nuclearRepulsionSeparate &&
              h2->metadata.energyUnit == "Hartree" &&
              heh->metadata.energyUnit == "Hartree",
          "catalog defines electronic-only Hamiltonians and Hartree units");

  const BenchmarkResult h2Result = runBenchmark("h2");
  require(h2Result.converged && h2Result.parameterCount == 1 &&
              h2Result.energyEvaluations > 0,
          "H2 catalog runner reports successful convergence and one parameter");
  require(std::abs(h2Result.variationalError -
                   (h2Result.electronicEnergy -
                    h2Result.exactElectronicEnergy)) < 1e-14 &&
              h2Result.referenceEnergy == H2_ELECTRONIC_REFERENCE_HARTREE,
          "H2 common result preserves exact and external reference energies");
  require(std::abs(h2Result.totalEnergy -
                   (h2Result.electronicEnergy +
                    h2Result.nuclearRepulsionEnergy)) < 1e-14,
          "H2 common result reports nuclear repulsion separately");

  const BenchmarkResult hehResult = runBenchmark("heh_plus");
  require(hehResult.converged && hehResult.parameterCount == 2 &&
              hehResult.energyEvaluations > 0,
          "HeH+ catalog runner reports convergence and two parameters");
  require(std::abs(hehResult.variationalError -
                   (hehResult.electronicEnergy -
                    hehResult.exactElectronicEnergy)) < 1e-14 &&
              hehResult.referenceEnergy == HEH_ELECTRONIC_REFERENCE_HARTREE,
          "HeH+ common result preserves exact and external reference energies");
  require(std::abs(hehResult.totalEnergy -
                   (hehResult.electronicEnergy +
                    hehResult.nuclearRepulsionEnergy)) < 1e-14,
          "HeH+ common result reports nuclear repulsion separately");

  bool rejected = false;
  try {
    runBenchmark("lih");
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "running an unknown benchmark ID is rejected");
}

void testH2Energy() {
  const double exact = exactH2GroundEnergy();
  require(std::abs(exact - (-1.8572750292882292)) < 1e-12,
          "independent diagonalization returns rounded-Hamiltonian minimum");
  require(std::abs(exact - H2_ELECTRONIC_REFERENCE_HARTREE) < 1e-8,
          "rounded Hamiltonian matches published electronic reference");

  const double atZero = evaluateH2Energy(0.0);
  const auto trialInstructions = h2TrialInstructions(0.5);
  require(trialInstructions.size() == 8 &&
              trialInstructions[0].gateName == "X" &&
              trialInstructions[1].gateName == "CNOT" &&
              trialInstructions[4].gateName == "RZ" &&
              trialInstructions[4].qubitIndices == std::vector<int>{0} &&
              trialInstructions[4].parameterValue == 0.5 &&
              trialInstructions[7].gateName == "CNOT",
          "exposed H2 trial instructions preserve the executed gate sequence and parameter");
  const double atHalf = evaluateH2Energy(0.5);
  require(std::abs(atZero - atHalf) > 1e-3,
          "changing trial parameter changes the energy");
  const double halfSine = std::sin(0.25);
  const double halfCosine = std::cos(0.25);
  const double expectedAtHalf =
      -1.83696799 * halfCosine * halfCosine -
      0.24521831 * halfSine * halfSine -
      2.0 * 0.18093120 * halfSine * halfCosine;
  require(std::abs(atZero - (-1.83696799)) < 1e-12,
          "zero parameter prepares the encoded reference state");
  require(std::abs(atHalf - expectedAtHalf) < 1e-12,
          "trial circuit prepares the analytically expected rotated state");
  bool rejected = false;
  try {
    evaluateH2Energy(std::numeric_limits<double>::quiet_NaN());
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, "non-finite variational parameter is rejected");

  const H2Result first = optimizeH2();
  const H2Result second = optimizeH2();
  require(first.parameter == second.parameter &&
              first.electronicEnergy == second.electronicEnergy &&
              first.energyEvaluations == second.energyEvaluations,
          "bounded golden-section optimizer is deterministic");
  require(std::abs(first.electronicEnergy - exact) < 1e-10,
          "optimized energy matches independent exact minimum");
  require(std::abs(first.totalEnergy -
                   (first.electronicEnergy + H2_NUCLEAR_REPULSION_HARTREE)) <
              1e-14,
          "total energy adds nuclear repulsion exactly once");
  require(std::abs(first.totalEnergy - H2_TOTAL_REFERENCE_HARTREE) < 1e-8,
          "total energy matches published reference");
  require(first.nuclearRepulsionEnergy == H2_NUCLEAR_REPULSION_HARTREE,
          "electronic and nuclear energies remain separately reported");
  require(std::abs(first.parameter - 0.22353700115632177) < 1e-7,
          "optimized parameter is near the exact ground-state mixing angle");
}

void testGeneratedH2() {
  const std::string fixture0735 = generatedH2Fixture(0.735);
  const std::string fixture1400 = generatedH2Fixture(1.4);
  const GeneratedH2 model0735 = loadGeneratedH2Text(fixture0735);
  const GeneratedH2 model1400 = loadGeneratedH2Text(fixture1400);
  require(std::abs(model0735.distanceAngstrom - 0.735) < 1e-15 &&
              std::abs(model1400.distanceAngstrom - 1.4) < 1e-15,
          "generated H2 loader accepts both geometry fixtures");
  require(model0735.electronicHamiltonian.qubitCount() == 2 &&
              model0735.electronicHamiltonian.terms().size() == 5,
          "generated H2 loader preserves two-qubit Pauli Hamiltonian");

  bool foundQ0Z = false;
  for (const PauliTerm &term : model0735.electronicHamiltonian.terms()) {
    if (term.operators == std::vector<Pauli>{Pauli::Z, Pauli::I}) {
      foundQ0Z = std::abs(term.coefficient - 0.39793742484317934) < 1e-15;
    }
  }
  require(foundQ0Z, "Qiskit q1q0 labels reverse into q0-as-LSB operator order");

  const double exact0735 = exactGeneratedH2GroundEnergy(model0735);
  const double exact1400 = exactGeneratedH2GroundEnergy(model1400);
  require(std::abs(exact0735 - model0735.referenceElectronicEnergy) < 1e-10 &&
              std::abs(exact1400 - model1400.referenceElectronicEnergy) < 1e-10,
          "generated physical-sector exact energies match FCI references");
  require(std::abs(exact0735 - exact1400) > 0.4,
          "generated H2 energies respond to the requested geometry");
  require(std::abs(evaluateGeneratedH2Energy(model0735, 0.0) -
                   evaluateGeneratedH2Energy(model0735, 0.5)) > 1e-3,
          "generated-Hamiltonian energy evaluation uses the supplied coefficients");

  const GeneratedH2EnergyResult result0735 = optimizeGeneratedH2(model0735);
  const GeneratedH2EnergyResult result1400 = optimizeGeneratedH2(model1400);
  for (const GeneratedH2EnergyResult *result : {&result0735, &result1400}) {
    require(result->converged && result->energyEvaluations > 0,
            "generated H2 optimizer reports convergence and evaluations");
    require(result->electronicEnergy >=
                result->exactElectronicEnergy - 1e-10 &&
                std::abs(result->electronicEnergy -
                         result->exactElectronicEnergy) < 1e-8,
            "generated H2 variational energy matches its sector minimum");
    require(std::abs(result->totalEnergy -
                     (result->electronicEnergy +
                      result->nuclearRepulsionEnergy)) < 1e-14 &&
                std::abs(result->exactTotalEnergy -
                         (result->exactElectronicEnergy +
                          result->nuclearRepulsionEnergy)) < 1e-14,
            "generated H2 total energies add nuclear repulsion exactly once");
    require(std::abs(result->variationalError -
                     (result->electronicEnergy -
                      result->exactElectronicEnergy)) < 1e-14,
            "generated H2 variational error is calculated minus exact");
  }

  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "QCOMPILE-H2-PAULI 1", "QCOMPILE-H2-PAULI 2")),
          "generated H2 loader rejects unknown schema versions");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "basis STO-3G\n", "")),
          "generated H2 loader rejects missing metadata");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "distance_angstrom 0.735", "distance_angstrom 0")),
          "generated H2 loader rejects invalid geometry");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "charge 0", "charge 1")),
          "generated H2 loader rejects incompatible charge");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "spin 0", "spin 1")),
          "generated H2 loader rejects incompatible spin");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "basis STO-3G", "basis 6-31G")),
          "generated H2 loader rejects incompatible basis");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "mapping PARITY", "mapping JORDAN_WIGNER")),
          "generated H2 loader rejects unsupported mapping");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "qcompile_operator_order q0_to_q1_q0_is_lsb",
              "qcompile_operator_order q1_to_q0")),
          "generated H2 loader rejects incompatible qubit order");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "term_count 5", "term_count 6")),
          "generated H2 loader rejects unsupported term counts");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "ZI -0.39793742484317934\n", "")),
          "generated H2 loader rejects a missing declared Pauli term");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "distance_bohr 1.3889487015553204",
              "distance_bohr 1.4")),
          "generated H2 loader rejects inconsistent distance units");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "qubits 2", "qubits 4")),
          "generated H2 loader rejects incompatible qubit counts");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "ZI -0.39793742484317934",
              "XX -0.39793742484317934")),
          "generated H2 loader rejects duplicate terms");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "II -1.0523732457728596",
              "QQ -1.0523732457728596")),
          "generated H2 loader rejects unknown Pauli symbols");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "II -1.0523732457728596",
              "III -1.0523732457728596")),
          "generated H2 loader rejects invalid Pauli widths");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "II -1.0523732457728596", "II nan")),
          "generated H2 loader rejects non-finite coefficients");
  require(generatedH2TextRejected(replaceOnce(
              fixture0735, "nuclear_repulsion_hartree 0.71996899444897966",
              "nuclear_repulsion_hartree 0")),
          "generated H2 loader rejects invalid nuclear repulsion");
  require(generatedH2TextRejected(fixture0735 + "unexpected trailing data\n"),
          "generated H2 loader rejects trailing content");
  bool badParameterRejected = false;
  try {
    evaluateGeneratedH2Energy(
        model0735, std::numeric_limits<double>::quiet_NaN());
  } catch (const std::invalid_argument &) {
    badParameterRejected = true;
  }
  require(badParameterRejected,
          "generated H2 evaluator rejects non-finite circuit parameters");
}

void testHeHEnergy() {
  using State16 = std::array<Complex, 16>;
  using Dense16 = std::array<std::array<Complex, 16>, 16>;
  const QubitHamiltonian hamiltonian = hehHamiltonian();
  require(hamiltonian.qubitCount() == 4 &&
              hamiltonian.terms().size() == 27,
          "frozen HeH+ Hamiltonian has four qubits and all 27 Pauli terms");

  State16 state{};
  std::vector<Complex> stateVector;
  double normSquared = 0.0;
  for (std::size_t index = 0; index < state.size(); ++index) {
    state[index] = Complex(
        static_cast<double>(static_cast<int>(index) - 7) / 11.0,
        static_cast<double>(static_cast<int>(index % 5) - 2) / 13.0
    );
    normSquared += std::norm(state[index]);
  }
  for (Complex &amplitude : state) {
    amplitude /= std::sqrt(normSquared);
    stateVector.push_back(amplitude);
  }
  Dense16 dense{};
  for (const PauliTerm &pauliTerm : hamiltonian.terms()) {
    for (std::size_t column = 0; column < state.size(); ++column) {
      std::size_t row = column;
      Complex phase(1.0, 0.0);
      for (std::size_t qubit = 0; qubit < 4; ++qubit) {
        const std::size_t mask = std::size_t{1} << qubit;
        const bool one = (column & mask) != 0;
        switch (pauliTerm.operators[qubit]) {
        case Pauli::I: break;
        case Pauli::X: row ^= mask; break;
        case Pauli::Y:
          row ^= mask;
          phase *= one ? Complex(0.0, -1.0) : Complex(0.0, 1.0);
          break;
        case Pauli::Z: if (one) phase = -phase; break;
        }
      }
      dense[row][column] += pauliTerm.coefficient * phase;
    }
  }
  Complex denseEnergy(0.0, 0.0);
  double termEnergy = 0.0;
  for (std::size_t row = 0; row < state.size(); ++row) {
    for (std::size_t column = 0; column < state.size(); ++column) {
      denseEnergy += std::conj(state[row]) * dense[row][column] * state[column];
    }
  }
  for (const PauliTerm &pauliTerm : hamiltonian.terms()) {
    termEnergy += pauliTerm.coefficient *
                  pauliExpectation(stateVector, pauliTerm.operators);
  }
  require(std::abs(denseEnergy.imag()) < 1e-12,
          "dense HeH+ Hamiltonian expectation is real");
  require(std::abs(termEnergy - denseEnergy.real()) < 1e-12,
          "four-qubit termwise Pauli expectation matches dense matrix");

  const double parameters[][2] = {
      {0.0, 0.0}, {0.6, 1.7}, {1.2, 4.1},
      {1.5707963267948966, 3.4}
  };
  for (const auto &parameter : parameters) {
    Simulator simulator(4);
    simulator.run(hehSingletTrialCircuit(parameter[0], parameter[1]));
    const auto &prepared = simulator.stateVector();
    double norm = 0.0;
    for (const Complex &amplitude : prepared) norm += std::norm(amplitude);
    require(std::abs(norm - 1.0) < 1e-12,
            "HeH+ trial circuit remains normalized");
    for (std::size_t basis = 0; basis < prepared.size(); ++basis) {
      if (basis != 5 && basis != 6 && basis != 9 && basis != 10) {
        require(std::abs(prepared[basis]) < 1e-12,
                "HeH+ trial state remains in the two-electron sector");
      }
    }
    require(std::abs(prepared[6] - prepared[9]) < 1e-12,
            "HeH+ trial state has zero-spin open-shell symmetry");
    const double expected[3] = {
        std::cos(parameter[0]),
        std::sin(parameter[0]) * std::cos(parameter[1]) / std::sqrt(2.0),
        std::sin(parameter[0]) * std::cos(parameter[1]) / std::sqrt(2.0)
    };
    const double expectedD1 =
        std::sin(parameter[0]) * std::sin(parameter[1]);
    require(std::abs(prepared[5].real() - expected[0]) < 1e-11 &&
                std::abs(prepared[6].real() - expected[1]) < 1e-11 &&
                std::abs(prepared[9].real() - expected[2]) < 1e-11 &&
                std::abs(prepared[10].real() - expectedD1) < 1e-11,
            "HeH+ circuit realizes the documented singlet CSF amplitudes");
    require(std::abs(hamiltonian.expectation(simulator) -
                     evaluateHeHEnergy(parameter[0], parameter[1])) < 1e-12,
            "HeH+ parameter evaluations use fresh circuit state");
  }

  const double exact = exactHeHSingletEnergy();
  require(std::abs(exact - HEH_ELECTRONIC_REFERENCE_HARTREE) < 1e-10,
          "HeH+ singlet-sector minimum matches independent Psi4 FCI");
  const HeHResult first = optimizeHeH();
  const HeHResult second = optimizeHeH();
  require(first.theta == second.theta && first.phi == second.phi &&
              first.electronicEnergy == second.electronicEnergy &&
              first.energyEvaluations == second.energyEvaluations,
          "bounded multi-start HeH+ optimizer is deterministic");
  require(first.energyEvaluations <= 10000 &&
              std::abs(first.electronicEnergy - exact) < 1e-8,
          "HeH+ variational energy converges to the exact singlet minimum");
  require(first.electronicEnergy >= exact - 1e-10,
          "HeH+ variational energy does not undercut exact reference");
  require(first.nuclearRepulsionEnergy ==
              HEH_NUCLEAR_REPULSION_HARTREE &&
              std::abs(first.totalEnergy -
                       (first.electronicEnergy + first.nuclearRepulsionEnergy)) <
                  1e-14,
          "HeH+ total adds nuclear repulsion exactly once");
  require(std::abs(first.totalEnergy - HEH_TOTAL_REFERENCE_HARTREE) < 1e-8,
          "HeH+ total energy matches independent Psi4 reference");
}

void testLexerParser() {
  std::vector<Token> tokens = Lexer("qubit q0\nRZ(-2 * (3 + 1)) q0\n").tokenize();
  require(tokens.front().type == TokenType::KeywordQubit, "qubit keyword tokenized");
  require(tokens.back().type == TokenType::EndOfFile, "EOF token emitted");
  Program program = Parser(tokens).parse();
  require(program.statements.size() == 2, "two statements parsed");
  const Expression &expression = program.statements[1].parameter;
  require(expression.type == ExpressionType::BinaryOperation &&
              expression.operation == '*',
          "multiplication is outer precedence level");
  require(expression.children[0].type == ExpressionType::UnaryOperation,
          "unary minus parsed as a factor");
  const Expression &parenthesized = expression.children[1];
  require(parenthesized.operation == '+', "parentheses preserve addition grouping");

  Program precedence = parse("qubit q0\nRZ(2 + 3 * 4) q0\n");
  const Expression &sum = precedence.statements[1].parameter;
  require(sum.operation == '+', "addition is expression root");
  require(sum.children[1].operation == '*', "multiplication binds before addition");

  bool malformedRejected = false;
  try {
    parse("qubit q0\nRZ(2 +) q0\n");
  } catch (const std::exception &) {
    malformedRejected = true;
  }
  require(malformedRejected, "malformed expression rejected");
  bool syntaxRejected = false;
  try {
    parse("qubit q0\nH\n");
  } catch (const std::exception &) {
    syntaxRejected = true;
  }
  require(syntaxRejected, "invalid syntax rejected");
}

void testSemantics() {
  Program program;
  auto result = check("qubit q0\nqubit q0\n", program);
  require(!result.errors.empty(), "duplicate qubit declaration rejected");
  result = check("H missing\n", program);
  require(!result.errors.empty(), "undeclared qubit rejected");
  result = check("RZ(theta) q0\nqubit q0\n", program);
  require(!result.errors.empty(), "qubit and variable use before declaration rejected");
  result = check("qubit q0, q1\nH q0, q1\n", program);
  require(!result.errors.empty(), "invalid gate arity rejected");
  result = check("qubit q0, q1\nCNOT q0\n", program);
  require(!result.errors.empty(), "invalid CNOT arity rejected");
  result = check("qubit q0\nCNOT q0, q0\n", program);
  require(!result.errors.empty(), "identical CNOT operands rejected");
  result = check("qubit q0\nRZ(1 / 0) q0\n", program);
  require(!result.errors.empty(), "division by zero rejected");
  result = check("let theta = 2\nlet theta = 3\n", program);
  require(!result.errors.empty(), "duplicate variable declaration rejected");
  result = check("let later = earlier\n", program);
  require(!result.errors.empty(), "variable use before declaration rejected");
}

Instruction makeInstruction(InstructionType type, const std::string &gate,
                            std::vector<int> qubits, double parameter = 0.0) {
  Instruction instruction;
  instruction.type = type;
  instruction.gateName = gate;
  instruction.qubitIndices = std::move(qubits);
  instruction.parameterValue = parameter;
  return instruction;
}

bool sameInstructions(const std::vector<Instruction> &left,
                      const std::vector<Instruction> &right) {
  if (left.size() != right.size()) return false;
  for (std::size_t index = 0; index < left.size(); ++index) {
    if (left[index].type != right[index].type ||
        left[index].gateName != right[index].gateName ||
        left[index].qubitIndices != right[index].qubitIndices ||
        left[index].parameterValue != right[index].parameterValue) {
      return false;
    }
  }
  return true;
}

void requireSameState(const std::vector<Complex> &left,
                      const std::vector<Complex> &right,
                      const std::string &message) {
  // State-vector equivalence tolerance: absolute complex-amplitude error < 1e-12.
  require(left.size() == right.size(), message + " (dimension)");
  for (std::size_t index = 0; index < left.size(); ++index) {
    require(std::abs(left[index] - right[index]) < 1e-12,
            message + " (amplitude " + std::to_string(index) + ")");
  }
}

void testInstructionOptimization() {
  Program declaredCircuit;
  const SemanticResult declaredSemantics =
      check("qubit q0, q1\n", declaredCircuit);
  require(declaredSemantics.errors.empty(),
          "optimizer test circuit has validated qubit declarations");
  const std::size_t qubitCount = declaredSemantics.qubitNames.size();
  const auto optimize = [qubitCount](
                           const std::vector<Instruction> &instructions) {
    return optimizeInstructions(instructions, qubitCount);
  };

  const Instruction h = makeInstruction(InstructionType::ApplyGate, "H", {0});
  const Instruction x = makeInstruction(InstructionType::ApplyGate, "X", {0});
  const Instruction z = makeInstruction(InstructionType::ApplyGate, "Z", {0});
  const Instruction cnot =
      makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1});
  for (const Instruction &gate : {h, x, z, cnot}) {
    require(optimize({gate, gate}).empty(),
            gate.gateName + " identical adjacent pair cancels");
  }
  const int lastValidQubit = static_cast<int>(qubitCount) - 1;
  const Instruction boundaryH =
      makeInstruction(InstructionType::ApplyGate, "H", {lastValidQubit});
  require(optimize({boundaryH, boundaryH}).empty(),
          "gate at the last declared qubit remains eligible");
  const int firstOutOfRangeQubit = static_cast<int>(qubitCount);
  const Instruction outOfRangeH =
      makeInstruction(InstructionType::ApplyGate, "H", {firstOutOfRangeQubit});
  require(sameInstructions(optimize({outOfRangeH, outOfRangeH}),
                           {outOfRangeH, outOfRangeH}),
          "identical H gates at a positive out-of-range index are preserved");
  const Instruction outOfRangeCnot = makeInstruction(
      InstructionType::ApplyGate, "CNOT", {0, firstOutOfRangeQubit});
  require(sameInstructions(optimize({outOfRangeCnot, outOfRangeCnot}),
                           {outOfRangeCnot, outOfRangeCnot}),
          "CNOT with a positive out-of-range operand is preserved");
  require(sameInstructions(optimize({h, h, h}), {h}),
          "repeated H-H-H cancellation leaves H");
  require(sameInstructions(
              optimize({h, makeInstruction(
                                           InstructionType::ApplyGate, "H", {1})}),
              {h, makeInstruction(InstructionType::ApplyGate, "H", {1})}),
          "single-qubit gates on different qubits do not cancel");
  const Instruction reversedCnot =
      makeInstruction(InstructionType::ApplyGate, "CNOT", {1, 0});
  require(sameInstructions(optimize({cnot, reversedCnot}),
                           {cnot, reversedCnot}),
          "CNOT with reversed operands does not cancel");
  require(sameInstructions(optimize({h, x, h}), {h, x, h}),
          "non-adjacent gates do not cancel across an intervening gate");

  const Instruction rz =
      makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 0.75);
  require(sameInstructions(optimize({rz, rz}), {rz, rz}),
          "RZ instructions are preserved");

  const Instruction measurement =
      makeInstruction(InstructionType::Measure, "measure", {1, 0});
  const Instruction secondMeasurement =
      makeInstruction(InstructionType::Measure, "measure", {0});
  const std::vector<Instruction> aroundMeasurement{
      h, h, measurement, x, x, secondMeasurement};
  const std::vector<Instruction> optimizedMeasurement =
      optimize(aroundMeasurement);
  require(sameInstructions(optimizedMeasurement,
                           {measurement, secondMeasurement}),
          "cancellations preserve measurement barriers, order, and qubit lists");
  require(std::count_if(optimizedMeasurement.begin(), optimizedMeasurement.end(),
                        [](const Instruction &instruction) {
                          return instruction.type == InstructionType::Measure;
                        }) == 2,
          "measurement count is preserved");

  const Instruction unknown =
      makeInstruction(InstructionType::ApplyGate, "UNKNOWN", {0});
  const Instruction malformedH =
      makeInstruction(InstructionType::ApplyGate, "H", {0, 1});
  const Instruction malformedCnot =
      makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 0});
  const Instruction malformedParameter =
      makeInstruction(InstructionType::ApplyGate, "H", {0}, 1.0);
  const Instruction malformedIndex =
      makeInstruction(InstructionType::ApplyGate, "X", {-1});
  for (const Instruction &invalid : {unknown, malformedH, malformedCnot,
                                     malformedParameter, malformedIndex}) {
    require(sameInstructions(optimize({invalid, invalid}),
                             {invalid, invalid}),
            "unknown or malformed instructions are preserved");
  }

  const std::vector<Instruction> fixedPoint{
      h, x, makeInstruction(InstructionType::ApplyGate, "Z", {1}), cnot, rz,
      measurement, h, unknown, malformedH};
  const std::vector<Instruction> optimized =
      optimize(fixedPoint);
  require(sameInstructions(optimize(optimized), optimized),
          "optimization pass is idempotent");
  for (std::size_t index = 1; index < optimized.size(); ++index) {
    require(!optimize(
                 {optimized[index - 1], optimized[index]}).empty(),
            "no eligible adjacent pair remains");
  }

  const std::vector<std::vector<Instruction>> contexts{
      {makeInstruction(InstructionType::ApplyGate, "X", {1}), h, h, cnot, cnot},
      {makeInstruction(InstructionType::ApplyGate, "H", {0}), x, x, rz,
       rz, makeInstruction(InstructionType::ApplyGate, "H", {1}), z, z}};
  for (const std::vector<Instruction> &circuit : contexts) {
    const std::vector<Instruction> reduced = optimize(circuit);
    Simulator originalSimulator(2);
    Simulator optimizedSimulator(2);
    originalSimulator.run(circuit);
    optimizedSimulator.run(reduced);
    requireSameState(originalSimulator.stateVector(),
                     optimizedSimulator.stateVector(),
                     "optimized and original final state vectors agree");
  }
}

void testCodeGeneration() {
  Program program;
  SemanticResult result = check(
      "qubit second, first\nlet theta = 2 + 3 * 4\nRZ(theta) first\n"
      "CNOT second, first\nmeasure first, second\n",
      program);
  require(result.errors.empty(), "valid code-generation input accepted");
  std::vector<Instruction> instructions = generateInstructions(program, result);
  require(instructions.size() == 3, "three executable instructions generated");
  require(instructions[0].type == InstructionType::ApplyGate &&
              instructions[0].gateName == "RZ" &&
              instructions[0].qubitIndices == std::vector<int>{1} &&
              std::abs(instructions[0].parameterValue - 14.0) < 1e-12,
          "RZ instruction has evaluated parameter and declaration-order operand");
  require(instructions[1].gateName == "CNOT" &&
              instructions[1].qubitIndices == std::vector<int>({0, 1}),
          "CNOT operands use declaration-order indices");
  require(instructions[2].type == InstructionType::Measure &&
              instructions[2].qubitIndices == std::vector<int>({1, 0}),
          "measurement instruction and operand order preserved");
}

void testSimulator() {
  std::string initial = run("qubit q0\nmeasure q0\n");
  require(initial.find("0: 1.000000") != std::string::npos,
          "initial state has all probability at zero");
  std::string x = run("qubit q0\nX q0\nmeasure q0\n");
  require(x.find("1: 1.000000") != std::string::npos, "X flips initial basis state");
  std::string h = run("qubit q0\nH q0\nmeasure q0\n");
  require(h.find("0: 0.500000") != std::string::npos &&
              h.find("1: 0.500000") != std::string::npos,
          "H produces equal basis probabilities");
  std::string bell = run("qubit q0, q1\nH q0\nCNOT q0, q1\nmeasure q0, q1\n");
  require(bell.find("00: 0.500000") != std::string::npos &&
              bell.find("11: 0.500000") != std::string::npos &&
              bell.find("01:") == std::string::npos &&
              bell.find("10:") == std::string::npos,
          "CNOT creates expected two-qubit correlated state");
  std::string rz = run("qubit q0\nH q0\nRZ(1.2) q0\nmeasure q0\n");
  require(rz.find("0: 0.500000") != std::string::npos &&
              rz.find("1: 0.500000") != std::string::npos,
          "RZ preserves probabilities and normalization");
  require(run("qubit q0\nH q0\nmeasure q0\n") ==
              run("qubit q0\nH q0\nmeasure q0\n"),
          "fixed-seed measurement output is repeatable");
}

void testCapacity() {
  for (int count : {-1, 25}) {
    bool rejected = false;
    try {
      Simulator simulator(count);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    require(rejected, "negative and over-limit qubit counts rejected");
  }

  auto rejectsInstruction = [](const Instruction &instruction) {
    try {
      Simulator simulator(2);
      simulator.run(std::vector<Instruction>{instruction});
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "X", {-1}, 0.0}),
          "negative single-qubit index rejected");
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "X", {2}, 0.0}),
          "index equal to qubit count rejected");
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "CNOT", {-1, 1}, 0.0}),
          "negative CNOT control index rejected");
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "CNOT", {0, 2}, 0.0}),
          "out-of-range CNOT target index rejected");
  require(rejectsInstruction(
              {InstructionType::Measure, "measure", {0, -1}, 0.0}),
          "invalid measurement index rejected");
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "CNOT", {0}, 0.0}),
          "missing CNOT operand rejected");
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "UNKNOWN", {0}, 0.0}),
          "unknown gate instruction rejected");

  Simulator smallCircuit(2);
  std::vector<Instruction> instructions;
  instructions.push_back({InstructionType::ApplyGate, "X", {1}, 0.0});
  instructions.push_back({InstructionType::Measure, "measure", {0, 1}, 0.0});
  std::ostringstream output;
  std::streambuf *original = std::cout.rdbuf(output.rdbuf());
  smallCircuit.run(instructions);
  std::cout.rdbuf(original);
  require(output.str().find("01: 1.000000") != std::string::npos,
          "supported small circuit executes normally");
}

} // namespace

int main() {
  try {
    testPauliExpectations();
    testHamiltonianExpectations();
    testSymmetricEigensolverValidation();
    testBenchmarkCatalog();
    testHeHEnergy();
    testH2Energy();
    testGeneratedH2();
    testLexerParser();
    testSemantics();
    testInstructionOptimization();
    testCodeGeneration();
    testSimulator();
    testCapacity();
  } catch (const std::exception &error) {
    std::cerr << "FAIL after " << assertions << " assertions: " << error.what() << '\n';
    return 1;
  }
  std::cout << "Passed " << assertions << " assertions.\n";
  return 0;
}
