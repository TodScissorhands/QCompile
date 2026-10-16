#include "chemistry/h2_geometry.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kAngstromToBohr = 1.0 / 0.52917721092;
constexpr double kMetadataTolerance = 1e-12;
constexpr double kNuclearTolerance = 2e-12;
constexpr std::array<const char *, 5> kExpectedLabels{
    "II", "IZ", "XX", "ZI", "ZZ"
};

std::string nextLine(std::istream &input, std::size_t lineNumber) {
  std::string line;
  if (!std::getline(input, line)) {
    throw std::invalid_argument("Generated H2 file ended before line " +
                                std::to_string(lineNumber));
  }
  if (line.empty()) {
    throw std::invalid_argument("Generated H2 file has an empty line at " +
                                std::to_string(lineNumber));
  }
  return line;
}

template <typename T>
T parseValueLine(std::istream &input, const std::string &key,
                 std::size_t lineNumber) {
  const std::string line = nextLine(input, lineNumber);
  std::istringstream stream(line);
  std::string parsedKey;
  T value{};
  std::string extra;
  if (!(stream >> parsedKey >> value) || parsedKey != key || (stream >> extra)) {
    throw std::invalid_argument("Invalid generated H2 field at line " +
                                std::to_string(lineNumber) + ": expected " + key);
  }
  return value;
}

void expectLine(std::istream &input, const std::string &expected,
                std::size_t lineNumber) {
  if (nextLine(input, lineNumber) != expected) {
    throw std::invalid_argument("Unsupported generated H2 metadata at line " +
                                std::to_string(lineNumber));
  }
}

double parseFiniteDoubleLine(std::istream &input, const std::string &key,
                             std::size_t lineNumber) {
  const double value = parseValueLine<double>(input, key, lineNumber);
  if (!std::isfinite(value)) {
    throw std::invalid_argument("Generated H2 field must be finite: " + key);
  }
  return value;
}

std::vector<Pauli> qiskitLabelToQCompile(const std::string &label) {
  if (label.size() != 2) {
    throw std::invalid_argument("Generated H2 Pauli labels must contain two symbols");
  }
  std::vector<Pauli> operators(2);
  for (std::size_t q = 0; q < 2; ++q) {
    switch (label[1 - q]) {
    case 'I': operators[q] = Pauli::I; break;
    case 'X': operators[q] = Pauli::X; break;
    case 'Y': operators[q] = Pauli::Y; break;
    case 'Z': operators[q] = Pauli::Z; break;
    default:
      throw std::invalid_argument("Generated H2 Pauli label contains an unknown symbol");
    }
  }
  return operators;
}

} // namespace

GeneratedH2::GeneratedH2(double distanceAngstromValue,
                         double distanceBohrValue,
                         double nuclearRepulsionEnergyValue,
                         double referenceElectronicEnergyValue,
                         QubitHamiltonian hamiltonian)
    : distanceAngstrom(distanceAngstromValue),
      distanceBohr(distanceBohrValue),
      nuclearRepulsionEnergy(nuclearRepulsionEnergyValue),
      referenceElectronicEnergy(referenceElectronicEnergyValue),
      electronicHamiltonian(std::move(hamiltonian)) {}

GeneratedH2 loadGeneratedH2(const std::string &path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Cannot open generated H2 file: " + path);
  }
  expectLine(input, "QCOMPILE-H2-PAULI 1", 1);
  const double distanceAngstrom = parseFiniteDoubleLine(input, "distance_angstrom", 2);
  const double distanceBohr = parseFiniteDoubleLine(input, "distance_bohr", 3);
  if (distanceAngstrom <= 0.0 || distanceBohr <= 0.0) {
    throw std::invalid_argument("Generated H2 bond distance must be positive");
  }
  if (std::abs(distanceBohr - distanceAngstrom * kAngstromToBohr) >
      kMetadataTolerance * std::max(1.0, distanceBohr)) {
    throw std::invalid_argument("Generated H2 Angstrom/bohr distances disagree");
  }

  expectLine(input, "charge 0", 4);
  expectLine(input, "electrons 2", 5);
  expectLine(input, "spin 0", 6);
  expectLine(input, "multiplicity 1", 7);
  expectLine(input, "basis STO-3G", 8);
  expectLine(input, "energy_unit Hartree", 9);
  expectLine(input, "spatial_orbitals 2", 10);
  expectLine(input, "spin_orbitals 4", 11);
  expectLine(input, "qubits 2", 12);
  expectLine(input, "mapping PARITY", 13);
  expectLine(input, "particle_counts_alpha_beta 1 1", 14);
  expectLine(input, "spin_orbital_order alpha-block-then-beta-block", 15);
  expectLine(input, "parity_reduction fixed_p1=1_fixed_p3=0", 16);
  expectLine(input, "pauli_label_order q1q0_rightmost_is_q0", 17);
  expectLine(input, "qcompile_operator_order q0_to_q1_q0_is_lsb", 18);
  expectLine(input, "electronic_energy_convention nuclear_repulsion_excluded", 19);
  const double nuclear = parseFiniteDoubleLine(input, "nuclear_repulsion_hartree", 20);
  const double reference = parseFiniteDoubleLine(
      input, "singlet_fci_electronic_reference_hartree", 21);
  if (nuclear <= 0.0 || reference >= 0.0) {
    throw std::invalid_argument("Generated H2 nuclear/reference energy is invalid");
  }
  const double expectedNuclear = 1.0 / distanceBohr;
  if (!std::isfinite(expectedNuclear) || expectedNuclear <= 0.0 ||
      std::abs(nuclear - expectedNuclear) > kNuclearTolerance) {
    throw std::invalid_argument("Generated H2 nuclear repulsion is inconsistent with geometry");
  }
  const unsigned int termCount = parseValueLine<unsigned int>(input, "term_count", 22);
  if (termCount != kExpectedLabels.size()) {
    throw std::invalid_argument("Generated H2 term count is unsupported");
  }
  expectLine(input, "TERMS", 23);

  std::vector<PauliTerm> terms;
  terms.reserve(termCount);
  std::set<std::string> labels;
  for (unsigned int index = 0; index < termCount; ++index) {
    const std::size_t lineNumber = 24 + index;
    const std::string line = nextLine(input, lineNumber);
    std::istringstream stream(line);
    std::string label;
    double coefficient = 0.0;
    std::string extra;
    if (!(stream >> label >> coefficient) || (stream >> extra)) {
      throw std::invalid_argument("Malformed generated H2 Pauli term at line " +
                                  std::to_string(lineNumber));
    }
    if (label.size() != 2 || !std::isfinite(coefficient)) {
      throw std::invalid_argument("Invalid generated H2 Pauli label or coefficient");
    }
    if (!labels.insert(label).second) {
      throw std::invalid_argument("Duplicate generated H2 Pauli label: " + label);
    }
    terms.push_back({coefficient, qiskitLabelToQCompile(label)});
  }
  for (const char *expectedLabel : kExpectedLabels) {
    if (labels.find(expectedLabel) == labels.end()) {
      throw std::invalid_argument("Missing required generated H2 Pauli term: " +
                                  std::string(expectedLabel));
    }
  }
  std::string trailing;
  while (std::getline(input, trailing)) {
    throw std::invalid_argument("Unexpected trailing content in generated H2 file");
  }
  if (!input.eof()) {
    throw std::runtime_error("Error reading generated H2 file: " + path);
  }

  return GeneratedH2(distanceAngstrom, distanceBohr, nuclear, reference,
                     QubitHamiltonian(2, std::move(terms)));
}
