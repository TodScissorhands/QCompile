#include "chemistry/h2_energy.h"

#include "chemistry/symmetric_eigensolver.h"
#include "compiler/code_generator.h"
#include "simulator/simulator.h"

#include <algorithm>
#include <cstddef>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kSearchTolerance = 1e-12;
constexpr int kMaximumIterations = 128;
constexpr double kReferenceTolerance = 1e-8;

std::vector<Instruction> h2TrialCircuit(double parameter) {
  return {
      {InstructionType::ApplyGate, "X", {0}, 0.0},
      {InstructionType::ApplyGate, "CNOT", {0, 1}, 0.0},
      {InstructionType::ApplyGate, "RZ", {0}, -kPi / 2.0},
      {InstructionType::ApplyGate, "H", {0}, 0.0},
      {InstructionType::ApplyGate, "RZ", {0}, parameter},
      {InstructionType::ApplyGate, "H", {0}, 0.0},
      {InstructionType::ApplyGate, "RZ", {0}, kPi / 2.0},
      {InstructionType::ApplyGate, "CNOT", {0, 1}, 0.0}
  };
}

std::array<std::array<double, 4>, 4> buildMatrix() {
  const QubitHamiltonian hamiltonian = h2Hamiltonian();
  std::array<std::array<double, 4>, 4> matrix{};
  for (const PauliTerm &term : hamiltonian.terms()) {
    for (std::size_t column = 0; column < 4; ++column) {
      std::size_t row = column;
      double phase = 1.0;
      for (std::size_t qubit = 0; qubit < term.operators.size(); ++qubit) {
        const std::size_t mask = std::size_t{1} << qubit;
        const bool bitIsOne = (column & mask) != 0;
        switch (term.operators[qubit]) {
        case Pauli::I:
          break;
        case Pauli::X:
          row ^= mask;
          break;
        case Pauli::Z:
          if (bitIsOne) {
            phase = -phase;
          }
          break;
        case Pauli::Y:
          throw std::logic_error("H2 reference matrix does not use Y terms");
        }
      }
      matrix[row][column] += term.coefficient * phase;
    }
  }
  return matrix;
}


} // namespace

std::vector<Instruction> h2TrialInstructions(double parameter) {
  if (!std::isfinite(parameter) || parameter < 0.0 || parameter > kPi / 2.0) {
    throw std::invalid_argument("H2 trial parameter must be in [0, pi/2]");
  }
  return h2TrialCircuit(parameter);
}

QubitHamiltonian h2Hamiltonian() {
  return QubitHamiltonian(2, {
      {-1.05237325, {Pauli::I, Pauli::I}},
      {0.39793742, {Pauli::Z, Pauli::I}},
      {-0.39793742, {Pauli::I, Pauli::Z}},
      {-0.01128010, {Pauli::Z, Pauli::Z}},
      {0.18093120, {Pauli::X, Pauli::X}}
  });
}

double evaluateH2Energy(double parameter) {
  if (!std::isfinite(parameter) || parameter < 0.0 || parameter > kPi / 2.0) {
    throw std::invalid_argument("H2 trial parameter must be in [0, pi/2]");
  }
  Simulator simulator(2);
  simulator.run(h2TrialInstructions(parameter));
  return h2Hamiltonian().expectation(simulator);
}

double exactH2GroundEnergy() {
  return lowestSymmetricEigenvalue(buildMatrix());
}

H2Result optimizeH2(VqeCallback callback) {
  const double exactEnergy = exactH2GroundEnergy();
  if (std::abs(exactEnergy - H2_ELECTRONIC_REFERENCE_HARTREE) >
      kReferenceTolerance) {
    throw std::runtime_error("H2 exact energy does not match the documented reference");
  }

  int evaluations = 0;
  double bestEnergy = std::numeric_limits<double>::infinity();

  auto evalWithCallback = [&](double p) {
    const double elEnergy = evaluateH2Energy(p);
    ++evaluations;
    const bool isBest = (elEnergy < bestEnergy);
    if (isBest) {
      bestEnergy = elEnergy;
    }
    if (callback) {
      callback({
          evaluations,
          p,
          0.0,
          elEnergy,
          elEnergy + H2_NUCLEAR_REPULSION_HARTREE,
          isBest
      });
    }
    return elEnergy;
  };

  constexpr double inverseGoldenRatio = 0.6180339887498948482;
  double lower = 0.0;
  double upper = kPi / 2.0;
  double left = upper - inverseGoldenRatio * (upper - lower);
  double right = lower + inverseGoldenRatio * (upper - lower);
  double leftEnergy = evalWithCallback(left);
  double rightEnergy = evalWithCallback(right);
  int iterations = 0;

  while (upper - lower > kSearchTolerance &&
         iterations < kMaximumIterations) {
    if (leftEnergy < rightEnergy) {
      upper = right;
      right = left;
      rightEnergy = leftEnergy;
      left = upper - inverseGoldenRatio * (upper - lower);
      leftEnergy = evalWithCallback(left);
    } else {
      lower = left;
      left = right;
      leftEnergy = rightEnergy;
      right = lower + inverseGoldenRatio * (upper - lower);
      rightEnergy = evalWithCallback(right);
    }
    ++iterations;
  }
  if (upper - lower > kSearchTolerance) {
    throw std::runtime_error("H2 optimizer did not converge within iteration limit");
  }

  const double parameter = 0.5 * (lower + upper);
  const double electronicEnergy = evalWithCallback(parameter);
  if (std::abs(electronicEnergy - exactEnergy) > 1e-10) {
    throw std::runtime_error("Optimized H2 state failed exact-energy validation");
  }

  return {
      parameter,
      electronicEnergy,
      exactEnergy,
      H2_NUCLEAR_REPULSION_HARTREE,
      electronicEnergy + H2_NUCLEAR_REPULSION_HARTREE,
      evaluations
  };
}

namespace {

using Complex = std::complex<double>;
using ComplexMatrix4 = std::array<std::array<Complex, 4>, 4>;

ComplexMatrix4 generatedH2Matrix(const QubitHamiltonian &hamiltonian) {
  if (hamiltonian.qubitCount() != 2) {
    throw std::invalid_argument("Generated H2 Hamiltonian must have two qubits");
  }
  ComplexMatrix4 matrix{};
  for (const PauliTerm &term : hamiltonian.terms()) {
    for (std::size_t column = 0; column < 4; ++column) {
      std::size_t row = column;
      Complex phase(1.0, 0.0);
      for (std::size_t qubit = 0; qubit < term.operators.size(); ++qubit) {
        const std::size_t mask = std::size_t{1} << qubit;
        const bool bitIsOne = (column & mask) != 0;
        switch (term.operators[qubit]) {
        case Pauli::I:
          break;
        case Pauli::X:
          row ^= mask;
          break;
        case Pauli::Y:
          row ^= mask;
          phase *= bitIsOne ? Complex(0.0, -1.0) : Complex(0.0, 1.0);
          break;
        case Pauli::Z:
          if (bitIsOne) {
            phase = -phase;
          }
          break;
        }
      }
      matrix[row][column] += term.coefficient * phase;
    }
  }
  return matrix;
}

std::array<std::array<double, 2>, 2> generatedH2PhysicalSector(
    const QubitHamiltonian &hamiltonian) {
  const ComplexMatrix4 matrix = generatedH2Matrix(hamiltonian);
  // Particle-number-reduced parity encoding leaves the H2 Nalpha=Nbeta=1
  // physical block at basis indices 1 and 2 (q0 is the least-significant bit).
  constexpr std::array<std::size_t, 2> sector{1, 2};
  constexpr double matrixTolerance = 1e-11;
  std::array<std::array<double, 2>, 2> projected{};
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t column = 0; column < 4; ++column) {
      const bool rowInSector = row == sector[0] || row == sector[1];
      const bool columnInSector = column == sector[0] || column == sector[1];
      if (rowInSector != columnInSector &&
          std::abs(matrix[row][column]) > matrixTolerance) {
        throw std::invalid_argument(
            "Generated H2 Hamiltonian couples the physical sector to its complement");
      }
    }
  }
  for (std::size_t row = 0; row < 2; ++row) {
    for (std::size_t column = 0; column < 2; ++column) {
      const Complex value = matrix[sector[row]][sector[column]];
      if (std::abs(value.imag()) > matrixTolerance) {
        throw std::invalid_argument(
            "Generated H2 physical-sector matrix has an imaginary component");
      }
      projected[row][column] = value.real();
    }
  }
  if (std::abs(projected[0][1] - projected[1][0]) > matrixTolerance) {
    throw std::invalid_argument("Generated H2 physical-sector matrix is not symmetric");
  }
  return projected;
}

void validateGeneratedH2TrialState(const Simulator &simulator) {
  const std::vector<Complex> &state = simulator.stateVector();
  if (state.size() != 4) {
    throw std::runtime_error("Generated H2 trial state has an invalid dimension");
  }
  double norm = 0.0;
  for (std::size_t index = 0; index < state.size(); ++index) {
    norm += std::norm(state[index]);
    if (index != 1 && index != 2 && std::abs(state[index]) > 1e-12) {
      throw std::runtime_error(
          "Generated H2 trial circuit left the physical parity sector");
    }
  }
  if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1e-12) {
    throw std::runtime_error("Generated H2 trial state is not normalized");
  }
}

} // namespace

double evaluateGeneratedH2Energy(const GeneratedH2 &model, double parameter) {
  if (!std::isfinite(parameter) || parameter < 0.0 || parameter > kPi / 2.0) {
    throw std::invalid_argument("Generated H2 trial parameter must be in [0, pi/2]");
  }
  Simulator simulator(2);
  simulator.run(h2TrialInstructions(parameter));
  validateGeneratedH2TrialState(simulator);
  return model.electronicHamiltonian.expectation(simulator);
}

double exactGeneratedH2GroundEnergy(const GeneratedH2 &model) {
  const auto matrix = generatedH2PhysicalSector(model.electronicHamiltonian);
  return lowestSymmetricEigenvalue(matrix);
}

GeneratedH2EnergyResult optimizeGeneratedH2(const GeneratedH2 &model) {
  constexpr double inverseGoldenRatio = 0.6180339887498948482;
  constexpr double parameterTolerance = 1e-12;
  constexpr int maximumIterations = 128;
  constexpr double referenceTolerance = 1e-8;
  constexpr double variationalTolerance = 1e-10;

  const double exactEnergy = exactGeneratedH2GroundEnergy(model);
  if (std::abs(exactEnergy - model.referenceElectronicEnergy) >
      referenceTolerance) {
    throw std::runtime_error(
        "Generated H2 physical-sector minimum disagrees with file reference");
  }

  double lower = 0.0;
  double upper = kPi / 2.0;
  double left = upper - inverseGoldenRatio * (upper - lower);
  double right = lower + inverseGoldenRatio * (upper - lower);
  double leftEnergy = evaluateGeneratedH2Energy(model, left);
  double rightEnergy = evaluateGeneratedH2Energy(model, right);
  int evaluations = 2;
  int iterations = 0;
  while (upper - lower > parameterTolerance &&
         iterations < maximumIterations) {
    if (leftEnergy < rightEnergy) {
      upper = right;
      right = left;
      rightEnergy = leftEnergy;
      left = upper - inverseGoldenRatio * (upper - lower);
      leftEnergy = evaluateGeneratedH2Energy(model, left);
      ++evaluations;
    } else {
      lower = left;
      left = right;
      leftEnergy = rightEnergy;
      right = lower + inverseGoldenRatio * (upper - lower);
      rightEnergy = evaluateGeneratedH2Energy(model, right);
      ++evaluations;
    }
    ++iterations;
  }
  if (upper - lower > parameterTolerance) {
    throw std::runtime_error(
        "Generated H2 optimizer did not converge within iteration limit");
  }
  const double parameter = 0.5 * (lower + upper);
  const double electronicEnergy = evaluateGeneratedH2Energy(model, parameter);
  ++evaluations;
  if (electronicEnergy < exactEnergy - variationalTolerance ||
      std::abs(electronicEnergy - exactEnergy) > referenceTolerance) {
    throw std::runtime_error(
        "Generated H2 trial state did not reach the exact physical-sector minimum");
  }
  const double totalEnergy = electronicEnergy + model.nuclearRepulsionEnergy;
  const double exactTotalEnergy = exactEnergy + model.nuclearRepulsionEnergy;
  return {
      parameter,
      electronicEnergy,
      exactEnergy,
      model.referenceElectronicEnergy,
      model.nuclearRepulsionEnergy,
      totalEnergy,
      exactTotalEnergy,
      electronicEnergy - exactEnergy,
      evaluations,
      true
  };
}
