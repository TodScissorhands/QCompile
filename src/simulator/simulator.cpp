#include "simulator/simulator.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <stdexcept>

#include "chemistry/hamiltonian.h"

Simulator::Simulator(int numberOfQubits)
    : qubitCount(numberOfQubits), randomGenerator(12345) {
  if (qubitCount < 0 || qubitCount > MAX_SIMULATOR_QUBITS) {
    throw std::invalid_argument(
        "Simulator supports between 0 and " +
        std::to_string(MAX_SIMULATOR_QUBITS) + " qubits"
    );
  }
  const std::size_t stateSize =
      std::size_t{1} << static_cast<unsigned int>(qubitCount);
  state.resize(stateSize, std::complex<double>(0.0, 0.0));
  state[0] = std::complex<double>(1.0, 0.0);
}

void Simulator::applySingleQubitGate(
    int qubit,
    std::complex<double> matrix00,
    std::complex<double> matrix01,
    std::complex<double> matrix10,
    std::complex<double> matrix11
) {
  const std::size_t mask = std::size_t{1} << static_cast<unsigned int>(qubit);

  for (std::size_t i = 0; i < state.size(); i = i + 1) {
    if ((i & mask) == 0) {
      const std::size_t pairedIndex = i | mask;
      std::complex<double> zeroAmplitude = state[i];
      std::complex<double> oneAmplitude = state[pairedIndex];

      state[i] = matrix00 * zeroAmplitude + matrix01 * oneAmplitude;
      state[pairedIndex] = matrix10 * zeroAmplitude + matrix11 * oneAmplitude;
    }
  }
}

void Simulator::applyCNOT(int control, int target) {
  const std::size_t controlMask =
      std::size_t{1} << static_cast<unsigned int>(control);
  const std::size_t targetMask =
      std::size_t{1} << static_cast<unsigned int>(target);

  for (std::size_t i = 0; i < state.size(); i = i + 1) {
    if ((i & controlMask) != 0 && (i & targetMask) == 0) {
      const std::size_t pairedIndex = i | targetMask;
      std::complex<double> temporary = state[i];
      state[i] = state[pairedIndex];
      state[pairedIndex] = temporary;
    }
  }
}

std::string Simulator::bitString(std::size_t stateIndex) const {
  std::string result;

  for (int i = 0; i < qubitCount; i = i + 1) {
    const std::size_t mask = std::size_t{1} << static_cast<unsigned int>(i);
    if ((stateIndex & mask) == 0) {
      result = result + "0";
    } else {
      result = result + "1";
    }
  }

  return result;
}

void Simulator::measure(const std::vector<int> &measuredQubits) {
  std::vector<double> probabilities;

  std::cout << "\nProbabilities:\n";
  for (std::size_t i = 0; i < state.size(); i = i + 1) {
    double probability = std::norm(state[i]);
    probabilities.push_back(probability);

    if (probability > 0.000001) {
      std::cout << "  " << bitString(i) << ": " << std::fixed
                << std::setprecision(6) << probability << "\n";
    }
  }

  std::discrete_distribution<int> distribution(
      probabilities.begin(),
      probabilities.end()
  );
  const int outcome = distribution(randomGenerator);
  std::string fullResult = bitString(static_cast<std::size_t>(outcome));
  std::string measuredResult;

  for (int i = 0; i < static_cast<int>(measuredQubits.size()); i = i + 1) {
    measuredResult = measuredResult + fullResult[measuredQubits[i]];
  }

  std::cout << "Result: " << measuredResult << "\n";
}

void Simulator::run(const std::vector<Instruction> &instructions) {
  double inverseSquareRootOfTwo = 1.0 / std::sqrt(2.0);

  for (int i = 0; i < static_cast<int>(instructions.size()); i = i + 1) {
    const Instruction &instruction = instructions[i];

    if (instruction.type != InstructionType::ApplyGate &&
        instruction.type != InstructionType::Measure) {
      throw std::invalid_argument("Unsupported instruction type");
    }

    if (instruction.type == InstructionType::ApplyGate) {
      const bool singleQubitGate =
          instruction.gateName == "H" || instruction.gateName == "X" ||
          instruction.gateName == "Z" || instruction.gateName == "RZ";
      const bool cnotGate = instruction.gateName == "CNOT";
      if (!singleQubitGate && !cnotGate) {
        throw std::invalid_argument("Unsupported gate instruction");
      }
      if ((singleQubitGate && instruction.qubitIndices.size() != 1) ||
          (cnotGate && instruction.qubitIndices.size() != 2)) {
        throw std::invalid_argument("Instruction has an invalid qubit count");
      }
    }

    for (std::size_t j = 0; j < instruction.qubitIndices.size(); j = j + 1) {
      const int qubit = instruction.qubitIndices[j];
      if (qubit < 0 || qubit >= qubitCount) {
        throw std::invalid_argument("Instruction qubit index out of range");
      }
    }

    if (instruction.type == InstructionType::Measure) {
      measure(instruction.qubitIndices);
    } else if (instruction.gateName == "H") {
      applySingleQubitGate(
          instruction.qubitIndices[0],
          inverseSquareRootOfTwo,
          inverseSquareRootOfTwo,
          inverseSquareRootOfTwo,
          -inverseSquareRootOfTwo
      );
    } else if (instruction.gateName == "X") {
      applySingleQubitGate(
          instruction.qubitIndices[0],
          0.0,
          1.0,
          1.0,
          0.0
      );
    } else if (instruction.gateName == "Z") {
      applySingleQubitGate(
          instruction.qubitIndices[0],
          1.0,
          0.0,
          0.0,
          -1.0
      );
    } else if (instruction.gateName == "RZ") {
      double halfAngle = instruction.parameterValue / 2.0;
      std::complex<double> zeroPhase(
          std::cos(-halfAngle),
          std::sin(-halfAngle)
      );
      std::complex<double> onePhase(std::cos(halfAngle), std::sin(halfAngle));
      applySingleQubitGate(
          instruction.qubitIndices[0],
          zeroPhase,
          0.0,
          0.0,
          onePhase
      );
    } else if (instruction.gateName == "CNOT") {
      applyCNOT(instruction.qubitIndices[0], instruction.qubitIndices[1]);
    }
  }
}

const std::vector<std::complex<double>> &Simulator::stateVector() const {
  return state;
}

double Simulator::expectation(const std::vector<Pauli> &operators) const {
  if (operators.size() != static_cast<std::size_t>(qubitCount)) {
    throw std::invalid_argument(
        "Pauli operator qubit count does not match simulator"
    );
  }
  return pauliExpectation(state, operators);
}
