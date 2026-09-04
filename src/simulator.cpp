#include "simulator.h"

#include <cmath>
#include <iomanip>
#include <iostream>

Simulator::Simulator(int numberOfQubits)
    : qubitCount(numberOfQubits), randomGenerator(12345) {
  int stateSize = 1 << qubitCount;
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
  int mask = 1 << qubit;

  for (int i = 0; i < static_cast<int>(state.size()); i = i + 1) {
    if ((i & mask) == 0) {
      int pairedIndex = i | mask;
      std::complex<double> zeroAmplitude = state[i];
      std::complex<double> oneAmplitude = state[pairedIndex];

      state[i] = matrix00 * zeroAmplitude + matrix01 * oneAmplitude;
      state[pairedIndex] = matrix10 * zeroAmplitude + matrix11 * oneAmplitude;
    }
  }
}

void Simulator::applyCNOT(int control, int target) {
  int controlMask = 1 << control;
  int targetMask = 1 << target;

  for (int i = 0; i < static_cast<int>(state.size()); i = i + 1) {
    if ((i & controlMask) != 0 && (i & targetMask) == 0) {
      int pairedIndex = i | targetMask;
      std::complex<double> temporary = state[i];
      state[i] = state[pairedIndex];
      state[pairedIndex] = temporary;
    }
  }
}

std::string Simulator::bitString(int stateIndex) const {
  std::string result;

  for (int i = 0; i < qubitCount; i = i + 1) {
    int mask = 1 << i;
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
  for (int i = 0; i < static_cast<int>(state.size()); i = i + 1) {
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
  int outcome = distribution(randomGenerator);
  std::string fullResult = bitString(outcome);
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
