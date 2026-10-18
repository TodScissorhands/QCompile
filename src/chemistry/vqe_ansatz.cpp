#include "chemistry/vqe_ansatz.h"
#include "simulator/simulator.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace compiler {

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void appendCnot(std::vector<Instruction> &circuit, int control, int target) {
  circuit.push_back({InstructionType::ApplyGate, "CNOT", {control, target}, 0.0});
}

void appendGate(
    std::vector<Instruction> &circuit,
    const std::string &gate,
    int qubit,
    double parameter = 0.0
) {
  circuit.push_back({InstructionType::ApplyGate, gate, {qubit}, parameter});
}

void appendMultiControlledRY(
    std::vector<Instruction> &circuit,
    int target,
    const std::vector<int> &controls,
    unsigned int activePattern,
    double angle
) {
  const std::size_t count = std::size_t{1} << controls.size();
  std::vector<unsigned int> gray(count);
  for (std::size_t index = 0; index < count; ++index) {
    const unsigned int binary = static_cast<unsigned int>(index);
    gray[index] = binary ^ (binary >> 1U);
  }
  std::vector<double> patternAngles(count, 0.0);
  patternAngles[activePattern] = angle;
  for (std::size_t index = 0; index < count; ++index) {
    double rotation = 0.0;
    for (std::size_t pattern = 0; pattern < count; ++pattern) {
      unsigned int parityValue = gray[index] & static_cast<unsigned int>(pattern);
      unsigned int parity = 0;
      while (parityValue != 0U) {
        parity ^= parityValue & 1U;
        parityValue >>= 1U;
      }
      rotation += (parity == 0U) ? patternAngles[pattern] : -patternAngles[pattern];
    }
    appendDecomposedRY(circuit, target, rotation / static_cast<double>(count));
    const unsigned int next = gray[(index + 1) % count];
    const unsigned int changed = gray[index] ^ next;
    unsigned int changedBits = changed;
    unsigned int controlBit = 0;
    while ((changedBits & 1U) == 0U) {
      ++controlBit;
      changedBits >>= 1U;
    }
    appendCnot(circuit, controls[controlBit], target);
  }
}

void appendGivens(
    std::vector<Instruction> &circuit,
    std::size_t first,
    std::size_t second,
    double angle
) {
  int target = -1;
  std::vector<int> changedBits;
  for (int qubit = 0; qubit < 4; ++qubit) {
    if (((first ^ second) & (std::size_t{1} << qubit)) != 0) {
      if (target < 0) target = qubit;
      changedBits.push_back(qubit);
    }
  }
  if (target < 0) {
    throw std::logic_error("Givens rotation requires distinct basis states");
  }

  std::vector<int> flippedBits;
  for (int qubit = 0; qubit < 4; ++qubit) {
    if ((first & (std::size_t{1} << qubit)) != 0) {
      appendGate(circuit, "X", qubit);
      flippedBits.push_back(qubit);
    }
  }
  std::vector<int> cxControls;
  for (int qubit : changedBits) {
    if (qubit != target) {
      appendCnot(circuit, target, qubit);
      cxControls.push_back(qubit);
    }
  }
  std::vector<int> controls;
  for (int qubit = 0; qubit < 4; ++qubit) {
    if (qubit != target) controls.push_back(qubit);
  }
  appendMultiControlledRY(circuit, target, controls, 0U, 2.0 * angle);
  for (auto iter = cxControls.rbegin(); iter != cxControls.rend(); ++iter) {
    appendCnot(circuit, target, *iter);
  }
  for (auto iter = flippedBits.rbegin(); iter != flippedBits.rend(); ++iter) {
    appendGate(circuit, "X", *iter);
  }
}

void appendSingletRotation(
    std::vector<Instruction> &circuit,
    std::size_t first,
    std::size_t second,
    double angle
) {
  // B maps |9> to (|9>+|6>)/sqrt(2); execute B†, G, B.
  appendGivens(circuit, 9, 6, -kPi / 4.0);
  appendGivens(circuit, first, second, angle);
  appendGivens(circuit, 9, 6, kPi / 4.0);
}

} // namespace

void appendDecomposedRY(std::vector<Instruction> &circuit, int target, double angle) {
  appendGate(circuit, "RZ", target, -kPi / 2.0);
  appendGate(circuit, "H", target);
  appendGate(circuit, "RZ", target, angle);
  appendGate(circuit, "H", target);
  appendGate(circuit, "RZ", target, kPi / 2.0);
}

std::vector<Instruction> createHardwareEfficientAnsatz(
    int qubitCount,
    int depth,
    const std::vector<double> &parameters
) {
  if (qubitCount <= 0) {
    throw std::invalid_argument("Hardware-efficient ansatz requires positive qubit count");
  }
  if (depth < 0) {
    throw std::invalid_argument("Hardware-efficient ansatz requires non-negative depth");
  }

  const std::size_t expectedParams = static_cast<std::size_t>(qubitCount * (depth + 1));
  if (parameters.size() != expectedParams) {
    throw std::invalid_argument("Hardware-efficient ansatz parameter count mismatch (expected " +
                                std::to_string(expectedParams) + ", got " +
                                std::to_string(parameters.size()) + ")");
  }

  std::vector<Instruction> circuit;
  std::size_t paramIdx = 0;

  for (int d = 0; d <= depth; ++d) {
    // Single-qubit RY layer
    for (int q = 0; q < qubitCount; ++q) {
      appendDecomposedRY(circuit, q, parameters[paramIdx++]);
    }

    // Linear entangling CNOT chain between layers
    if (d < depth && qubitCount > 1) {
      for (int q = 0; q + 1 < qubitCount; ++q) {
        appendCnot(circuit, q, q + 1);
      }
    }
  }

  return circuit;
}

std::vector<Instruction> createLiHActiveSpaceAnsatz(double theta, double phi) {
  std::vector<Instruction> circuit;
  // Hartree-Fock reference state |0101> (Li 2s^2 alpha & beta in active orbital 0)
  circuit.push_back({InstructionType::ApplyGate, "X", {0}, 0.0});
  circuit.push_back({InstructionType::ApplyGate, "X", {2}, 0.0});

  // Givens singlet rotations:
  // theta: double excitation between |0101> (5) and |1010> (10)
  appendGivens(circuit, 5, 10, theta);
  // phi: single excitations coupled to singlet CSF 1 (|0110> + |1001>)/sqrt(2)
  appendSingletRotation(circuit, 5, 9, phi);

  return circuit;
}

std::vector<Instruction> createLiHParityAnsatz(double theta) {
  std::vector<Instruction> circuit;
  // Parity-reduced (1,1) electron sector on 2 qubits:
  // Starts in |00> (representing |0101>), applies RY(theta) on q0, then CNOT(q0, q1)
  appendDecomposedRY(circuit, 0, theta);
  appendCnot(circuit, 0, 1);
  return circuit;
}

double evaluateAnsatzEnergy(
    const std::vector<Instruction> &circuit,
    int qubitCount,
    const QubitHamiltonian &hamiltonian
) {
  Simulator sim(qubitCount);
  sim.run(circuit);
  return hamiltonian.expectation(sim);
}

} // namespace compiler
