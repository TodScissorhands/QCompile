#pragma once

#include "compiler/code_generator.h"

#include <complex>
#include <cstddef>
#include <optional>
#include <random>
#include <string>
#include <vector>
enum class Pauli;

constexpr int MAX_SIMULATOR_QUBITS = 24;

struct StepResult {
  std::optional<std::size_t> measurementOutcome;
  std::string measurementResult;
};

class Simulator {
public:
  Simulator(int qubitCount);

  void run(const std::vector<Instruction> &instructions);
  StepResult applyInstruction(
      const Instruction &instruction,
      std::optional<std::size_t> fixedOutcome = std::nullopt,
      bool printOutput = false
  );
  void setStateVector(const std::vector<std::complex<double>> &newState);
  int getQubitCount() const;
  double expectation(const std::vector<Pauli> &operators) const;
  const std::vector<std::complex<double>> &stateVector() const;
  std::string bitString(std::size_t stateIndex) const;
private:
  int qubitCount;
  std::vector<std::complex<double>> state;
  std::mt19937 randomGenerator;

  void applySingleQubitGate(
      int qubit,
      std::complex<double> matrix00,
      std::complex<double> matrix01,
      std::complex<double> matrix10,
      std::complex<double> matrix11
  );
  void applyCNOT(int control, int target);
  StepResult measure(
      const std::vector<int> &measuredQubits,
      std::optional<std::size_t> fixedOutcome = std::nullopt,
      bool printOutput = true
  );
};
