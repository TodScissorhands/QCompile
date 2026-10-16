#pragma once

#include "compiler/code_generator.h"

#include <complex>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

enum class Pauli;

constexpr int MAX_SIMULATOR_QUBITS = 24;

class Simulator {
public:
  Simulator(int qubitCount);

  void run(const std::vector<Instruction> &instructions);
  double expectation(const std::vector<Pauli> &operators) const;
  const std::vector<std::complex<double>> &stateVector() const;

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
  void measure(const std::vector<int> &measuredQubits);
  std::string bitString(std::size_t stateIndex) const;
};
