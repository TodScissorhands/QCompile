#pragma once

#include "code_generator.h"

#include <complex>
#include <random>
#include <vector>

class Simulator {
public:
  Simulator(int qubitCount);

  void run(const std::vector<Instruction> &instructions);

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
  std::string bitString(int stateIndex) const;
};
