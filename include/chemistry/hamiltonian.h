#pragma once

#include <complex>
#include <vector>

enum class Pauli { I, X, Y, Z };

struct PauliTerm {
  double coefficient;
  // operators[qubit] acts on that qubit; qubit 0 is the least-significant bit.
  std::vector<Pauli> operators;
};

// Computes <state|P|state> for a normalized state vector.
double pauliExpectation(
    const std::vector<std::complex<double>> &state,
    const std::vector<Pauli> &operators
);

class Simulator;

class QubitHamiltonian {
public:
  QubitHamiltonian(int qubitCount, std::vector<PauliTerm> terms);

  int qubitCount() const;
  const std::vector<PauliTerm> &terms() const;
  double expectation(const Simulator &simulator) const;

private:
  int qubitCount_;
  std::vector<PauliTerm> terms_;
};
