#include "chemistry/hamiltonian.h"

#include "simulator/simulator.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace {

std::size_t expectedDimension(std::size_t qubitCount) {
  if (qubitCount >= sizeof(std::size_t) * 8) {
    throw std::invalid_argument("Pauli operator has too many qubits");
  }
  return std::size_t{1} << static_cast<unsigned int>(qubitCount);
}

} // namespace

double pauliExpectation(
    const std::vector<std::complex<double>> &state,
    const std::vector<Pauli> &operators
) {
  if (state.size() != expectedDimension(operators.size())) {
    throw std::invalid_argument(
        "State-vector size does not match Pauli operator qubit count"
    );
  }

  double normSquared = 0.0;
  for (const std::complex<double> &amplitude : state) {
    if (!std::isfinite(amplitude.real()) || !std::isfinite(amplitude.imag())) {
      throw std::invalid_argument("State vector contains a non-finite amplitude");
    }
    normSquared += std::norm(amplitude);
  }
  if (!std::isfinite(normSquared) || std::abs(normSquared - 1.0) > 1e-10) {
    throw std::invalid_argument("State vector must be normalized");
  }

  std::complex<double> expectation(0.0, 0.0);
  for (std::size_t basis = 0; basis < state.size(); ++basis) {
    std::size_t output = basis;
    std::complex<double> phase(1.0, 0.0);
    for (std::size_t qubit = 0; qubit < operators.size(); ++qubit) {
      const std::size_t mask = std::size_t{1} << static_cast<unsigned int>(qubit);
      const bool bitIsOne = (basis & mask) != 0;
      switch (operators[qubit]) {
      case Pauli::I:
        break;
      case Pauli::X:
        output ^= mask;
        break;
      case Pauli::Y:
        output ^= mask;
        phase *= bitIsOne ? std::complex<double>(0.0, -1.0)
                          : std::complex<double>(0.0, 1.0);
        break;
      case Pauli::Z:
        if (bitIsOne) {
          phase = -phase;
        }
        break;
      default:
        throw std::invalid_argument("Unsupported Pauli operator");
      }
    }
    expectation += std::conj(state[output]) * phase * state[basis];
  }

  if (std::abs(expectation.imag()) > 1e-10) {
    throw std::runtime_error("Hermitian Pauli expectation is not real");
  }
  return expectation.real();
}

QubitHamiltonian::QubitHamiltonian(int qubitCount, std::vector<PauliTerm> terms)
    : qubitCount_(qubitCount), terms_(std::move(terms)) {
  if (qubitCount_ < 0 || qubitCount_ > MAX_SIMULATOR_QUBITS) {
    throw std::invalid_argument("Hamiltonian qubit count is out of range");
  }
  for (const PauliTerm &term : terms_) {
    if (!std::isfinite(term.coefficient)) {
      throw std::invalid_argument("Hamiltonian coefficient must be finite");
    }
    if (term.operators.size() != static_cast<std::size_t>(qubitCount_)) {
      throw std::invalid_argument(
          "Hamiltonian Pauli term has an invalid qubit count"
      );
    }
  }
}

int QubitHamiltonian::qubitCount() const { return qubitCount_; }

const std::vector<PauliTerm> &QubitHamiltonian::terms() const { return terms_; }

double QubitHamiltonian::expectation(const Simulator &simulator) const {
  double energy = 0.0;
  for (const PauliTerm &term : terms_) {
    energy += term.coefficient * simulator.expectation(term.operators);
  }
  return energy;
}
