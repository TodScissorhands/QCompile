#pragma once

#include "chemistry/hamiltonian.h"
#include "compiler/code_generator.h"

#include <cstddef>
#include <vector>

namespace compiler {

// Decomposes single-qubit RY(theta) into native RZ and H gates:
// RY(theta) = RZ(-pi/2) H RZ(theta) H RZ(pi/2)
void appendDecomposedRY(std::vector<Instruction> &circuit, int target, double angle);

// Standard parameterized Hardware-Efficient Ansatz (RY-CX linear entangler)
// Parameter count: qubitCount * (depth + 1)
std::vector<Instruction> createHardwareEfficientAnsatz(
    int qubitCount,
    int depth,
    const std::vector<double> &parameters
);

// LiH active-space 4-qubit singlet trial circuit (Jordan-Wigner, 2 parameters)
// Starts in HF determinant |0101> and executes Givens singlet rotations
std::vector<Instruction> createLiHActiveSpaceAnsatz(double theta, double phi);

// LiH active-space 2-qubit singlet trial circuit (Parity-reduced (1,1) sector, 1 parameter)
std::vector<Instruction> createLiHParityAnsatz(double theta);

// Evaluates the expectation value of a QubitHamiltonian for a given circuit
double evaluateAnsatzEnergy(
    const std::vector<Instruction> &circuit,
    int qubitCount,
    const QubitHamiltonian &hamiltonian
);

} // namespace compiler
