#pragma once

#include "compiler/ast.h"
#include "compiler/code_generator.h"
#include "compiler/lexer.h"
#include "compiler/parser.h"
#include "compiler/semantic_checker.h"
#include "simulator/simulator.h"

#include <complex>
#include <string>
#include <vector>

namespace test {

extern int assertions;

void require(bool condition, const std::string &message);

Program parse(const std::string &source);

SemanticResult check(const std::string &source, Program &program);

std::string run(const std::string &source);

using Complex = std::complex<double>;

Instruction makeInstruction(InstructionType type, const std::string &gate,
                            std::vector<int> qubits, double parameter = 0.0);

bool sameInstructions(const std::vector<Instruction> &left,
                      const std::vector<Instruction> &right);

void requireSameState(const std::vector<Complex> &left,
                      const std::vector<Complex> &right,
                      const std::string &message);

} // namespace test
