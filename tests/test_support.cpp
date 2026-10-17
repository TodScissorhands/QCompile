#include "test_support.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace test {

int assertions = 0;

void require(bool condition, const std::string &message) {
  ++assertions;
  if (!condition) {
    throw std::runtime_error(message);
  }
}

Program parse(const std::string &source) {
  std::vector<Token> tokens = Lexer(source).tokenize();
  return Parser(tokens).parse();
}

SemanticResult check(const std::string &source, Program &program) {
  program = parse(source);
  return checkSemantics(program);
}

std::string run(const std::string &source) {
  Program program = parse(source);
  SemanticResult semantic = checkSemantics(program);
  require(semantic.errors.empty(), "simulator input must pass semantic checking");
  std::vector<Instruction> instructions = generateInstructions(program, semantic);
  std::ostringstream output;
  std::streambuf *original = std::cout.rdbuf(output.rdbuf());
  try {
    Simulator simulator(static_cast<int>(semantic.qubitNames.size()));
    simulator.run(instructions);
  } catch (...) {
    std::cout.rdbuf(original);
    throw;
  }
  std::cout.rdbuf(original);
  return output.str();
}

Instruction makeInstruction(InstructionType type, const std::string &gate,
                            std::vector<int> qubits, double parameter) {
  Instruction instruction;
  instruction.type = type;
  instruction.gateName = gate;
  instruction.qubitIndices = std::move(qubits);
  instruction.parameterValue = parameter;
  return instruction;
}

bool sameInstructions(const std::vector<Instruction> &left,
                      const std::vector<Instruction> &right) {
  if (left.size() != right.size()) return false;
  for (std::size_t index = 0; index < left.size(); ++index) {
    if (left[index].type != right[index].type ||
        left[index].gateName != right[index].gateName ||
        left[index].qubitIndices != right[index].qubitIndices) {
      return false;
    }
    const double lp = left[index].parameterValue;
    const double rp = right[index].parameterValue;
    if (std::isnan(lp) && std::isnan(rp)) {
      continue;
    }
    if (lp != rp) {
      return false;
    }
  }
  return true;
}

void requireSameState(const std::vector<Complex> &left,
                      const std::vector<Complex> &right,
                      const std::string &message) {
  // State-vector equivalence tolerance: absolute complex-amplitude error < 1e-12.
  require(left.size() == right.size(), message + " (dimension)");
  for (std::size_t index = 0; index < left.size(); ++index) {
    require(std::abs(left[index] - right[index]) < 1e-12,
            message + " (amplitude " + std::to_string(index) + ")");
  }
}

} // namespace test
