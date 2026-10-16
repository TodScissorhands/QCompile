#include "compiler/code_generator.h"

#include <iostream>

namespace {

int findQubitIndex(
    const std::vector<std::string> &names,
    const std::string &name
) {
  for (int i = 0; i < static_cast<int>(names.size()); i = i + 1) {
    if (names[i] == name) {
      return i;
    }
  }

  return -1;
}

} // namespace

namespace {

bool canCancel(const Instruction &left, const Instruction &right,
               std::size_t qubitCount) {
  if (left.type != InstructionType::ApplyGate ||
      right.type != InstructionType::ApplyGate ||
      left.gateName != right.gateName ||
      left.qubitIndices != right.qubitIndices ||
      left.parameterValue != 0.0 ||
      right.parameterValue != 0.0) {
    return false;
  }

  if (left.gateName == "H" || left.gateName == "X" ||
      left.gateName == "Z") {
    return left.qubitIndices.size() == 1 &&
           left.qubitIndices[0] >= 0 &&
           static_cast<std::size_t>(left.qubitIndices[0]) < qubitCount;
  }

  if (left.gateName == "CNOT") {
    return left.qubitIndices.size() == 2 &&
           left.qubitIndices[0] >= 0 &&
           left.qubitIndices[1] >= 0 &&
           static_cast<std::size_t>(left.qubitIndices[0]) < qubitCount &&
           static_cast<std::size_t>(left.qubitIndices[1]) < qubitCount &&
           left.qubitIndices[0] != left.qubitIndices[1];
  }

  return false;
}

} // namespace

std::vector<Instruction> optimizeInstructions(
    const std::vector<Instruction> &instructions,
    std::size_t qubitCount
) {
  std::vector<Instruction> optimized;
  optimized.reserve(instructions.size());

  for (const Instruction &instruction : instructions) {
    if (!optimized.empty() &&
        canCancel(optimized.back(), instruction, qubitCount)) {
      optimized.pop_back();
    } else {
      optimized.push_back(instruction);
    }
  }

  return optimized;
}

std::vector<Instruction> generateInstructions(
    const Program &program,
    const SemanticResult &result
) {
  std::vector<Instruction> instructions;

  for (int i = 0; i < static_cast<int>(program.statements.size()); i = i + 1) {
    const Statement &statement = program.statements[i];

    if (statement.type == StatementType::QubitDeclaration ||
        statement.type == StatementType::VariableDeclaration) {
      continue;
    }

    Instruction instruction;
    instruction.gateName = statement.gate;
    instruction.parameterValue = result.parameterValues[i];

    if (statement.gate == "measure") {
      instruction.type = InstructionType::Measure;
    } else {
      instruction.type = InstructionType::ApplyGate;
    }

    for (int j = 0; j < static_cast<int>(statement.qubits.size()); j = j + 1) {
      instruction.qubitIndices.push_back(
          findQubitIndex(result.qubitNames, statement.qubits[j]));
    }

    instructions.push_back(instruction);
  }

  return instructions;
}

void printInstructions(const std::vector<Instruction> &instructions) {
  std::cout << "\nInstructions:\n";

  for (int i = 0; i < static_cast<int>(instructions.size()); i = i + 1) {
    const Instruction &instruction = instructions[i];

    if (instruction.type == InstructionType::ApplyGate) {
      std::cout << "APPLY " << instruction.gateName;

      if (instruction.gateName == "RZ") {
        std::cout << "(" << instruction.parameterValue << ")";
      }
    } else {
      std::cout << "MEASURE";
    }

    std::cout << " [";
    for (int j = 0; j < static_cast<int>(instruction.qubitIndices.size());
         j = j + 1) {
      if (j != 0) {
        std::cout << ", ";
      }
      std::cout << instruction.qubitIndices[j];
    }
    std::cout << "]\n";
  }
}
