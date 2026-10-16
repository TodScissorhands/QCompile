#pragma once

#include "compiler/ast.h"
#include "compiler/semantic_checker.h"

#include <cstddef>
#include <string>
#include <vector>

enum class InstructionType { ApplyGate, Measure };

struct Instruction {
  InstructionType type;
  std::string gateName;
  std::vector<int> qubitIndices;
  double parameterValue;
};

std::vector<Instruction> generateInstructions(
    const Program &program,
    const SemanticResult &result
);

void printInstructions(const std::vector<Instruction> &instructions);
std::vector<Instruction> optimizeInstructions(
    const std::vector<Instruction> &instructions,
    std::size_t qubitCount
);
