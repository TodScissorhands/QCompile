#pragma once

#include "ast.h"
#include "semantic_checker.h"

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
