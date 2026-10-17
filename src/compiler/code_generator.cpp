#include "compiler/code_generator.h"
#include "compiler/circuit_dag.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

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

std::vector<Instruction> optimizeInstructions(
    const std::vector<Instruction> &instructions,
    std::size_t qubitCount
) {
  return compiler::optimizeWithDag(instructions, qubitCount);
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

void printExpression(
    std::ostream &out,
    const Expression &expression,
    int parentPrecedence,
    bool isRightChild
) {
  // If this is a binary expression, unroll left spine iteratively
  if (expression.type == ExpressionType::BinaryOperation) {
    int precedence = (expression.operation == '+' || expression.operation == '-') ? 1 : 2;
    // Parenthesize if lower precedence than parent, OR if equal precedence as the right operand
    // in left-associative syntax (e.g. 10 - (3 - 2), 24 / (4 / 2), 2 + (3 - 1)).
    bool needsParentheses = (precedence < parentPrecedence) || (isRightChild && precedence == parentPrecedence);
    if (needsParentheses) {
      out << "(";
    }

    std::vector<const Expression *> spine;
    const Expression *curr = &expression;
    while (curr->type == ExpressionType::BinaryOperation &&
           ((curr->operation == '+' || curr->operation == '-') ? 1 : 2) == precedence) {
      spine.push_back(curr);
      curr = &curr->children[0];
    }

    printExpression(out, *curr, precedence, false);

    for (auto it = spine.rbegin(); it != spine.rend(); ++it) {
      const Expression *node = *it;
      out << " " << node->operation << " ";
      printExpression(out, node->children[1], precedence, true);
    }

    if (needsParentheses) {
      out << ")";
    }
    return;
  }

  int precedence = (expression.type == ExpressionType::UnaryOperation) ? 3 : 4;
  bool needsParentheses = precedence < parentPrecedence;
  if (needsParentheses) {
    out << "(";
  }

  if (expression.type == ExpressionType::NumberLiteral ||
      expression.type == ExpressionType::Identifier) {
    out << expression.text;
  } else if (expression.type == ExpressionType::UnaryOperation) {
    out << expression.operation;
    printExpression(out, expression.children[0], precedence, false);
  }

  if (needsParentheses) {
    out << ")";
  }
}

std::string formatExpression(
    const Expression &expression,
    int parentPrecedence,
    bool isRightChild
) {
  std::ostringstream out;
  printExpression(out, expression, parentPrecedence, isRightChild);
  return out.str();
}
