#pragma once

#include <string>
#include <vector>

// The enum tags the plain Statement record without using a class hierarchy.
enum class StatementType { QubitDeclaration, GateOperation };

enum class ExpressionType {
  NumberLiteral,
  Identifier,
  BinaryOperation,
  UnaryOperation
};

struct Expression {
  ExpressionType type;
  std::string text;
  char operation;
  // A binary expression has two children; unary minus has one.
  std::vector<Expression> children;
};

struct Statement {
  StatementType type;
  std::vector<std::string> names;
  std::string gate;
  std::vector<std::string> qubits;
  bool hasParameter;
  Expression parameter;
};

struct Program {
  std::vector<Statement> statements;
};
