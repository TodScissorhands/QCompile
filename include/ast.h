#pragma once

#include <string>
#include <vector>

// The enum tags the plain Statement record without using a class hierarchy.
enum class StatementType {
  QubitDeclaration,
  VariableDeclaration,
  GateOperation
};

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
  // binary expression -> two children; unary minus -> one children.
  std::vector<Expression> children;
};

struct Statement {
  StatementType type;
  int line;
  std::vector<std::string> names;
  std::string variableName;
  Expression variableExpression;
  std::string gate;
  std::vector<std::string> qubits;
  bool hasParameter;
  Expression parameter;
};

struct Program {
  std::vector<Statement> statements;
};
