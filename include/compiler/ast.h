#pragma once

#include <iosfwd>
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
  ExpressionType type = ExpressionType::NumberLiteral;
  std::string text = "";
  char operation = '\0';
  // binary expression -> two children; unary minus -> one children.
  std::vector<Expression> children = {};
};

struct Statement {
  StatementType type = StatementType::QubitDeclaration;
  int line = 0;
  std::vector<std::string> names = {};
  std::string variableName = "";
  Expression variableExpression = {};
  std::string gate = "";
  std::vector<std::string> qubits = {};
  bool hasParameter = false;
  Expression parameter = {};
};

struct Program {
  std::vector<Statement> statements;
};

void printExpression(
    std::ostream &out,
    const Expression &expression,
    int parentPrecedence = 0,
    bool isRightChild = false
);

std::string formatExpression(
    const Expression &expression,
    int parentPrecedence = 0,
    bool isRightChild = false
);
