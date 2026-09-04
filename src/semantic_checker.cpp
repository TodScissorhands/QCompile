#include "semantic_checker.h"

#include <stdexcept>

namespace {

int findNameIndex(
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

bool evaluateExpression(
    const Expression &expression,
    int line,
    double &value,
    const std::vector<std::string> &variableNames,
    const std::vector<double> &variableValues,
    std::vector<std::string> &errors
) {
  if (expression.type == ExpressionType::NumberLiteral) {
    value = std::stod(expression.text);
    return true;
  }

  if (expression.type == ExpressionType::Identifier) {
    int variableIndex = findNameIndex(variableNames, expression.text);
    if (variableIndex == -1) {
      errors.push_back("Error: undefined variable '" + expression.text +
                       "' used at line " + std::to_string(line));
      return false;
    }

    value = variableValues[variableIndex];
    return true;
  }

  if (expression.type == ExpressionType::UnaryOperation) {
    double childValue = 0.0;
    bool isValid = evaluateExpression(
        expression.children[0],
        line,
        childValue,
        variableNames,
        variableValues,
        errors
    );
    value = -childValue;
    return isValid;
  }

  double leftValue = 0.0;
  double rightValue = 0.0;
  bool leftIsValid = evaluateExpression(
      expression.children[0],
      line,
      leftValue,
      variableNames,
      variableValues,
      errors
  );
  bool rightIsValid = evaluateExpression(
      expression.children[1],
      line,
      rightValue,
      variableNames,
      variableValues,
      errors
  );

  if (!leftIsValid || !rightIsValid) {
    return false;
  }

  if (expression.operation == '+') {
    value = leftValue + rightValue;
  } else if (expression.operation == '-') {
    value = leftValue - rightValue;
  } else if (expression.operation == '*') {
    value = leftValue * rightValue;
  } else if (expression.operation == '/') {
    if (rightValue == 0.0) {
      errors.push_back("Error: division by zero in gate parameter at line " +
                       std::to_string(line));
      return false;
    }
    value = leftValue / rightValue;
  }

  return true;
}

} // namespace

SemanticResult checkSemantics(const Program &program) {
  SemanticResult result;
  result.parameterValues.resize(program.statements.size(), 0.0);

  for (int i = 0; i < static_cast<int>(program.statements.size()); i = i + 1) {
    const Statement &statement = program.statements[i];

    if (statement.type == StatementType::QubitDeclaration) {
      for (int j = 0; j < static_cast<int>(statement.names.size()); j = j + 1) {
        if (findNameIndex(result.qubitNames, statement.names[j]) != -1) {
          result.errors.push_back("Error: qubit '" + statement.names[j] +
                                  "' is already declared at line " +
                                  std::to_string(statement.line));
        } else {
          result.qubitNames.push_back(statement.names[j]);
        }
      }
      continue;
    }

    if (statement.type == StatementType::VariableDeclaration) {
      bool alreadyDefined =
          findNameIndex(result.variableNames, statement.variableName) != -1;
      if (alreadyDefined) {
        result.errors.push_back("Error: variable '" + statement.variableName +
                                "' is already declared at line " +
                                std::to_string(statement.line));
      }

      double variableValue = 0.0;
      bool isValid = evaluateExpression(
          statement.variableExpression,
          statement.line,
          variableValue,
          result.variableNames,
          result.variableValues,
          result.errors
      );

      if (isValid && !alreadyDefined) {
        result.variableNames.push_back(statement.variableName);
        result.variableValues.push_back(variableValue);
      }
      continue;
    }

    for (int j = 0; j < static_cast<int>(statement.qubits.size()); j = j + 1) {
      if (findNameIndex(result.qubitNames, statement.qubits[j]) == -1) {
        result.errors.push_back("Error: undeclared qubit '" +
                                statement.qubits[j] + "' at line " +
                                std::to_string(statement.line));
      }
    }

    int requiredArity = -1;
    if (statement.gate == "H" || statement.gate == "X" ||
        statement.gate == "Z" || statement.gate == "RZ") {
      requiredArity = 1;
    } else if (statement.gate == "CNOT") {
      requiredArity = 2;
    }

    if (requiredArity != -1 &&
        static_cast<int>(statement.qubits.size()) != requiredArity) {
      result.errors.push_back("Error: gate '" + statement.gate + "' needs " +
                              std::to_string(requiredArity) +
                              " qubit(s) at line " +
                              std::to_string(statement.line));
    }

    if (statement.gate == "CNOT" && statement.qubits.size() == 2 &&
        statement.qubits[0] == statement.qubits[1]) {
      result.errors.push_back(
          "Error: CNOT control and target must differ at line " +
          std::to_string(statement.line));
    }

    if (statement.gate == "RZ" && !statement.hasParameter) {
      result.errors.push_back("Error: gate 'RZ' needs a parameter at line " +
                              std::to_string(statement.line));
    }

    if (statement.gate != "RZ" && statement.hasParameter) {
      result.errors.push_back("Error: gate '" + statement.gate +
                              "' does not take a parameter at line " +
                              std::to_string(statement.line));
    }

    if (statement.hasParameter) {
      double parameterValue = 0.0;
      if (evaluateExpression(
          statement.parameter,
          statement.line,
          parameterValue,
          result.variableNames,
          result.variableValues,
          result.errors
      )) {
        result.parameterValues[i] = parameterValue;
      }
    }
  }

  return result;
}
