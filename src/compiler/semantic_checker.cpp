#include "compiler/semantic_checker.h"

#include <cmath>
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
    std::vector<std::string> &errors,
    const std::string &context
) {
  if (expression.type == ExpressionType::NumberLiteral) {
    try {
      std::size_t processed = 0;
      value = std::stod(expression.text, &processed);
      if (processed != expression.text.size()) {
        errors.push_back("Error: invalid number literal '" + expression.text +
                         "' at line " + std::to_string(line));
        return false;
      }
    } catch (const std::out_of_range &) {
      errors.push_back("Error: number literal '" + expression.text +
                       "' out of range at line " + std::to_string(line));
      return false;
    } catch (const std::invalid_argument &) {
      errors.push_back("Error: invalid number literal '" + expression.text +
                       "' at line " + std::to_string(line));
      return false;
    }
    if (!std::isfinite(value)) {
      errors.push_back("Error: non-finite number literal '" + expression.text +
                       "' at line " + std::to_string(line));
      return false;
    }
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
        errors,
        context
    );
    if (!isValid) {
      return false;
    }
    value = -childValue;
    if (!std::isfinite(value)) {
      errors.push_back("Error: non-finite arithmetic result at line " +
                       std::to_string(line));
      return false;
    }
    return true;
  }

  // Binary expression: unroll left-associative spines iteratively to prevent stack overflow.
  std::vector<const Expression *> operations;
  const Expression *curr = &expression;
  while (curr->type == ExpressionType::BinaryOperation) {
    operations.push_back(curr);
    curr = &curr->children[0];
  }

  double accumulator = 0.0;
  bool accumValid = evaluateExpression(
      *curr,
      line,
      accumulator,
      variableNames,
      variableValues,
      errors,
      context
  );

  bool allValid = accumValid;

  for (auto it = operations.rbegin(); it != operations.rend(); ++it) {
    const Expression *opNode = *it;
    double rightValue = 0.0;
    bool rightValid = evaluateExpression(
        opNode->children[1],
        line,
        rightValue,
        variableNames,
        variableValues,
        errors,
        context
    );

    if (!rightValid) {
      allValid = false;
    }

    if (allValid) {
      if (opNode->operation == '+') {
        accumulator = accumulator + rightValue;
      } else if (opNode->operation == '-') {
        accumulator = accumulator - rightValue;
      } else if (opNode->operation == '*') {
        accumulator = accumulator * rightValue;
      } else if (opNode->operation == '/') {
        if (rightValue == 0.0) {
          errors.push_back("Error: division by zero in " + context + " at line " +
                           std::to_string(line));
          allValid = false;
        } else {
          accumulator = accumulator / rightValue;
        }
      }

      if (allValid && !std::isfinite(accumulator)) {
        errors.push_back("Error: non-finite arithmetic result at line " +
                         std::to_string(line));
        allValid = false;
      }
    }
  }

  if (!allValid) {
    return false;
  }

  value = accumulator;
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
          result.errors,
          "declaration"
      );

      if (isValid && !alreadyDefined) {
        result.variableNames.push_back(statement.variableName);
        result.variableValues.push_back(variableValue);
      }
      continue;
    }
    if (statement.gate == "measure") {
      std::vector<std::string> seenQubits;
      for (int j = 0; j < static_cast<int>(statement.qubits.size()); j = j + 1) {
        if (findNameIndex(seenQubits, statement.qubits[j]) != -1) {
          result.errors.push_back("Error: duplicate qubit '" +
                                  statement.qubits[j] +
                                  "' in measure statement at line " +
                                  std::to_string(statement.line));
        } else {
          seenQubits.push_back(statement.qubits[j]);
        }
      }
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
          result.errors,
          "gate parameter"
      )) {
        result.parameterValues[i] = parameterValue;
      }
    }
  }

  return result;
}
