#include "compiler/semantic_checker.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace {
void appendError(std::vector<std::string> &errors,
                 std::vector<SourceSpan> &spans, std::string message,
                 SourceSpan span) {
  errors.push_back(std::move(message));
  spans.push_back(span.isValid() ? span : SourceSpan{});
}

SourceSpan expressionOperationSpan(const Expression &expression) {
  return expression.operatorSpan.isValid() ? expression.operatorSpan
                                           : expression.sourceSpan;
}

SourceSpan indexedSpan(const std::vector<SourceSpan> &spans,
                       std::size_t index, const SourceSpan &fallback) {
  return index < spans.size() && spans[index].isValid() ? spans[index]
                                                        : fallback;
}

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
    std::vector<SourceSpan> &errorSpans,
    const std::string &context
) {
  if (expression.type == ExpressionType::NumberLiteral) {
    try {
      std::size_t processed = 0;
      value = std::stod(expression.text, &processed);
      if (processed != expression.text.size()) {
        appendError(errors, errorSpans,
                    "Error: invalid number literal '" + expression.text +
                        "' at line " + std::to_string(line),
                    expression.sourceSpan);
        return false;
      }
    } catch (const std::out_of_range &) {
      appendError(errors, errorSpans,
                  "Error: number literal '" + expression.text +
                      "' out of range at line " + std::to_string(line),
                  expression.sourceSpan);
      return false;
    } catch (const std::invalid_argument &) {
      appendError(errors, errorSpans,
                  "Error: invalid number literal '" + expression.text +
                      "' at line " + std::to_string(line),
                  expression.sourceSpan);
      return false;
    }
    if (!std::isfinite(value)) {
      appendError(errors, errorSpans,
                  "Error: non-finite number literal '" + expression.text +
                      "' at line " + std::to_string(line),
                  expression.sourceSpan);
      return false;
    }
    return true;
  }

  if (expression.type == ExpressionType::Identifier) {
    const int variableIndex = findNameIndex(variableNames, expression.text);
    if (variableIndex == -1) {
      appendError(errors, errorSpans,
                  "Error: undefined variable '" + expression.text +
                      "' used at line " + std::to_string(line),
                  expression.sourceSpan);
      return false;
    }
    value = variableValues[static_cast<std::size_t>(variableIndex)];
    return true;
  }

  if (expression.type == ExpressionType::UnaryOperation) {
    double childValue = 0.0;
    const bool isValid = evaluateExpression(
        expression.children[0], line, childValue, variableNames,
        variableValues, errors, errorSpans, context);
    if (!isValid) return false;
    value = -childValue;
    if (!std::isfinite(value)) {
      appendError(errors, errorSpans,
                  "Error: non-finite arithmetic result at line " +
                      std::to_string(line),
                  expressionOperationSpan(expression));
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
  bool allValid = evaluateExpression(
      *curr, line, accumulator, variableNames, variableValues, errors,
      errorSpans, context);

  for (auto it = operations.rbegin(); it != operations.rend(); ++it) {
    const Expression *opNode = *it;
    double rightValue = 0.0;
    const bool rightValid = evaluateExpression(
        opNode->children[1], line, rightValue, variableNames, variableValues,
        errors, errorSpans, context);
    if (!rightValid) allValid = false;

    if (allValid) {
      if (opNode->operation == '+') {
        accumulator += rightValue;
      } else if (opNode->operation == '-') {
        accumulator -= rightValue;
      } else if (opNode->operation == '*') {
        accumulator *= rightValue;
      } else if (opNode->operation == '/') {
        if (rightValue == 0.0) {
          appendError(errors, errorSpans,
                      "Error: division by zero in " + context + " at line " +
                          std::to_string(line),
                      expressionOperationSpan(*opNode));
          allValid = false;
        } else {
          accumulator /= rightValue;
        }
      }

      if (allValid && !std::isfinite(accumulator)) {
        appendError(errors, errorSpans,
                    "Error: non-finite arithmetic result at line " +
                        std::to_string(line),
                    expressionOperationSpan(*opNode));
        allValid = false;
      }
    }
  }

  if (!allValid) return false;
  value = accumulator;
  return true;
}

} // namespace

SemanticResult checkSemantics(const Program &program) {
  SemanticResult result;
  result.parameterValues.resize(program.statements.size(), 0.0);

  for (std::size_t i = 0; i < program.statements.size(); ++i) {
    const Statement &statement = program.statements[i];

    if (statement.type == StatementType::QubitDeclaration) {
      for (std::size_t j = 0; j < statement.names.size(); ++j) {
        if (findNameIndex(result.qubitNames, statement.names[j]) != -1) {
          appendError(result.errors, result.errorSpans,
                      "Error: qubit '" + statement.names[j] +
                          "' is already declared at line " +
                          std::to_string(statement.line),
                      indexedSpan(statement.nameSpans, j, statement.sourceSpan));
        } else {
          result.qubitNames.push_back(statement.names[j]);
        }
      }
      continue;
    }

    if (statement.type == StatementType::VariableDeclaration) {
      const bool alreadyDefined =
          findNameIndex(result.variableNames, statement.variableName) != -1;
      if (alreadyDefined) {
        appendError(result.errors, result.errorSpans,
                    "Error: variable '" + statement.variableName +
                        "' is already declared at line " +
                        std::to_string(statement.line),
                    statement.variableNameSpan);
      }

      double variableValue = 0.0;
      const bool isValid = evaluateExpression(
          statement.variableExpression, statement.line, variableValue,
          result.variableNames, result.variableValues, result.errors,
          result.errorSpans, "declaration");
      if (isValid && !alreadyDefined) {
        result.variableNames.push_back(statement.variableName);
        result.variableValues.push_back(variableValue);
      }
      continue;
    }

    if (statement.gate == "measure") {
      std::vector<std::string> seenQubits;
      for (std::size_t j = 0; j < statement.qubits.size(); ++j) {
        if (findNameIndex(seenQubits, statement.qubits[j]) != -1) {
          appendError(result.errors, result.errorSpans,
                      "Error: duplicate qubit '" + statement.qubits[j] +
                          "' in measure statement at line " +
                          std::to_string(statement.line),
                      indexedSpan(statement.qubitSpans, j, statement.sourceSpan));
        } else {
          seenQubits.push_back(statement.qubits[j]);
        }
      }
    }

    for (std::size_t j = 0; j < statement.qubits.size(); ++j) {
      if (findNameIndex(result.qubitNames, statement.qubits[j]) == -1) {
        appendError(result.errors, result.errorSpans,
                    "Error: undeclared qubit '" + statement.qubits[j] +
                        "' at line " + std::to_string(statement.line),
                    indexedSpan(statement.qubitSpans, j, statement.sourceSpan));
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
      appendError(result.errors, result.errorSpans,
                  "Error: gate '" + statement.gate + "' needs " +
                      std::to_string(requiredArity) +
                      " qubit(s) at line " + std::to_string(statement.line),
                  statement.sourceSpan);
    }

    if (statement.gate == "CNOT" && statement.qubits.size() == 2 &&
        statement.qubits[0] == statement.qubits[1]) {
      appendError(result.errors, result.errorSpans,
                  "Error: CNOT control and target must differ at line " +
                      std::to_string(statement.line),
                  indexedSpan(statement.qubitSpans, 1, statement.sourceSpan));
    }

    if (statement.gate == "RZ" && !statement.hasParameter) {
      appendError(result.errors, result.errorSpans,
                  "Error: gate 'RZ' needs a parameter at line " +
                      std::to_string(statement.line),
                  statement.gateSpan);
    }

    if (statement.gate != "RZ" && statement.hasParameter) {
      appendError(result.errors, result.errorSpans,
                  "Error: gate '" + statement.gate +
                      "' does not take a parameter at line " +
                      std::to_string(statement.line),
                  statement.gateSpan);
    }

    if (statement.hasParameter) {
      double parameterValue = 0.0;
      if (evaluateExpression(
              statement.parameter, statement.line, parameterValue,
              result.variableNames, result.variableValues, result.errors,
              result.errorSpans, "gate parameter")) {
        result.parameterValues[i] = parameterValue;
      }
    }
  }

  return result;
}
