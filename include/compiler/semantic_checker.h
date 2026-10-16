#pragma once

#include "compiler/ast.h"

#include <string>
#include <vector>

struct SemanticResult {
  std::vector<std::string> errors;
  std::vector<std::string> qubitNames;
  std::vector<std::string> variableNames;
  std::vector<double> variableValues;
  std::vector<double> parameterValues;
};

SemanticResult checkSemantics(const Program &program);
