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
  // Kept parallel to errors; invalid spans mark unlocated/global errors.
  std::vector<SourceSpan> errorSpans;
};

SemanticResult checkSemantics(const Program &program);
