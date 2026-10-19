#pragma once

#include "compiler/ast.h"
#include "compiler/code_generator.h"
#include "compiler/coupling_map.h"
#include "compiler/router.h"

#include <optional>
#include <cstddef>
#include <string>
#include <vector>

namespace compiler {

enum class DiagnosticSeverity {
  Error,
  Warning,
  Info
};

struct Diagnostic {
  DiagnosticSeverity severity = DiagnosticSeverity::Error;
  std::string message;
  int line = 0;   // 1-based; 0 if unlocated/global
  int column = 0; // 1-based UTF-8 byte column; 0 if unlocated/global
  SourceSpan sourceSpan; // Half-open, with a zero-length insertion point allowed
};

struct CompilerOptions {
  bool optimize = true;
  std::optional<CouplingMap> targetTopology = std::nullopt;
};
struct CompilationResult {
  bool success = false;
  std::vector<Diagnostic> diagnostics;
  std::vector<Instruction> rawInstructions;
  std::vector<Instruction> instructions; // optimized if options.optimize == true
  std::size_t qubitCount = 0;
  Program program;
  std::vector<std::string> qubitNames;

  // Milestone 6 Circuit DAG & scheduling metrics
  std::size_t rawGateCount = 0;
  std::size_t optimizedGateCount = 0;
  std::size_t rawTwoQubitCount = 0;
  std::size_t optimizedTwoQubitCount = 0;
  std::size_t rawDepth = 0;
  std::size_t optimizedDepth = 0;

  // Milestone 7 Hardware Topology & Routing metrics
  bool routed = false;
  std::vector<Instruction> routedInstructions;
  std::size_t routedGateCount = 0;
  std::size_t routedTwoQubitCount = 0;
  std::size_t insertedSwapCount = 0;
  std::size_t routedDepth = 0;
  double routingOverhead = 0.0;
  std::optional<QubitMapping> initialMapping = std::nullopt;
  std::optional<QubitMapping> finalMapping = std::nullopt;

  bool hasErrors() const { return !success; }
  std::vector<std::string> errorMessages() const {
    std::vector<std::string> result;
    result.reserve(diagnostics.size());
    for (const auto &diag : diagnostics) {
      result.push_back(diag.message);
    }
    return result;
  }
};

// Programmatic compiler entry point.
// Reuses the established Lexer, Parser, SemanticChecker, and CodeGenerator pipelines.
// Emits no routine output to stdout or stderr; never terminates the process.
CompilationResult compileProgram(
    const std::string &sourceText,
    const CompilerOptions &options = {}
);

} // namespace compiler
