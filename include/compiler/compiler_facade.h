#pragma once

#include "compiler/ast.h"
#include "compiler/code_generator.h"

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
  int line = 0;   // 1-based source line; 0 if unlocated/global
  int column = 0; // 1-based source column; 0 if unlocated
};

struct CompilerOptions {
  bool optimize = true;
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
