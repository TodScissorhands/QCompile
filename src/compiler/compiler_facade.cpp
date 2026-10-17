#include "compiler/compiler_facade.h"

#include "compiler/lexer.h"
#include "compiler/parser.h"
#include "compiler/semantic_checker.h"

#include <cctype>
#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace compiler {

namespace {

int extractLineNumber(const std::string &message) {
  const std::string pattern = " at line ";
  const std::size_t pos = message.rfind(pattern);
  if (pos != std::string::npos) {
    const std::size_t start = pos + pattern.size();
    if (start < message.size() &&
        std::isdigit(static_cast<unsigned char>(message[start]))) {
      try {
        std::size_t processed = 0;
        int line = std::stoi(message.substr(start), &processed);
        return line;
      } catch (...) {
        return 0;
      }
    }
  }
  return 0;
}

} // namespace

CompilationResult compileProgram(
    const std::string &sourceText,
    const CompilerOptions &options
) {
  CompilationResult result;

  // 1. Lexical analysis
  std::vector<Token> tokens;
  try {
    Lexer lexer(sourceText);
    tokens = lexer.tokenize();
  } catch (const std::exception &error) {
    const std::string msg = error.what();
    result.diagnostics.push_back({
        DiagnosticSeverity::Error,
        msg,
        extractLineNumber(msg),
        0
    });
    result.success = false;
    return result;
  }

  // 2. Syntax analysis (Parsing)
  try {
    Parser parser(std::move(tokens));
    result.program = parser.parse();
  } catch (const std::exception &error) {
    const std::string msg = error.what();
    result.diagnostics.push_back({
        DiagnosticSeverity::Error,
        msg,
        extractLineNumber(msg),
        0
    });
    result.success = false;
    return result;
  }

  // 3. Semantic checking
  SemanticResult semanticResult = checkSemantics(result.program);
  if (!semanticResult.errors.empty()) {
    result.diagnostics.reserve(semanticResult.errors.size());
    for (const auto &err : semanticResult.errors) {
      result.diagnostics.push_back({
          DiagnosticSeverity::Error,
          err,
          extractLineNumber(err),
          0
      });
    }
    result.success = false;
    return result;
  }

  // 4. Code generation & Optimization
  result.qubitCount = semanticResult.qubitNames.size();
  result.qubitNames = semanticResult.qubitNames;
  result.rawInstructions = generateInstructions(result.program, semanticResult);

  if (options.optimize) {
    result.instructions = optimizeInstructions(result.rawInstructions, result.qubitCount);
  } else {
    result.instructions = result.rawInstructions;
  }

  result.success = true;
  return result;
}

} // namespace compiler
