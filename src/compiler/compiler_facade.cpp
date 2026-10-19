#include "compiler/compiler_facade.h"
#include "compiler/circuit_dag.h"

#include "compiler/lexer.h"
#include "compiler/parser.h"
#include "compiler/semantic_checker.h"

#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace compiler {

namespace {

Diagnostic makeDiagnostic(const std::string &message, SourceSpan span = {}) {
  Diagnostic diagnostic;
  diagnostic.severity = DiagnosticSeverity::Error;
  diagnostic.message = message;
  if (span.isValid()) {
    diagnostic.sourceSpan = span;
    diagnostic.line = span.start.line;
    diagnostic.column = span.start.column;
  }
  return diagnostic;
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
  } catch (const SourceLocatedError &error) {
    result.diagnostics.push_back(makeDiagnostic(error.what(), error.span()));
    result.success = false;
    return result;
  } catch (const std::exception &error) {
    result.diagnostics.push_back(makeDiagnostic(error.what()));
    result.success = false;
    return result;
  }

  // 2. Syntax analysis (Parsing)
  try {
    Parser parser(std::move(tokens));
    result.program = parser.parse();
  } catch (const SourceLocatedError &error) {
    result.diagnostics.push_back(makeDiagnostic(error.what(), error.span()));
    result.success = false;
    return result;
  } catch (const std::exception &error) {
    result.diagnostics.push_back(makeDiagnostic(error.what()));
    result.success = false;
    return result;
  }

  // 3. Semantic checking
  SemanticResult semanticResult = checkSemantics(result.program);
  if (!semanticResult.errors.empty()) {
    result.diagnostics.reserve(semanticResult.errors.size());
    for (std::size_t i = 0; i < semanticResult.errors.size(); ++i) {
      const SourceSpan span = i < semanticResult.errorSpans.size()
                                  ? semanticResult.errorSpans[i]
                                  : SourceSpan{};
      result.diagnostics.push_back(
          makeDiagnostic(semanticResult.errors[i], span));
    }
    result.success = false;
    return result;
  }

  // 4. Code generation & Optimization
  result.qubitCount = semanticResult.qubitNames.size();
  result.qubitNames = semanticResult.qubitNames;
  result.rawInstructions = generateInstructions(result.program, semanticResult);

  CircuitDag rawDag = CircuitDag::fromInstructions(result.rawInstructions, result.qubitCount);
  result.rawGateCount = rawDag.gateCount(InstructionType::ApplyGate);
  result.rawTwoQubitCount = rawDag.twoQubitGateCount();
  result.rawDepth = rawDag.criticalPathDepth();

  if (options.optimize) {
    CircuitDag optDag = rawDag.optimize();
    result.instructions = optDag.toInstructions();
    result.optimizedGateCount = optDag.gateCount(InstructionType::ApplyGate);
    result.optimizedTwoQubitCount = optDag.twoQubitGateCount();
    result.optimizedDepth = optDag.criticalPathDepth();
  } else {
    result.instructions = result.rawInstructions;
    result.optimizedGateCount = result.rawGateCount;
    result.optimizedTwoQubitCount = result.rawTwoQubitCount;
    result.optimizedDepth = result.rawDepth;
  }

  // 5. Topology Routing (if requested)
  if (options.targetTopology.has_value()) {
    RoutingResult routing = routeCircuit(result.instructions, result.qubitCount, *options.targetTopology);
    if (!routing.success) {
      Diagnostic diag;
      diag.severity = DiagnosticSeverity::Error;
      diag.message = "Topology routing failed: " + routing.errorMessage;
      diag.line = 0;
      diag.column = 0;
      result.diagnostics.push_back(diag);
      result.success = false;
      return result;
    }

    result.routed = true;
    result.routedInstructions = std::move(routing.routedInstructions);
    result.routedGateCount = routing.routedGateCount;
    result.routedTwoQubitCount = routing.routedTwoQubitCount;
    result.insertedSwapCount = routing.insertedSwapCount;
    result.routedDepth = routing.routedDepth;
    result.routingOverhead = routing.routingOverhead;
    result.initialMapping = routing.initialMapping;
    result.finalMapping = routing.finalMapping;
  }

  result.success = true;
  return result;
}

} // namespace compiler
