#include "test_support.h"
#include "compiler/compiler_facade.h"

#include <iostream>
#include <sstream>
#include <string>

namespace test {

void testCompilerFacade() {
  // 1. Successful compilation of a standard Bell state circuit
  {
    const std::string source =
        "# Prepare a Bell state\n"
        "qubit q0, q1\n"
        "H q0\n"
        "CNOT q0, q1\n"
        "measure q0, q1\n";

    // Capture stdout and stderr to guarantee no unsolicited output
    std::stringstream outCapture;
    std::stringstream errCapture;
    auto *oldCout = std::cout.rdbuf(outCapture.rdbuf());
    auto *oldCerr = std::cerr.rdbuf(errCapture.rdbuf());

    compiler::CompilationResult result = compiler::compileProgram(source);

    std::cout.rdbuf(oldCout);
    std::cerr.rdbuf(oldCerr);

    require(result.success, "valid Bell state program compiles successfully");
    require(!result.hasErrors(), "hasErrors() is false on valid compilation");
    require(result.diagnostics.empty(), "no diagnostics on clean compilation");
    require(result.qubitCount == 2, "declared qubit count is 2");
    require(result.qubitNames.size() == 2 &&
                result.qubitNames[0] == "q0" && result.qubitNames[1] == "q1",
            "qubit names preserved in declaration order");
    require(result.rawInstructions.size() == 3, "three raw instructions generated");
    require(result.instructions.size() == 3, "three optimized instructions preserved");
    require(outCapture.str().empty(), "facade writes zero output to stdout on success");
    require(errCapture.str().empty(), "facade writes zero output to stderr on success");
  }

  // 2. Optimization toggle verification
  {
    const std::string redundantCircuit =
        "qubit q0\n"
        "H q0\n"
        "H q0\n"
        "X q0\n";

    compiler::CompilerOptions unoptOptions;
    unoptOptions.optimize = false;
    compiler::CompilationResult unoptResult =
        compiler::compileProgram(redundantCircuit, unoptOptions);
    require(unoptResult.success, "unoptimized compilation succeeds");
    require(unoptResult.rawInstructions.size() == 3, "raw instruction count is 3");
    require(unoptResult.instructions.size() == 3, "unoptimized instruction count is 3");

    compiler::CompilerOptions optOptions;
    optOptions.optimize = true;
    compiler::CompilationResult optResult =
        compiler::compileProgram(redundantCircuit, optOptions);
    require(optResult.success, "optimized compilation succeeds");
    require(optResult.rawInstructions.size() == 3, "raw instructions still retained");
    require(optResult.instructions.size() == 1,
            "adjacent H-H cancelled in optimized instructions");
    require(optResult.instructions[0].gateName == "X",
            "remaining instruction is X");
  }

  // 3. Empty program and whitespace/comment handling
  {
    compiler::CompilationResult emptyResult = compiler::compileProgram("");
    require(emptyResult.success, "empty program compiles successfully");
    require(emptyResult.qubitCount == 0, "empty program has 0 qubits");
    require(emptyResult.instructions.empty(), "empty program has 0 instructions");
    require(emptyResult.diagnostics.empty(), "empty program has no diagnostics");

    compiler::CompilationResult whitespaceResult =
        compiler::compileProgram("   \n\n  # standalone comment\n\n  ");
    require(whitespaceResult.success, "whitespace-only program compiles successfully");
    require(whitespaceResult.qubitCount == 0, "whitespace program has 0 qubits");
    require(whitespaceResult.instructions.empty(), "whitespace program has 0 instructions");
  }

  // 4. Lexical diagnostic reporting (unexpected character)
  {
    const std::string invalidLex =
        "qubit q0\n"
        "@invalid\n"
        "H q0\n";

    std::stringstream outCapture;
    std::stringstream errCapture;
    auto *oldCout = std::cout.rdbuf(outCapture.rdbuf());
    auto *oldCerr = std::cerr.rdbuf(errCapture.rdbuf());

    compiler::CompilationResult result = compiler::compileProgram(invalidLex);

    std::cout.rdbuf(oldCout);
    std::cerr.rdbuf(oldCerr);

    require(!result.success, "lexical error yields success = false");
    require(result.hasErrors(), "hasErrors() is true on lexical error");
    require(!result.diagnostics.empty(), "diagnostic emitted for lexical error");
    require(result.diagnostics[0].line == 2, "lexical error correctly localized to line 2");
    require(result.diagnostics[0].message.find("unexpected character '@'") != std::string::npos,
            "diagnostic message details unexpected character");
    require(outCapture.str().empty(), "facade writes zero output to stdout on lexical error");
    require(errCapture.str().empty(), "facade writes zero output to stderr on lexical error");
  }

  // 5. Parser diagnostic reporting (syntax error)
  {
    const std::string invalidSyntax =
        "qubit q0\n"
        "let x = 2 +\n"
        "H q0\n";

    compiler::CompilationResult result = compiler::compileProgram(invalidSyntax);
    require(!result.success, "syntax error yields success = false");
    require(!result.diagnostics.empty(), "diagnostic emitted for syntax error");
    require(result.diagnostics[0].line == 2, "syntax error correctly localized to line 2");
    require(result.instructions.empty(), "no instructions generated on syntax error");
  }

  // 6. Semantic diagnostic reporting (multiple semantic errors)
  {
    const std::string semanticErrors =
        "qubit q0\n"
        "H undeclared\n"
        "measure q0, q0\n";

    compiler::CompilationResult result = compiler::compileProgram(semanticErrors);
    require(!result.success, "semantic error yields success = false");
    require(result.diagnostics.size() >= 2, "all semantic errors collected in diagnostics");
    require(result.diagnostics[0].line == 2, "undeclared qubit localized to line 2");
    require(result.diagnostics[1].line == 3, "duplicate measure qubit localized to line 3");
    require(result.diagnostics[0].severity == compiler::DiagnosticSeverity::Error,
            "diagnostic severity is Error");
  }

  // 7. Division by zero in gate parameter and let declaration
  {
    const std::string divZero =
        "qubit q0\n"
        "RZ(1.0 / 0.0) q0\n";

    compiler::CompilationResult result = compiler::compileProgram(divZero);
    require(!result.success, "division by zero rejected");
    require(!result.diagnostics.empty(), "division by zero diagnostic emitted");
    require(result.diagnostics[0].line == 2, "division by zero localized to line 2");
    require(result.diagnostics[0].message.find("division by zero") != std::string::npos,
            "diagnostic message mentions division by zero");
  }
}

} // namespace test
