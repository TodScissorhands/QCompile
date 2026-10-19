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
  // 8. Exact UTF-8-byte source spans use 1-based line/column and half-open ends.
  {
    const std::string lexicalSource = "qubit q0\nH q0 @\n";
    const auto lexical = compiler::compileProgram(lexicalSource);
    const auto &lexicalDiagnostic = lexical.diagnostics.at(0);
    require(lexicalDiagnostic.line == 2 && lexicalDiagnostic.column == 6,
            "lexical diagnostic begins at line 2, column 6");
    require(lexicalDiagnostic.sourceSpan.start.offset == 14 &&
                lexicalDiagnostic.sourceSpan.end.offset == 15 &&
                lexicalDiagnostic.sourceSpan.end.line == 2 &&
                lexicalDiagnostic.sourceSpan.end.column == 7,
            "lexical diagnostic selects exactly the invalid character");

    const std::string syntaxSource = "qubit q0\nH q0 q1\n";
    const auto syntax = compiler::compileProgram(syntaxSource);
    const auto &syntaxDiagnostic = syntax.diagnostics.at(0);
    require(syntaxDiagnostic.line == 2 && syntaxDiagnostic.column == 6 &&
                syntaxDiagnostic.sourceSpan.start.offset == 14 &&
                syntaxDiagnostic.sourceSpan.end.offset == 16 &&
                syntaxDiagnostic.sourceSpan.end.column == 8,
            "parser diagnostic selects exactly the unexpected q1 token");

    const auto eof = compiler::compileProgram("qubit q0\nH");
    const auto &eofDiagnostic = eof.diagnostics.at(0);
    require(eofDiagnostic.line == 2 && eofDiagnostic.column == 2 &&
                eofDiagnostic.sourceSpan.start.offset == 10 &&
                eofDiagnostic.sourceSpan.end.offset == 10 &&
                eofDiagnostic.sourceSpan.end.column == 2,
            "missing operand at EOF is a zero-length insertion point");

    const std::string semanticSource =
        "qubit q0\nmeasure q0, q0\nH missing\n";
    const auto semantic = compiler::compileProgram(semanticSource);
    require(semantic.diagnostics.size() == 2,
            "same-line duplicate and later undeclared-qubit diagnostics are both retained");
    const auto &duplicate = semantic.diagnostics[0];
    require(duplicate.line == 2 && duplicate.column == 13 &&
                duplicate.sourceSpan.start.offset == 21 &&
                duplicate.sourceSpan.end.offset == 23 &&
                duplicate.sourceSpan.end.column == 15,
            "duplicate-measure diagnostic selects the second q0 on its line");
    const auto &undeclared = semantic.diagnostics[1];
    require(undeclared.line == 3 && undeclared.column == 3 &&
                undeclared.sourceSpan.start.offset == 26 &&
                undeclared.sourceSpan.end.offset == 33 &&
                undeclared.sourceSpan.end.column == 10,
            "undeclared-qubit diagnostic selects the offending identifier");

    const auto sameLine =
        compiler::compileProgram("qubit q0\nmeasure missing, missing\n");
    require(sameLine.diagnostics.size() == 3,
            "multiple semantic diagnostics on one line are preserved");
    require(sameLine.diagnostics[0].line == 2 &&
                sameLine.diagnostics[0].column == 18 &&
                sameLine.diagnostics[0].sourceSpan.start.offset == 26 &&
                sameLine.diagnostics[0].sourceSpan.end.offset == 33 &&
                sameLine.diagnostics[1].line == 2 &&
                sameLine.diagnostics[1].column == 9 &&
                sameLine.diagnostics[1].sourceSpan.start.offset == 17 &&
                sameLine.diagnostics[1].sourceSpan.end.offset == 24 &&
                sameLine.diagnostics[2].line == 2 &&
                sameLine.diagnostics[2].column == 18 &&
                sameLine.diagnostics[2].sourceSpan.start.offset == 26 &&
                sameLine.diagnostics[2].sourceSpan.end.offset == 33,
            "same-line semantic diagnostics retain exact per-token spans");

    const auto division =
        compiler::compileProgram("qubit q0\nlet theta = 1 / 0\n");
    const auto &divisionDiagnostic = division.diagnostics.at(0);
    require(divisionDiagnostic.line == 2 && divisionDiagnostic.column == 15 &&
                divisionDiagnostic.sourceSpan.start.offset == 23 &&
                divisionDiagnostic.sourceSpan.end.offset == 24 &&
                divisionDiagnostic.sourceSpan.end.column == 16,
            "division-by-zero diagnostic selects the slash operator");

    const std::string utf8Source = "#π\nH q0 @\n";
    const auto utf8 = compiler::compileProgram(utf8Source);
    const auto &utf8Diagnostic = utf8.diagnostics.at(0);
    require(utf8Diagnostic.line == 2 && utf8Diagnostic.column == 6 &&
                utf8Diagnostic.sourceSpan.start.offset == 9 &&
                utf8Diagnostic.sourceSpan.end.offset == 10,
            "source offsets count UTF-8 bytes while columns stay one-based");
    const auto utf8EofTokens = Lexer("#π").tokenize();
    const auto &utf8Eof = utf8EofTokens.back();
    require(utf8Eof.line == 1 && utf8Eof.column == 4 &&
                utf8Eof.sourceSpan.start.offset == 3,
            "UTF-8 byte columns and offsets count both bytes of a multibyte character");
    const std::string supplementaryCharacter = "\xF0\x9F\x98\x80";
    const auto supplementary = compiler::compileProgram(
        std::string("qubit q0\nH q0 ") + supplementaryCharacter + "\n");
    const auto &supplementaryDiagnostic = supplementary.diagnostics.at(0);
    require(supplementaryDiagnostic.line == 2 &&
                supplementaryDiagnostic.column == 6 &&
                supplementaryDiagnostic.sourceSpan.start.offset == 14 &&
                supplementaryDiagnostic.sourceSpan.end.offset == 18 &&
                supplementaryDiagnostic.sourceSpan.end.column == 10 &&
                supplementaryDiagnostic.message.find(supplementaryCharacter) !=
                    std::string::npos,
            "unsupported supplementary code point has a complete UTF-8 byte span");

    std::string malformedUtf8Source = "qubit q0\nH q0 ";
    malformedUtf8Source.push_back(static_cast<char>(0xF0));
    malformedUtf8Source.push_back(' ');
    const auto malformedUtf8 =
        compiler::compileProgram(malformedUtf8Source);
    const auto &malformedDiagnostic = malformedUtf8.diagnostics.at(0);
    require(malformedDiagnostic.sourceSpan.start.offset == 14 &&
                malformedDiagnostic.sourceSpan.end.offset == 15 &&
                malformedDiagnostic.line == 2 &&
                malformedDiagnostic.column == 6,
            "malformed UTF-8 preserves the existing first-byte error location");
  }

  // 9. Topology failures are global and remain explicitly unlocated.
  {
    compiler::CompilerOptions options;
    options.targetTopology = compiler::CouplingMap(2);
    const auto result = compiler::compileProgram(
        "qubit q0, q1\nCNOT q0, q1\n", options);
    require(!result.success && result.diagnostics.size() == 1,
            "disconnected topology produces one routing diagnostic");
    const auto &diagnostic = result.diagnostics[0];
    require(diagnostic.line == 0 && diagnostic.column == 0 &&
                !diagnostic.sourceSpan.isValid(),
            "global routing diagnostic remains unlocated");
  }
}

} // namespace test
