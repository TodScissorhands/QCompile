#include "test_support.h"
#include "compiler/lexer.h"
#include "compiler/parser.h"

#include <cmath>
#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace test {

void testParser() {
  std::vector<Token> tokens = Lexer("qubit q0\nRZ(-2 * (3 + 1)) q0\n").tokenize();
  require(tokens.front().type == TokenType::KeywordQubit, "qubit keyword tokenized");
  require(tokens.back().type == TokenType::EndOfFile, "EOF token emitted");
  Program program = Parser(tokens).parse();
  require(program.statements.size() == 2, "two statements parsed");
  const Expression &expression = program.statements[1].parameter;
  require(expression.type == ExpressionType::BinaryOperation &&
              expression.operation == '*',
          "multiplication is outer precedence level");
  require(expression.children[0].type == ExpressionType::UnaryOperation,
          "unary minus parsed as a factor");
  const Expression &parenthesized = expression.children[1];
  require(parenthesized.operation == '+', "parentheses preserve addition grouping");

  Program precedence = parse("qubit q0\nRZ(2 + 3 * 4) q0\n");
  const Expression &sum = precedence.statements[1].parameter;
  require(sum.operation == '+', "addition is expression root");
  require(sum.children[1].operation == '*', "multiplication binds before addition");

  bool malformedRejected = false;
  try {
    parse("qubit q0\nRZ(2 +) q0\n");
  } catch (const std::exception &) {
    malformedRejected = true;
  }
  require(malformedRejected, "malformed expression rejected");
  bool syntaxRejected = false;
  try {
    parse("qubit q0\nH\n");
  } catch (const std::exception &) {
    syntaxRejected = true;
  }
  require(syntaxRejected, "invalid syntax rejected");

  // Equivalence test: Bell state circuit with comments vs without comments
  std::string noCommentsSource =
      "qubit q0, q1\n"
      "H q0\n"
      "CNOT q0, q1\n"
      "measure q0, q1\n";
  std::string withCommentsSource =
      "# Quantum entanglement demo\n"
      "qubit q0, q1 # allocate pair\n"
      "# Create superposition\n"
      "H q0 # Hadamard on q0\n"
      "# Entangle control and target\n"
      "CNOT q0, q1 # Bell pair created\n"
      "# Measure both qubits\n"
      "measure q0, q1 # record outcomes\n"
      "# End of circuit";

  Program progClean;
  SemanticResult semClean = check(noCommentsSource, progClean);
  std::vector<Instruction> instrsClean = generateInstructions(progClean, semClean);

  Program progCommented;
  SemanticResult semCommented = check(withCommentsSource, progCommented);
  std::vector<Instruction> instrsCommented = generateInstructions(progCommented, semCommented);

  require(instrsClean.size() == instrsCommented.size(),
          "instruction counts match between commented and uncommented programs");
  for (std::size_t k = 0; k < instrsClean.size(); ++k) {
    require(instrsClean[k].type == instrsCommented[k].type &&
            instrsClean[k].gateName == instrsCommented[k].gateName &&
            instrsClean[k].qubitIndices == instrsCommented[k].qubitIndices &&
            instrsClean[k].parameterValue == instrsCommented[k].parameterValue,
            "commented program generates identical instructions to uncommented program");
  }

  std::string outClean = run(noCommentsSource);
  std::string outCommented = run(withCommentsSource);

  auto extractProbabilities = [](const std::string &output) -> std::string {
    std::size_t start = output.find("Pre-measurement probabilities:\n");
    if (start == std::string::npos) return "";
    std::size_t end = output.find("Result:", start);
    if (end == std::string::npos) return output.substr(start);
    return output.substr(start, end - start);
  };

  std::string probsClean = extractProbabilities(outClean);
  std::string probsCommented = extractProbabilities(outCommented);
  require(!probsClean.empty(), "clean program produces probabilities");
  require(probsClean == probsCommented,
          "commented Bell-state program produces identical probability distribution to uncommented program");

  // Phase 4: Parser and expression safety regression tests
  // 1. Valid nested expression at moderate depth (50 parentheses)
  std::string nested50 = std::string(50, '(') + "42.0" + std::string(50, ')');
  Program nestedProg = parse("qubit q0\nlet x = " + nested50 + "\n");
  SemanticResult nestedSem = checkSemantics(nestedProg);
  require(nestedSem.errors.empty(), "50-level nested expression accepted");
  require(!nestedSem.variableValues.empty() && std::abs(nestedSem.variableValues[0] - 42.0) < 1e-12,
          "50-level nested expression evaluates correctly");

  // 2. Expression exceeding configured nesting limit rejected cleanly
  std::string nestedExceed = std::string(300, '(') + "1.0" + std::string(300, ')');
  bool depthExceededRejected = false;
  std::string depthErrorMsg;
  try {
    parse("qubit q0\nlet x = " + nestedExceed + "\n");
  } catch (const std::exception &e) {
    depthExceededRejected = true;
    depthErrorMsg = e.what();
  }
  require(depthExceededRejected, "deeply nested expression exceeding limit rejected");
  require(depthErrorMsg.find("expression nesting exceeds maximum depth limit") != std::string::npos,
          "clean diagnostic for depth limit exceeded");

  // Unary operator chain exceeding depth limit
  std::string unaryExceed = std::string(300, '-') + "1.0";
  bool unaryDepthRejected = false;
  try {
    parse("qubit q0\nlet x = " + unaryExceed + "\n");
  } catch (const std::exception &e) {
    unaryDepthRejected = true;
  }
  require(unaryDepthRejected, "unary operator chain exceeding limit rejected");

  // 3. Long left-associative expression (2,000 terms) evaluates correctly
  std::string longExpr2000 = "1.0";
  for (int idx = 1; idx < 2000; ++idx) {
    longExpr2000 += " + 1.0";
  }
  Program longProg2000 = parse("qubit q0\nlet x = " + longExpr2000 + "\n");
  SemanticResult longSem2000 = checkSemantics(longProg2000);
  require(longSem2000.errors.empty(), "2000-term expression accepted");
  require(!longSem2000.variableValues.empty() && std::abs(longSem2000.variableValues[0] - 2000.0) < 1e-9,
          "2000-term expression evaluates to 2000.0");

  // 4. Substantially longer expression (10,000 terms) evaluated safely and efficiently
  std::string longExpr10000 = "1.0";
  for (int idx = 1; idx < 10000; ++idx) {
    longExpr10000 += " + 1.0";
  }
  Program longProg10000 = parse("qubit q0\nlet x = " + longExpr10000 + "\n");
  SemanticResult longSem10000 = checkSemantics(longProg10000);
  require(longSem10000.errors.empty(), "10000-term expression accepted");
  require(!longSem10000.variableValues.empty() && std::abs(longSem10000.variableValues[0] - 10000.0) < 1e-8,
          "10000-term expression evaluates to 10000.0");

  // 5. Explicit associativity and precedence evaluation regression tests
  // Left-associative subtraction: 10 - 3 - 2 = (10 - 3) - 2 = 5
  Program subProg = parse("qubit q0\nlet x = 10 - 3 - 2\n");
  SemanticResult subSem = checkSemantics(subProg);
  require(subSem.errors.empty(), "10 - 3 - 2 parses and checks");
  require(!subSem.variableValues.empty() && std::abs(subSem.variableValues[0] - 5.0) < 1e-12,
          "10 - 3 - 2 evaluates to 5");

  // Left-associative division: 24 / 4 / 2 = (24 / 4) / 2 = 3
  Program divProg = parse("qubit q0\nlet x = 24 / 4 / 2\n");
  SemanticResult divSem = checkSemantics(divProg);
  require(divSem.errors.empty(), "24 / 4 / 2 parses and checks");
  require(!divSem.variableValues.empty() && std::abs(divSem.variableValues[0] - 3.0) < 1e-12,
          "24 / 4 / 2 evaluates to 3");

  // Mixed precedence: 2 + 3 * 4 = 2 + (3 * 4) = 14
  Program mixProg1 = parse("qubit q0\nlet x = 2 + 3 * 4\n");
  SemanticResult mixSem1 = checkSemantics(mixProg1);
  require(mixSem1.errors.empty(), "2 + 3 * 4 parses and checks");
  require(!mixSem1.variableValues.empty() && std::abs(mixSem1.variableValues[0] - 14.0) < 1e-12,
          "2 + 3 * 4 evaluates to 14");

  // Unary and parentheses: -2 * (3 + 1) = -2 * 4 = -8
  Program mixProg2 = parse("qubit q0\nlet x = -2 * (3 + 1)\n");
  SemanticResult mixSem2 = checkSemantics(mixProg2);
  require(mixSem2.errors.empty(), "-2 * (3 + 1) parses and checks");
  require(!mixSem2.variableValues.empty() && std::abs(mixSem2.variableValues[0] - (-8.0)) < 1e-12,
          "-2 * (3 + 1) evaluates to -8");

  // 6. Expression printer grouping and precedence tests
  std::vector<std::pair<std::string, std::string>> printerCases = {
      {"10 - 3 - 2", "10 - 3 - 2"},
      {"10 - (3 - 2)", "10 - (3 - 2)"},
      {"24 / 4 / 2", "24 / 4 / 2"},
      {"24 / (4 / 2)", "24 / (4 / 2)"},
      {"2 * (3 / 2)", "2 * (3 / 2)"},
      {"2 * 3 / 2", "2 * 3 / 2"},
      {"2 + (3 - 1)", "2 + (3 - 1)"},
      {"2 + 3 - 1", "2 + 3 - 1"},
      {"2 + 3 * 4", "2 + 3 * 4"},
      {"(2 + 3) * 4", "(2 + 3) * 4"},
      {"-(2 + 3)", "-(2 + 3)"},
      {"2 * -3", "2 * -3"},
      {"(2 + 3 * 4) - 5 / 2", "2 + 3 * 4 - 5 / 2"},
      {"((2 + 3) * 4) / 5", "(2 + 3) * 4 / 5"}
  };

  for (const auto &pair : printerCases) {
    Program p = parse("qubit q0\nlet x = " + pair.first + "\n");
    std::string printed = formatExpression(p.statements[1].variableExpression);
    require(printed == pair.second,
            "Expression printer formatted '" + pair.first + "' as '" + printed + "', expected '" + pair.second + "'");
  }

  // Phase 6: Parser token ownership, lifetime, and bounds safety regression tests
  // 1. Parsing directly from an rvalue temporary vector
  Program rvalueProg = Parser(Lexer("qubit q0\nH q0\n").tokenize()).parse();
  require(rvalueProg.statements.size() == 2, "rvalue temporary token vector parsed successfully");

  // 2. Transferring/returning parser constructed from temporary token storage
  auto makeParser = [](const std::string &src) {
    return Parser(Lexer(src).tokenize());
  };
  Parser returnedParser = makeParser("qubit q0, q1\nCNOT q0, q1\n");
  Program returnedProg = returnedParser.parse();
  require(returnedProg.statements.size() == 2, "parser returned across stack frames parses valid AST");

  // 3. Empty token vector rejected with diagnostic
  bool emptyRejected = false;
  std::string emptyError;
  try {
    Parser emptyParser(std::vector<Token>{});
    emptyParser.parse();
  } catch (const std::exception &e) {
    emptyRejected = true;
    emptyError = e.what();
  }
  require(emptyRejected, "empty token vector rejected as malformed");
  require(emptyError.find("malformed token stream missing EndOfFile") != std::string::npos,
          "expected diagnostic for empty token stream");

  // 4. Missing EOF token (truncated token stream) rejected with diagnostic
  std::vector<Token> truncatedTokens = {
      {TokenType::KeywordQubit, "qubit", 1, 1, {}},
      {TokenType::Identifier, "q0", 1, 7, {}}
  }; // lacks EndOfFile or Newline
  bool truncatedRejected = false;
  std::string truncatedError;
  try {
    Parser truncatedParser(truncatedTokens);
    truncatedParser.parse();
  } catch (const std::exception &e) {
    truncatedRejected = true;
    truncatedError = e.what();
  }
  require(truncatedRejected, "truncated token stream missing EOF rejected as malformed");
  require(truncatedError.find("malformed token stream missing EndOfFile") != std::string::npos,
          "expected diagnostic for truncated token stream missing EOF");

  // Phase 7: AST field initialization safety regression tests
  // 1. Default-constructed Expression fields have deterministic defaults
  Expression defaultExpr;
  require(defaultExpr.type == ExpressionType::NumberLiteral, "default expression type is NumberLiteral");
  require(defaultExpr.operation == '\0', "default expression operation is null character");
  require(defaultExpr.text.empty(), "default expression text is empty");
  require(defaultExpr.children.empty(), "default expression children vector is empty");

  // 2. Default-constructed Statement fields have deterministic defaults, including nested expressions
  Statement defaultStmt;
  require(defaultStmt.type == StatementType::QubitDeclaration, "default statement type is QubitDeclaration");
  require(defaultStmt.line == 0, "default statement line is 0");
  require(defaultStmt.names.empty(), "default statement names vector is empty");
  require(defaultStmt.variableName.empty(), "default statement variableName is empty");
  require(defaultStmt.gate.empty(), "default statement gate is empty");
  require(defaultStmt.qubits.empty(), "default statement qubits vector is empty");
  require(!defaultStmt.hasParameter, "default statement hasParameter is false");
  require(defaultStmt.variableExpression.type == ExpressionType::NumberLiteral,
          "nested variableExpression type defaults to NumberLiteral");
  require(defaultStmt.variableExpression.operation == '\0',
          "nested variableExpression operation defaults to null character");
  require(defaultStmt.parameter.type == ExpressionType::NumberLiteral,
          "nested parameter type defaults to NumberLiteral");
  require(defaultStmt.parameter.operation == '\0',
          "nested parameter operation defaults to null character");

  // 3. Parser-produced literal expressions have deterministic operation values
  Program literalProg = parse("qubit q0\nlet x = 42.5\nRZ(x) q0\n");
  require(literalProg.statements.size() == 3, "three statements parsed");
  const Expression &letExpr = literalProg.statements[1].variableExpression;
  require(letExpr.type == ExpressionType::NumberLiteral, "parsed let expression is NumberLiteral");
  require(letExpr.operation == '\0', "parsed number literal operation is null character");
  require(letExpr.text == "42.5", "parsed number literal text matches");

  const Expression &paramExpr = literalProg.statements[2].parameter;
  require(paramExpr.type == ExpressionType::Identifier, "parsed gate parameter expression is Identifier");
  require(paramExpr.operation == '\0', "parsed identifier expression operation is null character");
  require(paramExpr.text == "x", "parsed identifier text matches");

  // Phase 9: Empty program parsing and semantic evaluation
  Program emptyValidProg = parse("");
  require(emptyValidProg.statements.empty(), "empty string parses cleanly to 0 statements");
  SemanticResult emptySem = checkSemantics(emptyValidProg);
  require(emptySem.errors.empty(), "empty program passes semantic checking");
  require(emptySem.qubitNames.empty(), "empty program has no declared qubits");
  std::vector<Instruction> emptyInstrs = generateInstructions(emptyValidProg, emptySem);
  require(emptyInstrs.empty(), "empty program generates zero instructions");
  // A truncated token stream reports a structured insertion point after its last token.
  auto missingEofTokens = Lexer("qubit q0").tokenize();
  missingEofTokens.pop_back();
  bool missingEofLocated = false;
  try {
    (void)Parser(std::move(missingEofTokens)).parse();
  } catch (const compiler::SourceLocatedError &error) {
    missingEofLocated = true;
    require(error.span().isValid() && error.span().start.offset == 8 &&
                error.span().end.offset == 8 &&
                error.span().start.line == 1 &&
                error.span().start.column == 9 &&
                error.span().end.column == 9,
            "missing EOF is located at the exact end of the last token");
  }
  require(missingEofLocated, "malformed token stream raises a located parser error");
}

} // namespace test
