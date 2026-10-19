#include "compiler/parser.h"

#include <utility>

namespace {
SourceSpan tokenSpan(const Token &token) {
  return token.sourceSpan.isValid() ? token.sourceSpan : SourceSpan{};
}
} // namespace

Parser::Parser(std::vector<Token> tokenList)
    : tokens(std::move(tokenList)),
      eofToken{TokenType::EndOfFile, "", 0, 0, {}},
      current(0),
      recursionDepth(0) {
  if (tokens.empty()) return;
  if (tokens.back().type == TokenType::EndOfFile) {
    eofToken = tokens.back();
    return;
  }

  const SourceSpan lastSpan = tokenSpan(tokens.back());
  if (lastSpan.isValid()) {
    const SourcePosition end = lastSpan.end;
    eofToken = {TokenType::EndOfFile, "", end.line, end.column,
                pointSpan(end)};
  }
}

const Token &Parser::peek() const {
  if (current >= static_cast<int>(tokens.size())) return eofToken;
  return tokens[static_cast<std::size_t>(current)];
}

const Token &Parser::advance() {
  if (current >= static_cast<int>(tokens.size())) return eofToken;
  const Token &token = tokens[static_cast<std::size_t>(current)];
  ++current;
  return token;
}

bool Parser::check(TokenType type) const { return peek().type == type; }

bool Parser::match(TokenType type) {
  if (!check(type)) return false;
  advance();
  return true;
}

void Parser::skipNewlines() {
  while (match(TokenType::Newline)) {
  }
}

void Parser::expect(TokenType type) {
  if (!check(type)) {
    const Token &unexpected = peek();
    const SourceSpan location = tokenSpan(unexpected);
    const SourceSpan span =
        unexpected.type == TokenType::Newline ||
                unexpected.type == TokenType::EndOfFile
            ? pointSpan(location.start)
            : location;
    throw SourceLocatedError(
        "Error: unexpected token '" + unexpected.text + "' at line " +
            std::to_string(unexpected.line),
        span);
  }
  advance();
}

std::vector<std::string> Parser::parseNameList(
    std::vector<SourceSpan> *nameSpans) {
  std::vector<std::string> names;

  expect(TokenType::Identifier);
  names.push_back(tokens[static_cast<std::size_t>(current - 1)].text);
  if (nameSpans) {
    nameSpans->push_back(tokenSpan(tokens[static_cast<std::size_t>(current - 1)]));
  }

  while (match(TokenType::Comma)) {
    expect(TokenType::Identifier);
    names.push_back(tokens[static_cast<std::size_t>(current - 1)].text);
    if (nameSpans) {
      nameSpans->push_back(tokenSpan(tokens[static_cast<std::size_t>(current - 1)]));
    }
  }

  return names;
}

Expression Parser::parseExpression() {
  Expression left = parseTerm();

  while (check(TokenType::Plus) || check(TokenType::Minus)) {
    const Token &operatorToken = peek();
    const char operation = operatorToken.text[0];
    advance();

    Expression right = parseTerm();
    Expression result;
    result.type = ExpressionType::BinaryOperation;
    result.operation = operation;
    result.operatorSpan = tokenSpan(operatorToken);
    result.sourceSpan = joinedSpan(left.sourceSpan, right.sourceSpan);
    result.children.push_back(std::move(left));
    result.children.push_back(std::move(right));
    left = std::move(result);
  }

  return left;
}

Expression Parser::parseTerm() {
  Expression left = parseFactor();

  while (check(TokenType::Star) || check(TokenType::Slash)) {
    const Token &operatorToken = peek();
    const char operation = operatorToken.text[0];
    advance();

    Expression right = parseFactor();
    Expression result;
    result.type = ExpressionType::BinaryOperation;
    result.operation = operation;
    result.operatorSpan = tokenSpan(operatorToken);
    result.sourceSpan = joinedSpan(left.sourceSpan, right.sourceSpan);
    result.children.push_back(std::move(left));
    result.children.push_back(std::move(right));
    left = std::move(result);
  }

  return left;
}

Expression Parser::parseFactor() {
  if (check(TokenType::Number)) {
    const Token &token = advance();
    Expression expression;
    expression.type = ExpressionType::NumberLiteral;
    expression.text = token.text;
    expression.sourceSpan = tokenSpan(token);
    return expression;
  }

  if (check(TokenType::Identifier)) {
    const Token &token = advance();
    Expression expression;
    expression.type = ExpressionType::Identifier;
    expression.text = token.text;
    expression.sourceSpan = tokenSpan(token);
    return expression;
  }

  if (check(TokenType::Minus)) {
    const Token &operatorToken = advance();
    const SourceSpan operatorSpan = tokenSpan(operatorToken);
    if (recursionDepth >= kMaxNestingDepth) {
      throw SourceLocatedError(
          "Error: expression nesting exceeds maximum depth limit of " +
              std::to_string(kMaxNestingDepth) + " at line " +
              std::to_string(operatorToken.line),
          operatorSpan);
    }
    ++recursionDepth;
    struct DepthGuard {
      int &depth;
      ~DepthGuard() { --depth; }
    } guard{recursionDepth};

    Expression expression;
    expression.type = ExpressionType::UnaryOperation;
    expression.operation = '-';
    expression.operatorSpan = operatorSpan;
    expression.children.push_back(parseFactor());
    expression.sourceSpan =
        joinedSpan(operatorSpan, expression.children.front().sourceSpan);
    return expression;
  }

  if (check(TokenType::LeftParen)) {
    const Token &openToken = advance();
    const SourceSpan openSpan = tokenSpan(openToken);
    if (recursionDepth >= kMaxNestingDepth) {
      throw SourceLocatedError(
          "Error: expression nesting exceeds maximum depth limit of " +
              std::to_string(kMaxNestingDepth) + " at line " +
              std::to_string(openToken.line),
          openSpan);
    }
    ++recursionDepth;
    struct DepthGuard {
      int &depth;
      ~DepthGuard() { --depth; }
    } guard{recursionDepth};

    Expression expression = parseExpression();
    expect(TokenType::RightParen);
    expression.sourceSpan = joinedSpan(
        openSpan, tokenSpan(tokens[static_cast<std::size_t>(current - 1)]));
    return expression;
  }

  const Token &unexpected = peek();
  const SourceSpan location = tokenSpan(unexpected);
  throw SourceLocatedError(
      "Error: unexpected token '" + unexpected.text + "' at line " +
          std::to_string(unexpected.line),
      unexpected.type == TokenType::Newline ||
              unexpected.type == TokenType::EndOfFile
          ? pointSpan(location.start)
          : location);
}

Statement Parser::parseStatement() {
  if (match(TokenType::KeywordQubit)) {
    const Token &keyword = tokens[static_cast<std::size_t>(current - 1)];
    Statement statement;
    statement.type = StatementType::QubitDeclaration;
    statement.line = keyword.line;
    statement.hasParameter = false;
    statement.names = parseNameList(&statement.nameSpans);
    statement.sourceSpan =
        joinedSpan(tokenSpan(keyword), tokenSpan(tokens[static_cast<std::size_t>(current - 1)]));
    return statement;
  }

  if (match(TokenType::KeywordLet)) {
    const Token &keyword = tokens[static_cast<std::size_t>(current - 1)];
    Statement statement;
    statement.type = StatementType::VariableDeclaration;
    statement.line = keyword.line;
    statement.hasParameter = false;
    expect(TokenType::Identifier);
    const Token &variableName = tokens[static_cast<std::size_t>(current - 1)];
    statement.variableName = variableName.text;
    statement.variableNameSpan = tokenSpan(variableName);
    expect(TokenType::Equal);
    statement.variableExpression = parseExpression();
    statement.sourceSpan =
        joinedSpan(tokenSpan(keyword), tokenSpan(tokens[static_cast<std::size_t>(current - 1)]));
    return statement;
  }

  if (check(TokenType::Gate) || check(TokenType::KeywordMeasure)) {
    const Token &gateToken = peek();
    const SourceSpan gateSpan = tokenSpan(gateToken);
    const int statementLine = gateToken.line;
    const std::string gateName = gateToken.text;
    advance();

    Statement statement;
    statement.type = StatementType::GateOperation;
    statement.line = statementLine;
    statement.gate = gateName;
    statement.gateSpan = gateSpan;

    if (match(TokenType::LeftParen)) {
      statement.hasParameter = true;
      statement.parameter = parseExpression();
      expect(TokenType::RightParen);
    } else {
      statement.hasParameter = false;
    }

    statement.qubits = parseNameList(&statement.qubitSpans);
    statement.sourceSpan =
        joinedSpan(gateSpan, tokenSpan(tokens[static_cast<std::size_t>(current - 1)]));
    return statement;
  }

  const Token &unexpected = peek();
  throw SourceLocatedError(
      "Error: unexpected token '" + unexpected.text + "' at line " +
          std::to_string(unexpected.line),
      tokenSpan(unexpected));
}

Program Parser::parse() {
  if (tokens.empty() || tokens.back().type != TokenType::EndOfFile) {
    throw SourceLocatedError(
        "Error: malformed token stream missing EndOfFile",
        tokenSpan(eofToken));
  }

  Program program;
  skipNewlines();

  while (!check(TokenType::EndOfFile)) {
    program.statements.push_back(parseStatement());
    // Allow the final statement to end at EOF without a trailing newline.
    if (!check(TokenType::EndOfFile)) expect(TokenType::Newline);
    skipNewlines();
  }

  return program;
}
