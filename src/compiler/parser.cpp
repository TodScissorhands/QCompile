#include "compiler/parser.h"

#include <stdexcept>
#include <utility>

Parser::Parser(std::vector<Token> tokenList)
    : tokens(std::move(tokenList)),
      eofToken{TokenType::EndOfFile, "", 1, 1},
      current(0),
      recursionDepth(0) {
  if (!tokens.empty()) {
    eofToken.line = tokens.back().line;
    eofToken.column = tokens.back().column;
  }
}

const Token &Parser::peek() const {
  if (current >= static_cast<int>(tokens.size())) {
    return eofToken;
  }
  return tokens[current];
}

const Token &Parser::advance() {
  if (current >= static_cast<int>(tokens.size())) {
    return eofToken;
  }
  const Token &token = tokens[current];
  current = current + 1;
  return token;
}
bool Parser::check(TokenType type) const { return peek().type == type; }

bool Parser::match(TokenType type) {
  if (!check(type)) {
    return false;
  }

  advance();
  return true;
}

void Parser::skipNewlines() {
  while (match(TokenType::Newline)) {
  }
}

void Parser::expect(TokenType type) {
  if (!check(type)) {
    throw std::runtime_error("Error: unexpected token '" + peek().text +
                             "' at line " + std::to_string(peek().line));
  }

  advance();
}

std::vector<std::string> Parser::parseNameList() {
  std::vector<std::string> names;

  expect(TokenType::Identifier);
  names.push_back(tokens[current - 1].text);

  while (match(TokenType::Comma)) {
    expect(TokenType::Identifier);
    names.push_back(tokens[current - 1].text);
  }

  return names;
}

Expression Parser::parseExpression() {
  Expression left = parseTerm();

  while (check(TokenType::Plus) || check(TokenType::Minus)) {
    char operation = peek().text[0];
    advance();

    Expression right = parseTerm();
    Expression result;
    result.type = ExpressionType::BinaryOperation;
    result.operation = operation;
    result.children.push_back(std::move(left));
    result.children.push_back(std::move(right));
    left = std::move(result);
  }

  return left;
}

Expression Parser::parseTerm() {
  Expression left = parseFactor();

  while (check(TokenType::Star) || check(TokenType::Slash)) {
    char operation = peek().text[0];
    advance();

    Expression right = parseFactor();
    Expression result;
    result.type = ExpressionType::BinaryOperation;
    result.operation = operation;
    result.children.push_back(std::move(left));
    result.children.push_back(std::move(right));
    left = std::move(result);
  }

  return left;
}

Expression Parser::parseFactor() {
  if (check(TokenType::Number)) {
    Expression expression;
    expression.type = ExpressionType::NumberLiteral;
    expression.text = advance().text;
    return expression;
  }

  if (check(TokenType::Identifier)) {
    Expression expression;
    expression.type = ExpressionType::Identifier;
    expression.text = advance().text;
    return expression;
  }

  if (match(TokenType::Minus)) {
    int operatorLine = peek().line;
    if (recursionDepth >= kMaxNestingDepth) {
      throw std::runtime_error("Error: expression nesting exceeds maximum depth limit of " +
                               std::to_string(kMaxNestingDepth) + " at line " +
                               std::to_string(operatorLine));
    }
    recursionDepth = recursionDepth + 1;
    struct DepthGuard {
      int &depth;
      ~DepthGuard() { depth = depth - 1; }
    } guard{recursionDepth};

    Expression expression;
    expression.type = ExpressionType::UnaryOperation;
    expression.operation = '-';
    expression.children.push_back(parseFactor());
    return expression;
  }

  if (match(TokenType::LeftParen)) {
    int parenLine = peek().line;
    if (recursionDepth >= kMaxNestingDepth) {
      throw std::runtime_error("Error: expression nesting exceeds maximum depth limit of " +
                               std::to_string(kMaxNestingDepth) + " at line " +
                               std::to_string(parenLine));
    }
    recursionDepth = recursionDepth + 1;
    struct DepthGuard {
      int &depth;
      ~DepthGuard() { depth = depth - 1; }
    } guard{recursionDepth};

    Expression expression = parseExpression();
    expect(TokenType::RightParen);
    return expression;
  }

  throw std::runtime_error("Error: unexpected token '" + peek().text +
                           "' at line " + std::to_string(peek().line));
}

Statement Parser::parseStatement() {
  if (match(TokenType::KeywordQubit)) {
    Statement statement;
    statement.type = StatementType::QubitDeclaration;
    statement.line = tokens[current - 1].line;
    statement.hasParameter = false;
    statement.names = parseNameList();
    return statement;
  }

  if (match(TokenType::KeywordLet)) {
    Statement statement;
    statement.type = StatementType::VariableDeclaration;
    statement.line = tokens[current - 1].line;
    statement.hasParameter = false;
    expect(TokenType::Identifier);
    statement.variableName = tokens[current - 1].text;
    expect(TokenType::Equal);
    statement.variableExpression = parseExpression();
    return statement;
  }

  if (check(TokenType::Gate) || check(TokenType::KeywordMeasure)) {
    int statementLine = peek().line;
    std::string gateName = peek().text;
    advance();

    Statement statement;
    statement.type = StatementType::GateOperation;
    statement.line = statementLine;
    statement.gate = gateName;

    if (match(TokenType::LeftParen)) {
      statement.hasParameter = true;
      statement.parameter = parseExpression();
      expect(TokenType::RightParen);
    } else {
      statement.hasParameter = false;
    }

    statement.qubits = parseNameList();
    return statement;
  }

  throw std::runtime_error("Error: unexpected token '" + peek().text +
                           "' at line " + std::to_string(peek().line));
}

Program Parser::parse() {
  if (tokens.empty() || tokens.back().type != TokenType::EndOfFile) {
    throw std::runtime_error("Error: malformed token stream missing EndOfFile");
  }

  Program program;
  skipNewlines();

  while (!check(TokenType::EndOfFile)) {
    program.statements.push_back(parseStatement());

    // Allow the final statement to end at EOF without a trailing newline.
    if (!check(TokenType::EndOfFile)) {
      expect(TokenType::Newline);
    }
    skipNewlines();
  }

  return program;
}
