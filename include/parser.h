#pragma once

#include "ast.h"
#include "lexer.h"

#include <string>
#include <vector>

class Parser {
public:
  Parser(const std::vector<Token> &tokenList);

  Program parse();

private:
  const std::vector<Token> &tokens;

  int current;

  const Token &peek() const;

  const Token &advance();

  bool check(TokenType type) const;

  bool match(TokenType type);

  void skipNewlines();

  void expect(TokenType type);

  std::vector<std::string> parseNameList();

  Expression parseExpression();

  Expression parseTerm();

  Expression parseFactor();

  Statement parseStatement();
};
