#pragma once

#include "compiler/ast.h"
#include "compiler/lexer.h"

#include <string>
#include <vector>

class Parser {
public:
  explicit Parser(std::vector<Token> tokenList);

  Program parse();

private:
  std::vector<Token> tokens;
  Token eofToken;

  int current;
  int recursionDepth;

  static constexpr int kMaxNestingDepth = 256;

  const Token &peek() const;
  const Token &advance();

  bool check(TokenType type) const;

  bool match(TokenType type);

  void skipNewlines();

  void expect(TokenType type);

  std::vector<std::string> parseNameList(
      std::vector<SourceSpan> *nameSpans = nullptr);

  Expression parseExpression();

  Expression parseTerm();

  Expression parseFactor();

  Statement parseStatement();
};
