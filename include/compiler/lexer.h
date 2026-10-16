#pragma once

#include <string>
#include <vector>

enum class TokenType {
  KeywordQubit,
  KeywordMeasure,
  KeywordLet,
  Gate,
  Identifier,
  Number,
  Comma,
  Newline,
  EndOfFile,
  LeftParen,
  RightParen,
  Plus,
  Minus,
  Star,
  Slash,
  Equal
};

struct Token {
  TokenType type;
  std::string text;
  int line;
  int column;
};

class Lexer {
public:
  Lexer(const std::string &sourceText);

  std::vector<Token> tokenize() const;

private:
  std::string source;
};

const char *tokenTypeName(TokenType type);
