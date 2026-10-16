#include "compiler/lexer.h"

#include <cctype>
#include <stdexcept>

Token classifyWord(const std::string &word, int line, int column) {
  if (word == "qubit") {
    return {TokenType::KeywordQubit, word, line, column};
  } else if (word == "measure") {
    return {TokenType::KeywordMeasure, word, line, column};
  } else if (word == "let") {
    return {TokenType::KeywordLet, word, line, column};
  } else if (word == "H" || word == "X" || word == "Z" || word == "CNOT" ||
             word == "RZ") {
    return {TokenType::Gate, word, line, column};
  } else {
    return {TokenType::Identifier, word, line, column};
  }
}

Lexer::Lexer(const std::string &sourceText) { source = sourceText; }

std::vector<Token> Lexer::tokenize() const {
  std::vector<Token> tokens;
  int line = 1;
  int column = 1;
  int position = 0;

  while (position < static_cast<int>(source.size())) {
    char current = source[position];

    if (current == ' ' || current == '\t' || current == '\r') {
      position = position + 1;
      column = column + 1;
      // Newlines are tokens because they end statements in this DSL.
    } else if (current == '\n') {
      tokens.push_back({TokenType::Newline, "\\n", line, column});
      line = line + 1;
      column = 1;
      position = position + 1;
    } else if (current == ',') {
      tokens.push_back({TokenType::Comma, ",", line, column});
      position = position + 1;
      column = column + 1;
    } else if (std::isdigit(static_cast<unsigned char>(current))) {
      int start = position;
      position = position + 1;

      while (position < static_cast<int>(source.size()) &&
             std::isdigit(static_cast<unsigned char>(source[position]))) {
        position = position + 1;
      }

      if (position < static_cast<int>(source.size()) &&
          source[position] == '.') {
        position = position + 1;

        while (position < static_cast<int>(source.size()) &&
               std::isdigit(static_cast<unsigned char>(source[position]))) {
          position = position + 1;
        }
      }

      std::string number = source.substr(start, position - start);
      tokens.push_back({TokenType::Number, number, line, column});
      column = column + (position - start);
    } else if (current == '(') {
      tokens.push_back({TokenType::LeftParen, "(", line, column});
      position = position + 1;
      column = column + 1;
    } else if (current == ')') {
      tokens.push_back({TokenType::RightParen, ")", line, column});
      position = position + 1;
      column = column + 1;
    } else if (current == '+') {
      tokens.push_back({TokenType::Plus, "+", line, column});
      position = position + 1;
      column = column + 1;
    } else if (current == '-') {
      tokens.push_back({TokenType::Minus, "-", line, column});
      position = position + 1;
      column = column + 1;
    } else if (current == '*') {
      tokens.push_back({TokenType::Star, "*", line, column});
      position = position + 1;
      column = column + 1;
    } else if (current == '/') {
      tokens.push_back({TokenType::Slash, "/", line, column});
      position = position + 1;
      column = column + 1;
    } else if (current == '=') {
      tokens.push_back({TokenType::Equal, "=", line, column});
      position = position + 1;
      column = column + 1;
    } else if (std::isalpha(static_cast<unsigned char>(current)) ||
               current == '_') {
      int start = position;
      position = position + 1;

      while (position < static_cast<int>(source.size())) {
        char next = source[position];
        bool isWordCharacter =
            std::isalnum(static_cast<unsigned char>(next)) || next == '_';

        if (!isWordCharacter) {
          break;
        }
        position = position + 1;
      }

      std::string word = source.substr(start, position - start);
      tokens.push_back(classifyWord(word, line, column));
      column = column + (position - start);
    } else {
      throw std::runtime_error("Error: unexpected character '" +
                               std::string(1, current) + "' at line " +
                               std::to_string(line));
    }
  }

  tokens.push_back({TokenType::EndOfFile, "", line, column});
  return tokens;
}

const char *tokenTypeName(TokenType type) {
  switch (type) {
  case TokenType::KeywordQubit:
    return "KeywordQubit";
  case TokenType::KeywordMeasure:
    return "KeywordMeasure";
  case TokenType::KeywordLet:
    return "KeywordLet";
  case TokenType::Gate:
    return "Gate";
  case TokenType::Identifier:
    return "Identifier";
  case TokenType::Number:
    return "Number";
  case TokenType::Comma:
    return "Comma";
  case TokenType::Newline:
    return "Newline";
  case TokenType::EndOfFile:
    return "EndOfFile";
  case TokenType::LeftParen:
    return "LeftParen";
  case TokenType::RightParen:
    return "RightParen";
  case TokenType::Plus:
    return "Plus";
  case TokenType::Minus:
    return "Minus";
  case TokenType::Star:
    return "Star";
  case TokenType::Slash:
    return "Slash";
  case TokenType::Equal:
    return "Equal";
  }

  return "Unknown";
}
