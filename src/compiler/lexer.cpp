#include "compiler/lexer.h"

#include <cctype>
#include <utility>

Token classifyWord(const std::string &word, int line, int column) {
  if (word == "qubit") {
    return {TokenType::KeywordQubit, word, line, column, {}};
  } else if (word == "measure") {
    return {TokenType::KeywordMeasure, word, line, column, {}};
  } else if (word == "let") {
    return {TokenType::KeywordLet, word, line, column, {}};
  } else if (word == "H" || word == "X" || word == "Z" || word == "CNOT" ||
             word == "RZ") {
    return {TokenType::Gate, word, line, column, {}};
  } else {
    return {TokenType::Identifier, word, line, column, {}};
  }
}

namespace {
bool isUtf8Continuation(unsigned char byte) {
  return (byte & 0xC0) == 0x80;
}

std::size_t unexpectedCharacterByteLength(const std::string &source,
                                          std::size_t offset) {
  const unsigned char lead = static_cast<unsigned char>(source[offset]);
  std::size_t length = 1;
  if (lead >= 0xC2 && lead <= 0xDF) {
    length = 2;
  } else if (lead >= 0xE0 && lead <= 0xEF) {
    length = 3;
  } else if (lead >= 0xF0 && lead <= 0xF4) {
    length = 4;
  } else {
    return 1;
  }

  if (length > source.size() - offset) return 1;
  const unsigned char second =
      static_cast<unsigned char>(source[offset + 1]);
  if (!isUtf8Continuation(second) ||
      (lead == 0xE0 && second < 0xA0) ||
      (lead == 0xED && second > 0x9F) ||
      (lead == 0xF0 && second < 0x90) ||
      (lead == 0xF4 && second > 0x8F)) {
    return 1;
  }
  for (std::size_t i = 2; i < length; ++i) {
    if (!isUtf8Continuation(
            static_cast<unsigned char>(source[offset + i]))) {
      return 1;
    }
  }
  return length;
}

Token makeToken(TokenType type, std::string text, SourcePosition start,
                SourcePosition end) {
  return {type, std::move(text), start.line, start.column, {start, end}};
}
} // namespace

Lexer::Lexer(const std::string &sourceText) { source = sourceText; }

std::vector<Token> Lexer::tokenize() const {
  std::vector<Token> tokens;
  int line = 1;
  int column = 1;
  int position = 0;

  while (position < static_cast<int>(source.size())) {
    const char current = source[static_cast<std::size_t>(position)];

    if (current == ' ' || current == '\t' || current == '\r') {
      ++position;
      ++column;
    } else if (current == '#') {
      while (position < static_cast<int>(source.size()) &&
             source[static_cast<std::size_t>(position)] != '\n') {
        ++position;
        ++column;
      }
      // Line remains unchanged; column is reset on '\n' or retained at EOF.
      // Do not consume '\n' here so statement-terminating newline tokens are preserved.
      // Newlines are tokens because they end statements in this DSL.
    } else if (current == '\n') {
      const SourcePosition start{static_cast<std::size_t>(position), line, column};
      tokens.push_back(makeToken(TokenType::Newline, "\\n", start,
                                 {start.offset + 1, line + 1, 1}));
      ++line;
      column = 1;
      ++position;
    } else if (current == ',') {
      const SourcePosition start{static_cast<std::size_t>(position), line, column};
      tokens.push_back(makeToken(TokenType::Comma, ",", start,
                                 {start.offset + 1, line, column + 1}));
      ++position;
      ++column;
    } else if (std::isdigit(static_cast<unsigned char>(current))) {
      const int startOffset = position;
      const int startColumn = column;
      ++position;

      while (position < static_cast<int>(source.size()) &&
             std::isdigit(static_cast<unsigned char>(source[static_cast<std::size_t>(position)]))) {
        ++position;
      }
      if (position < static_cast<int>(source.size()) &&
          source[static_cast<std::size_t>(position)] == '.') {
        ++position;
        while (position < static_cast<int>(source.size()) &&
               std::isdigit(static_cast<unsigned char>(source[static_cast<std::size_t>(position)]))) {
          ++position;
        }
      }

      const std::string number = source.substr(
          static_cast<std::size_t>(startOffset),
          static_cast<std::size_t>(position - startOffset));
      const SourcePosition start{static_cast<std::size_t>(startOffset), line, startColumn};
      const SourcePosition end{static_cast<std::size_t>(position), line,
                               startColumn + position - startOffset};
      tokens.push_back(makeToken(TokenType::Number, number, start, end));
      column = end.column;
    } else if (current == '(' || current == ')' || current == '+' ||
               current == '-' || current == '*' || current == '/' ||
               current == '=') {
      TokenType type = TokenType::LeftParen;
      switch (current) {
      case ')': type = TokenType::RightParen; break;
      case '+': type = TokenType::Plus; break;
      case '-': type = TokenType::Minus; break;
      case '*': type = TokenType::Star; break;
      case '/': type = TokenType::Slash; break;
      case '=': type = TokenType::Equal; break;
      default: break;
      }
      const SourcePosition start{static_cast<std::size_t>(position), line, column};
      tokens.push_back(makeToken(type, std::string(1, current), start,
                                 {start.offset + 1, line, column + 1}));
      ++position;
      ++column;
    } else if (std::isalpha(static_cast<unsigned char>(current)) ||
               current == '_') {
      const int startOffset = position;
      const int startColumn = column;
      ++position;
      while (position < static_cast<int>(source.size())) {
        const char next = source[static_cast<std::size_t>(position)];
        const bool isWordCharacter =
            std::isalnum(static_cast<unsigned char>(next)) || next == '_';
        if (!isWordCharacter) break;
        ++position;
      }

      const std::string word = source.substr(
          static_cast<std::size_t>(startOffset),
          static_cast<std::size_t>(position - startOffset));
      Token token = classifyWord(word, line, startColumn);
      token.sourceSpan = {
          {static_cast<std::size_t>(startOffset), line, startColumn},
          {static_cast<std::size_t>(position), line,
           startColumn + position - startOffset}};
      tokens.push_back(std::move(token));
      column = startColumn + position - startOffset;
    } else {
      const std::size_t characterLength = unexpectedCharacterByteLength(
          source, static_cast<std::size_t>(position));
      const SourcePosition start{static_cast<std::size_t>(position), line,
                                 column};
      const SourcePosition end{
          start.offset + characterLength, line,
          column + static_cast<int>(characterLength)};
      throw SourceLocatedError(
          "Error: unexpected character '" +
              source.substr(static_cast<std::size_t>(position),
                            characterLength) +
              "' at line " + std::to_string(line),
          {start, end});
    }
  }

  const SourcePosition eof{source.size(), line, column};
  tokens.push_back(makeToken(TokenType::EndOfFile, "", eof, eof));
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
