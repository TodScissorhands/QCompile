#include "test_support.h"
#include "compiler/lexer.h"
#include "compiler/parser.h"

namespace test {

void testLexer() {
  // Phase 3: '#' line comment regression tests
  // Standalone comment line
  Program standaloneProg = parse("# this is a standalone comment\nqubit q0\n");
  require(standaloneProg.statements.size() == 1, "standalone comment ignored");

  // Trailing comment after valid statement
  Program trailingProg = parse("qubit q0 # trailing comment\nH q0\n");
  require(trailingProg.statements.size() == 2, "trailing comment ignored");

  // Comment containing text that would otherwise be invalid syntax
  Program invalidSyntaxInComment = parse("qubit q0\n# @$%^! &*() ~` let = = 123..\nH q0\n");
  require(invalidSyntaxInComment.statements.size() == 2, "invalid syntax inside comment ignored");

  // Blank lines surrounding comments
  Program blankLinesProg = parse("\n\n# comment with blank lines\n\nqubit q0\n\n# another comment\n\n");
  require(blankLinesProg.statements.size() == 1, "comments surrounded by blank lines parsed");

  // Comment reaching end-of-file without a terminating newline
  std::vector<Token> eofCommentTokens = Lexer("qubit q0\n# 12345").tokenize();
  require(eofCommentTokens.size() == 4, "tokens include qubit, q0, newline, EOF");
  require(eofCommentTokens[3].type == TokenType::EndOfFile, "EOF token emitted");
  require(eofCommentTokens[3].line == 2, "EOF line is 2");
  require(eofCommentTokens[3].column == 8, "EOF column reflects scanned comment characters");
  Program eofCommentProg = Parser(eofCommentTokens).parse();
  require(eofCommentProg.statements.size() == 1, "comment at EOF without newline parsed");
}

} // namespace test
