#include "parser.h"

#include <stdexcept>

Parser::Parser(const std::vector<Token>& tokenList)
    : tokens(tokenList), current(0) {
}

const Token& Parser::peek() const {
    return tokens[current];
}

const Token& Parser::advance() {
    const Token& token = tokens[current];
    current = current + 1;
    return token;
}

bool Parser::check(TokenType type) const {
    return peek().type == type;
}

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
        throw std::runtime_error(
            "Error: unexpected token '" + peek().text + "' at line " +
            std::to_string(peek().line));
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

Statement Parser::parseStatement() {
    if (match(TokenType::KeywordQubit)) {
        Statement statement;
        statement.type = StatementType::QubitDeclaration;
        statement.names = parseNameList();
        return statement;
    }

    if (check(TokenType::Gate) || check(TokenType::KeywordMeasure)) {
        std::string gateName = peek().text;
        advance();

        Statement statement;
        statement.type = StatementType::GateOperation;
        statement.gate = gateName;
        statement.qubits = parseNameList();
        return statement;
    }

    throw std::runtime_error(
        "Error: unexpected token '" + peek().text + "' at line " +
        std::to_string(peek().line));
}

Program Parser::parse() {
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
