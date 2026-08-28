#include "lexer.h"
#include "parser.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

void printProgram(const Program& program) {
    for (int i = 0; i < static_cast<int>(program.statements.size()); i = i + 1) {
        const Statement& statement = program.statements[i];

        if (statement.type == StatementType::QubitDeclaration) {
            std::cout << "QubitDecl: ";

            for (int j = 0; j < static_cast<int>(statement.names.size()); j = j + 1) {
                if (j != 0) {
                    std::cout << ", ";
                }
                std::cout << statement.names[j];
            }

            std::cout << '\n';
        } else if (statement.type == StatementType::GateOperation) {
            std::cout << "GateOp: " << statement.gate << " [";

            for (int j = 0; j < static_cast<int>(statement.qubits.size()); j = j + 1) {
                if (j != 0) {
                    std::cout << ", ";
                }
                std::cout << statement.qubits[j];
            }

            std::cout << "]\n";
        }
    }
}

std::string readSourceFile(const std::string& filename) {
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        throw std::runtime_error("Could not open source file: " + filename);
    }

    std::string source;
    std::string line;

    while (std::getline(inputFile, line)) {
        // getline removes newlines, so restore them for the lexer.
        source = source + line + "\n";
    }

    return source;
}

int countSourceLines(const std::vector<Token>& tokens) {
    int largestLine = 0;

    for (int i = 0; i < static_cast<int>(tokens.size()); i = i + 1) {
        if (tokens[i].type != TokenType::EndOfFile && tokens[i].line > largestLine) {
            largestLine = tokens[i].line;
        }
    }

    return largestLine;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./qcompile <path-to-program.txt>\n";
        return 1;
    }

    std::string filename = argv[1];

    try {
        std::string source = readSourceFile(filename);

        Lexer lexer(source);
        std::vector<Token> tokens = lexer.tokenize();

        std::cout << "Tokens:\n";
        for (int i = 0; i < static_cast<int>(tokens.size()); i = i + 1) {
            const Token& token = tokens[i];
            std::cout << "  " << tokenTypeName(token.type) << "(\""
                      << token.text << "\", line " << token.line << ")\n";
        }

        Parser parser(tokens);
        Program program = parser.parse();

        std::cout << "\nAST:\n";
        printProgram(program);
        std::cout << "Parsed " << program.statements.size()
                  << " statements from " << countSourceLines(tokens) << " lines.\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        std::cerr << "Failed to parse - see error above.\n";
        return 1;
    }

    return 0;
}
