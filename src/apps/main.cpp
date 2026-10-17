#include "compiler/code_generator.h"
#include "compiler/lexer.h"
#include "compiler/parser.h"
#include "compiler/semantic_checker.h"
#include "simulator/simulator.h"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>

void printExpression(const Expression &expression, int parentPrecedence, bool isRightChild = false) {
  printExpression(std::cout, expression, parentPrecedence, isRightChild);
}

void printProgram(const Program &program) {
  for (int i = 0; i < static_cast<int>(program.statements.size()); i = i + 1) {
    const Statement &statement = program.statements[i];

    if (statement.type == StatementType::QubitDeclaration) {
      std::cout << "QubitDecl: ";

      for (int j = 0; j < static_cast<int>(statement.names.size()); j = j + 1) {
        if (j != 0) {
          std::cout << ", ";
        }
        std::cout << statement.names[j];
      }

      std::cout << '\n';
    } else if (statement.type == StatementType::VariableDeclaration) {
      std::cout << "VariableDecl: " << statement.variableName << " = ";
      printExpression(statement.variableExpression, 0);
      std::cout << '\n';
    } else if (statement.type == StatementType::GateOperation) {
      std::cout << "GateOp: " << statement.gate;

      if (statement.hasParameter) {
        std::cout << "(";
        printExpression(statement.parameter, 0);
        std::cout << ")";
      }

      std::cout << " [";

      for (int j = 0; j < static_cast<int>(statement.qubits.size());
           j = j + 1) {
        if (j != 0) {
          std::cout << ", ";
        }
        std::cout << statement.qubits[j];
      }

      std::cout << "]\n";
    }
  }
}

std::string readSourceFile(const std::string &filename) {
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

int countSourceLines(const std::vector<Token> &tokens) {
  int largestLine = 0;

  for (int i = 0; i < static_cast<int>(tokens.size()); i = i + 1) {
    if (tokens[i].type != TokenType::EndOfFile &&
        tokens[i].line > largestLine) {
      largestLine = tokens[i].line;
    }
  }

  return largestLine;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: ./qcompile <path-to-program.txt>\n";
    return 1;
  }

  std::string arg = argv[1];
  if (arg == "--help" || arg == "-h") {
    std::cout << "Usage: ./qcompile <path-to-program.txt>\n";
    return 0;
  }

  std::string filename = arg;
  try {
    std::string source = readSourceFile(filename);

    Lexer lexer(source);
    std::vector<Token> tokens = lexer.tokenize();

    std::cout << "Tokens:\n";
    for (int i = 0; i < static_cast<int>(tokens.size()); i = i + 1) {
      const Token &token = tokens[i];
      std::cout << "  " << tokenTypeName(token.type) << "(\"" << token.text
                << "\", line " << token.line << ")\n";
    }

    Parser parser(tokens);
    Program program = parser.parse();

    std::cout << "\nAST:\n";
    printProgram(program);
    std::cout << "Parsed " << program.statements.size() << " statements from "
              << countSourceLines(tokens) << " lines.\n";

    SemanticResult semanticResult = checkSemantics(program);
    if (!semanticResult.errors.empty()) {
      std::cerr << "\nSemantic errors:\n";
      for (int i = 0; i < static_cast<int>(semanticResult.errors.size()); i = i + 1) {
        std::cerr << semanticResult.errors[i] << '\n';
      }
      return 1;
    }

    const std::size_t qubitCount = semanticResult.qubitNames.size();
    const std::vector<Instruction> originalInstructions =
        generateInstructions(program, semanticResult);
    const std::vector<Instruction> instructions =
        optimizeInstructions(originalInstructions, qubitCount);

    const auto countInstructions = [](const std::vector<Instruction> &items,
                                      InstructionType type) {
      std::size_t count = 0;
      for (const Instruction &instruction : items) {
        if (instruction.type == type) {
          ++count;
        }
      }
      return count;
    };
    std::cout << "\nInstruction counts:\n"
              << "  Original unitary gates: "
              << countInstructions(originalInstructions,
                                   InstructionType::ApplyGate)
              << "\n  Optimized unitary gates: "
              << countInstructions(instructions, InstructionType::ApplyGate)
              << "\n  Original measurements: "
              << countInstructions(originalInstructions, InstructionType::Measure)
              << "\n  Optimized measurements: "
              << countInstructions(instructions, InstructionType::Measure)
              << '\n';
    printInstructions(instructions);

    if (qubitCount > MAX_SIMULATOR_QUBITS) {
      throw std::invalid_argument(
          "Simulator supports between 0 and " +
          std::to_string(MAX_SIMULATOR_QUBITS) + " qubits"
      );
    }

    Simulator simulator(static_cast<int>(qubitCount));
    simulator.run(instructions);
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    std::cerr << "Failed to parse - see error above.\n";
    return 1;
  }

  return 0;
}
