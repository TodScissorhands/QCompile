#include "compiler/coupling_map.h"
#include "compiler/router.h"
#include "compiler/circuit_dag.h"
#include "compiler/code_generator.h"
#include "compiler/lexer.h"
#include "compiler/parser.h"
#include "compiler/semantic_checker.h"
#include "simulator/simulator.h"

#include <cstddef>
#include <fstream>
#include <iomanip>
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
    std::cerr << "Usage: ./qcompile [--target <line|ring|star|grid|all>] <path-to-program.txt>\n";
    return 1;
  }

  std::string filename;
  std::string targetTopologyType;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      std::cout << "Usage: ./qcompile [--target <line|ring|star|grid|all>] <path-to-program.txt>\n";
      return 0;
    } else if (arg == "--target" && i + 1 < argc) {
      targetTopologyType = argv[++i];
    } else if (arg.rfind("--target=", 0) == 0) {
      targetTopologyType = arg.substr(9);
    } else if (!arg.empty() && arg[0] == '-') {
      std::cerr << "Unknown option: " << arg << "\n";
      return 1;
    } else {
      filename = arg;
    }
  }

  if (filename.empty()) {
    std::cerr << "Usage: ./qcompile [--target <line|ring|star|grid|all>] <path-to-program.txt>\n";
    return 1;
  }

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

    compiler::CircuitDag rawDag = compiler::CircuitDag::fromInstructions(originalInstructions, qubitCount);
    compiler::CircuitDag optDag = compiler::CircuitDag::fromInstructions(instructions, qubitCount);

    std::cout << "\nInstruction counts:\n"
              << "  Original unitary gates: "
              << rawDag.gateCount(InstructionType::ApplyGate)
              << "\n  Optimized unitary gates: "
              << optDag.gateCount(InstructionType::ApplyGate)
              << "\n  Original two-qubit gates (CNOT): "
              << rawDag.twoQubitGateCount()
              << "\n  Optimized two-qubit gates (CNOT): "
              << optDag.twoQubitGateCount()
              << "\n  Original circuit depth: "
              << rawDag.criticalPathDepth() << " layers"
              << "\n  Optimized circuit depth: "
              << optDag.criticalPathDepth() << " layers"
              << "\n  Original measurements: "
              << rawDag.gateCount(InstructionType::Measure)
              << "\n  Optimized measurements: "
              << optDag.gateCount(InstructionType::Measure)
              << '\n';

    if (!targetTopologyType.empty()) {
      std::size_t physicalQubits = (qubitCount < 2) ? 2 : qubitCount;
      compiler::CouplingMap topology;
      std::string topologyName;

      if (targetTopologyType == "line") {
        topology = compiler::CouplingMap::createLine(physicalQubits);
        topologyName = "Line (" + std::to_string(physicalQubits) + " physical qubits, " +
                       std::to_string(topology.edgeCount()) + " coupling edges)";
      } else if (targetTopologyType == "ring") {
        topology = compiler::CouplingMap::createRing(physicalQubits);
        topologyName = "Ring (" + std::to_string(physicalQubits) + " physical qubits, " +
                       std::to_string(topology.edgeCount()) + " coupling edges)";
      } else if (targetTopologyType == "star") {
        topology = compiler::CouplingMap::createStar(physicalQubits, 0);
        topologyName = "Star (" + std::to_string(physicalQubits) + " physical qubits, center 0, " +
                       std::to_string(topology.edgeCount()) + " coupling edges)";
      } else if (targetTopologyType == "grid") {
        std::size_t cols = 2;
        std::size_t rows = (physicalQubits + cols - 1) / cols;
        topology = compiler::CouplingMap::createGrid(rows, cols);
        topologyName = "Grid (" + std::to_string(rows) + "x" + std::to_string(cols) + ", " +
                       std::to_string(topology.qubitCount()) + " physical qubits, " +
                       std::to_string(topology.edgeCount()) + " coupling edges)";
      } else if (targetTopologyType == "all") {
        topology = compiler::CouplingMap::createFullyConnected(physicalQubits);
        topologyName = "FullyConnected (" + std::to_string(physicalQubits) + " physical qubits, " +
                       std::to_string(topology.edgeCount()) + " coupling edges)";
      } else {
        std::cerr << "Unknown target topology: " << targetTopologyType << "\n";
        return 1;
      }

      compiler::RoutingResult routing = compiler::routeCircuit(instructions, qubitCount, topology);
      if (!routing.success) {
        std::cerr << "Routing error: " << routing.errorMessage << "\n";
        return 1;
      }

      std::cout << "\nTopology Routing:\n"
                << "  Target topology: " << topologyName
                << "\n  Inserted SWAPs: " << routing.insertedSwapCount
                << " (decomposed into " << (3 * routing.insertedSwapCount) << " CNOTs)"
                << "\n  Logical 2Q gates: " << routing.logicalTwoQubitCount
                << " -> Routed 2Q gates: " << routing.routedTwoQubitCount
                << "\n  Routing 2Q overhead: +" << std::fixed << std::setprecision(1)
                << (routing.routingOverhead * 100.0) << "%"
                << "\n  Logical depth: " << routing.logicalDepth << " layers"
                << " -> Routed depth: " << routing.routedDepth << " layers\n";
    }
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
