#pragma once

#include <string>
#include <vector>

// The enum tags the plain Statement record without using a class hierarchy.
enum class StatementType {
    QubitDeclaration,
    GateOperation
};

struct Statement {
    StatementType type;
    std::vector<std::string> names;
    std::string gate;
    std::vector<std::string> qubits;
};

struct Program {
    std::vector<Statement> statements;
};
