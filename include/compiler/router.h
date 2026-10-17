#pragma once

#include "compiler/code_generator.h"
#include "compiler/coupling_map.h"

#include <cstddef>
#include <string>
#include <vector>

namespace compiler {

struct QubitMapping {
  std::vector<int> logicalToPhysical; // Size = logicalCount
  std::vector<int> physicalToLogical; // Size = physicalCount (-1 if idle)

  QubitMapping() = default;
  QubitMapping(std::size_t logicalCount, std::size_t physicalCount);

  // Apply a physical SWAP between physical wires u and v
  void applySwap(std::size_t u, std::size_t v);
};

struct RoutingResult {
  bool success = false;
  std::string errorMessage;
  std::vector<Instruction> routedInstructions;
  QubitMapping initialMapping;
  QubitMapping finalMapping;

  // Routing Metrics
  std::size_t logicalGateCount = 0;        // Unitary gates in post-opt input
  std::size_t logicalTwoQubitCount = 0;    // CNOTs in post-opt input
  std::size_t logicalInstructionCount = 0; // Total instructions in post-opt input
  std::size_t logicalDepth = 0;            // Critical-path depth of logical DAG

  std::size_t insertedSwapCount = 0;       // Conceptual SWAPs inserted
  std::size_t routedGateCount = 0;         // Unitary gates in routed output
  std::size_t routedTwoQubitCount = 0;     // CNOTs in routed output
  std::size_t routedInstructionCount = 0;  // Total instructions in routed output
  std::size_t routedDepth = 0;             // Critical-path depth of routed DAG
  double routingOverhead = 0.0;            // (routed2Q - log2Q) / max(1, log2Q)
};

// Top-level routing entry point
RoutingResult routeCircuit(
    const std::vector<Instruction> &logicalInstructions,
    std::size_t logicalQubitCount,
    const CouplingMap &topology
);

} // namespace compiler
