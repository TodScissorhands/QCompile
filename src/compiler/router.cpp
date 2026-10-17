#include "compiler/router.h"

#include "compiler/circuit_dag.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace compiler {

QubitMapping::QubitMapping(std::size_t logicalCount, std::size_t physicalCount) {
  if (logicalCount > physicalCount) {
    throw std::invalid_argument("QubitMapping: logical qubit count exceeds physical hardware capacity");
  }

  logicalToPhysical.resize(logicalCount);
  physicalToLogical.assign(physicalCount, -1);

  for (std::size_t l = 0; l < logicalCount; ++l) {
    logicalToPhysical[l] = static_cast<int>(l);
    physicalToLogical[l] = static_cast<int>(l);
  }
}

void QubitMapping::applySwap(std::size_t u, std::size_t v) {
  if (u >= physicalToLogical.size() || v >= physicalToLogical.size() || u == v) {
    return;
  }

  int l_u = physicalToLogical[u];
  int l_v = physicalToLogical[v];

  std::swap(physicalToLogical[u], physicalToLogical[v]);

  if (l_u >= 0 && static_cast<std::size_t>(l_u) < logicalToPhysical.size()) {
    logicalToPhysical[l_u] = static_cast<int>(v);
  }
  if (l_v >= 0 && static_cast<std::size_t>(l_v) < logicalToPhysical.size()) {
    logicalToPhysical[l_v] = static_cast<int>(u);
  }
}

RoutingResult routeCircuit(
    const std::vector<Instruction> &logicalInstructions,
    std::size_t logicalQubitCount,
    const CouplingMap &topology
) {
  RoutingResult result;

  // 1. Hardware capacity validation
  if (logicalQubitCount > topology.qubitCount()) {
    result.success = false;
    result.errorMessage = "logical qubit count (" + std::to_string(logicalQubitCount) +
                          ") exceeds physical hardware capacity (" +
                          std::to_string(topology.qubitCount()) + ")";
    return result;
  }

  // 2. Upfront instruction validation
  for (const auto &inst : logicalInstructions) {
    if (inst.type == InstructionType::ApplyGate) {
      if (inst.gateName != "H" && inst.gateName != "X" && inst.gateName != "Z" &&
          inst.gateName != "RZ" && inst.gateName != "CNOT") {
        result.success = false;
        result.errorMessage = "unsupported gate for routing: " + inst.gateName;
        return result;
      }

      for (int q : inst.qubitIndices) {
        if (q < 0 || static_cast<std::size_t>(q) >= logicalQubitCount) {
          result.success = false;
          result.errorMessage = "logical qubit index out of bounds: " + std::to_string(q);
          return result;
        }
      }

      if (inst.gateName == "CNOT") {
        if (inst.qubitIndices.size() != 2 || inst.qubitIndices[0] == inst.qubitIndices[1]) {
          result.success = false;
          result.errorMessage = "malformed CNOT operands";
          return result;
        }
      } else {
        if (inst.qubitIndices.size() != 1) {
          result.success = false;
          result.errorMessage = "malformed single-qubit gate operands";
          return result;
        }
      }
    } else if (inst.type == InstructionType::Measure) {
      for (int q : inst.qubitIndices) {
        if (q < 0 || static_cast<std::size_t>(q) >= logicalQubitCount) {
          result.success = false;
          result.errorMessage = "logical measurement qubit index out of bounds: " + std::to_string(q);
          return result;
        }
      }
    }
  }

  // 3. Compute logical input metrics
  result.logicalInstructionCount = logicalInstructions.size();
  for (const auto &inst : logicalInstructions) {
    if (inst.type == InstructionType::ApplyGate) {
      result.logicalGateCount++;
      if (inst.gateName == "CNOT") {
        result.logicalTwoQubitCount++;
      }
    }
  }
  result.logicalDepth = CircuitDag::fromInstructions(logicalInstructions, logicalQubitCount).criticalPathDepth();

  // 4. Routing execution
  QubitMapping mapping(logicalQubitCount, topology.qubitCount());
  result.initialMapping = mapping;

  std::vector<Instruction> routed;
  routed.reserve(logicalInstructions.size());
  std::size_t insertedSwaps = 0;

  auto emitSwap = [&](std::size_t u, std::size_t v) {
    // Conceptual SWAP decomposed into exactly 3 primitive CNOTs
    routed.push_back(Instruction{InstructionType::ApplyGate, "CNOT", {static_cast<int>(u), static_cast<int>(v)}, 0.0});
    routed.push_back(Instruction{InstructionType::ApplyGate, "CNOT", {static_cast<int>(v), static_cast<int>(u)}, 0.0});
    routed.push_back(Instruction{InstructionType::ApplyGate, "CNOT", {static_cast<int>(u), static_cast<int>(v)}, 0.0});
    ++insertedSwaps;
    mapping.applySwap(u, v);
  };

  for (const auto &inst : logicalInstructions) {
    if (inst.type == InstructionType::ApplyGate) {
      if (inst.qubitIndices.size() == 1) {
        int l = inst.qubitIndices[0];
        int p = mapping.logicalToPhysical[l];
        routed.push_back(Instruction{InstructionType::ApplyGate, inst.gateName, {p}, inst.parameterValue});
      } else if (inst.gateName == "CNOT") {
        int l_c = inst.qubitIndices[0];
        int l_t = inst.qubitIndices[1];
        int p_c = mapping.logicalToPhysical[l_c];
        int p_t = mapping.logicalToPhysical[l_t];

        if (topology.hasEdge(static_cast<std::size_t>(p_c), static_cast<std::size_t>(p_t))) {
          // Already adjacent and topology-compliant
          routed.push_back(Instruction{InstructionType::ApplyGate, "CNOT", {p_c, p_t}, 0.0});
        } else {
          // Non-adjacent: find shortest path via BFS
          std::vector<std::size_t> path;
          bool pathFound = topology.findShortestPath(static_cast<std::size_t>(p_c), static_cast<std::size_t>(p_t), path);
          if (!pathFound || path.size() < 2) {
            result.success = false;
            result.errorMessage = "no physical coupling path exists between physical qubits " +
                                  std::to_string(p_c) + " and " + std::to_string(p_t);
            return result;
          }

          // Insert decomposed SWAPs along the path to bring control adjacent to target
          const std::size_t k = path.size() - 1;
          for (std::size_t i = 0; i < k - 1; ++i) {
            emitSwap(path[i], path[i + 1]);
          }

          // Now control is at physical wire path[k - 1], target is at path[k] = p_t
          int new_p_c = static_cast<int>(path[k - 1]);
          int new_p_t = static_cast<int>(path[k]);

          // Emit physical CNOT interaction on adjacent coupling edge
          routed.push_back(Instruction{InstructionType::ApplyGate, "CNOT", {new_p_c, new_p_t}, 0.0});
        }
      }
    } else if (inst.type == InstructionType::Measure) {
      // Translate logical measurement operands to their current physical wires
      std::vector<int> physicalMeasured;
      physicalMeasured.reserve(inst.qubitIndices.size());
      for (int l : inst.qubitIndices) {
        physicalMeasured.push_back(mapping.logicalToPhysical[l]);
      }
      routed.push_back(Instruction{InstructionType::Measure, "measure", physicalMeasured, 0.0});
    }
  }

  // 5. Post-route validation and metrics calculation
  for (const auto &inst : routed) {
    if (inst.type == InstructionType::ApplyGate && inst.qubitIndices.size() == 2) {
      if (!topology.hasEdge(static_cast<std::size_t>(inst.qubitIndices[0]),
                            static_cast<std::size_t>(inst.qubitIndices[1]))) {
        result.success = false;
        result.errorMessage = "internal error: emitted non-compliant physical interaction";
        return result;
      }
    }
  }

  result.finalMapping = mapping;
  result.routedInstructions = std::move(routed);
  result.insertedSwapCount = insertedSwaps;
  result.routedInstructionCount = result.routedInstructions.size();

  for (const auto &inst : result.routedInstructions) {
    if (inst.type == InstructionType::ApplyGate) {
      result.routedGateCount++;
      if (inst.gateName == "CNOT") {
        result.routedTwoQubitCount++;
      }
    }
  }

  // Universal accounting identities
  if (result.logicalTwoQubitCount > 0) {
    result.routingOverhead = static_cast<double>(3 * result.insertedSwapCount) /
                             static_cast<double>(result.logicalTwoQubitCount);
  } else {
    result.routingOverhead = 0.0;
  }

  result.routedDepth = CircuitDag::fromInstructions(result.routedInstructions, topology.qubitCount()).criticalPathDepth();
  result.success = true;
  return result;
}

} // namespace compiler
