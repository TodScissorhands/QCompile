#pragma once

#include "compiler/code_generator.h"

#include <cstddef>
#include <string>
#include <vector>

namespace compiler {

struct DagNode {
  std::size_t id = 0;
  Instruction instruction;
  std::vector<std::size_t> predecessors;
  std::vector<std::size_t> successors;
  int layer = 0; // 0-based ASAP layer
};

class CircuitDag {
public:
  CircuitDag() = default;
  explicit CircuitDag(std::size_t qubitCount);

  // Construction from instruction sequence
  static CircuitDag fromInstructions(
      const std::vector<Instruction> &instructions,
      std::size_t qubitCount
  );

  // Graph inspection
  std::size_t qubitCount() const { return qubitCount_; }
  std::size_t nodeCount() const { return nodes_.size(); }
  bool empty() const { return nodes_.empty(); }
  const DagNode &node(std::size_t id) const { return nodes_.at(id); }
  const std::vector<DagNode> &nodes() const { return nodes_; }

  // Manual node/edge manipulation with cycle rejection
  std::size_t addNode(const Instruction &instruction);
  bool addEdge(std::size_t fromId, std::size_t toId);
  bool hasEdge(std::size_t fromId, std::size_t toId) const;

  // Topological sorting and cycle detection
  bool hasCycle() const;
  bool topologicalSort(std::vector<std::size_t> &order) const;

  // Circuit metrics
  std::size_t instructionCount() const;
  std::size_t gateCount(InstructionType type = InstructionType::ApplyGate) const;
  std::size_t twoQubitGateCount() const;
  std::size_t criticalPathDepth() const;

  // Parallel ASAP layer scheduling
  std::vector<std::vector<Instruction>> layerSchedule() const;

  // Serialization to instruction vector following topological order
  std::vector<Instruction> toInstructions() const;

  // Conservative commutation-aware optimization
  CircuitDag optimize() const;

private:
  void computeLayers();

  std::size_t qubitCount_ = 0;
  std::vector<DagNode> nodes_;
};

// Conservative commutation checking predicate
bool canCommute(
    const Instruction &first,
    const Instruction &second,
    std::size_t qubitCount
);

// Top-level optimizer function using CircuitDag commutation analysis
std::vector<Instruction> optimizeWithDag(
    const std::vector<Instruction> &instructions,
    std::size_t qubitCount
);

} // namespace compiler
