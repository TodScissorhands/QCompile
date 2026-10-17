#include "compiler/circuit_dag.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <unordered_map>

namespace compiler {

namespace {

bool isCandidateForCancellationOrMerge(const Instruction &inst, std::size_t qubitCount) {
  if (inst.type != InstructionType::ApplyGate) {
    return false;
  }

  if (inst.gateName == "H" || inst.gateName == "X" || inst.gateName == "Z") {
    return inst.qubitIndices.size() == 1 &&
           inst.qubitIndices[0] >= 0 &&
           static_cast<std::size_t>(inst.qubitIndices[0]) < qubitCount &&
           inst.parameterValue == 0.0;
  }

  if (inst.gateName == "CNOT") {
    return inst.qubitIndices.size() == 2 &&
           inst.qubitIndices[0] >= 0 &&
           inst.qubitIndices[1] >= 0 &&
           static_cast<std::size_t>(inst.qubitIndices[0]) < qubitCount &&
           static_cast<std::size_t>(inst.qubitIndices[1]) < qubitCount &&
           inst.qubitIndices[0] != inst.qubitIndices[1] &&
           inst.parameterValue == 0.0;
  }

  if (inst.gateName == "RZ") {
    return inst.qubitIndices.size() == 1 &&
           inst.qubitIndices[0] >= 0 &&
           static_cast<std::size_t>(inst.qubitIndices[0]) < qubitCount &&
           std::isfinite(inst.parameterValue);
  }

  return false;
}

bool canCancel(const Instruction &left, const Instruction &right, std::size_t qubitCount) {
  if (left.type != InstructionType::ApplyGate ||
      right.type != InstructionType::ApplyGate ||
      left.gateName != right.gateName ||
      left.qubitIndices != right.qubitIndices ||
      left.parameterValue != 0.0 ||
      right.parameterValue != 0.0) {
    return false;
  }

  if (left.gateName == "H" || left.gateName == "X" || left.gateName == "Z") {
    return left.qubitIndices.size() == 1 &&
           left.qubitIndices[0] >= 0 &&
           static_cast<std::size_t>(left.qubitIndices[0]) < qubitCount;
  }

  if (left.gateName == "CNOT") {
    return left.qubitIndices.size() == 2 &&
           left.qubitIndices[0] >= 0 &&
           left.qubitIndices[1] >= 0 &&
           static_cast<std::size_t>(left.qubitIndices[0]) < qubitCount &&
           static_cast<std::size_t>(left.qubitIndices[1]) < qubitCount &&
           left.qubitIndices[0] != left.qubitIndices[1];
  }

  return false;
}

bool canMergeRz(const Instruction &left, const Instruction &right, std::size_t qubitCount, double &mergedAngle) {
  if (left.type != InstructionType::ApplyGate ||
      right.type != InstructionType::ApplyGate ||
      left.gateName != "RZ" ||
      right.gateName != "RZ" ||
      left.qubitIndices != right.qubitIndices ||
      left.qubitIndices.size() != 1 ||
      left.qubitIndices[0] < 0 ||
      static_cast<std::size_t>(left.qubitIndices[0]) >= qubitCount) {
    return false;
  }

  if (!std::isfinite(left.parameterValue) || !std::isfinite(right.parameterValue)) {
    return false;
  }

  double sum = left.parameterValue + right.parameterValue;
  if (!std::isfinite(sum)) {
    return false;
  }

  mergedAngle = sum;
  return true;
}

} // namespace

bool canCommute(
    const Instruction &first,
    const Instruction &second,
    std::size_t qubitCount
) {
  // Measurements and non-ApplyGate operations never commute (strict semantic barrier)
  if (first.type != InstructionType::ApplyGate ||
      second.type != InstructionType::ApplyGate) {
    return false;
  }

  // Check for disjoint qubit support
  bool sharesQubit = false;
  for (int q1 : first.qubitIndices) {
    for (int q2 : second.qubitIndices) {
      if (q1 == q2) {
        sharesQubit = true;
        break;
      }
    }
    if (sharesQubit) {
      break;
    }
  }

  if (!sharesQubit) {
    // Operations on disjoint qubits commute unconditionally
    return true;
  }

  // If both are single-qubit gates on the same qubit:
  if (first.qubitIndices.size() == 1 && second.qubitIndices.size() == 1 &&
      first.qubitIndices[0] == second.qubitIndices[0]) {
    int q = first.qubitIndices[0];
    if (q < 0 || static_cast<std::size_t>(q) >= qubitCount) {
      return false; // out of range: do not commute
    }

    // Both are Z: commute
    if (first.gateName == "Z" && second.gateName == "Z") {
      return true;
    }

    // Both are RZ: commute if parameters are finite
    if (first.gateName == "RZ" && second.gateName == "RZ") {
      return std::isfinite(first.parameterValue) &&
             std::isfinite(second.parameterValue);
    }

    // Z and RZ commute on the same qubit if parameter is finite
    if ((first.gateName == "Z" && second.gateName == "RZ") ||
        (first.gateName == "RZ" && second.gateName == "Z")) {
      double rzParam = (first.gateName == "RZ") ? first.parameterValue : second.parameterValue;
      return std::isfinite(rzParam);
    }

    return false;
  }

  // Single-qubit gate and CNOT interaction
  const Instruction *single = nullptr;
  const Instruction *cnot = nullptr;
  if (first.qubitIndices.size() == 1 && second.gateName == "CNOT" && second.qubitIndices.size() == 2) {
    single = &first;
    cnot = &second;
  } else if (second.qubitIndices.size() == 1 && first.gateName == "CNOT" && first.qubitIndices.size() == 2) {
    single = &second;
    cnot = &first;
  }

  if (single != nullptr && cnot != nullptr) {
    int qSingle = single->qubitIndices[0];
    int control = cnot->qubitIndices[0];
    int target = cnot->qubitIndices[1];

    if (qSingle < 0 || control < 0 || target < 0 ||
        static_cast<std::size_t>(qSingle) >= qubitCount ||
        static_cast<std::size_t>(control) >= qubitCount ||
        static_cast<std::size_t>(target) >= qubitCount ||
        control == target) {
      return false;
    }

    // 1. Z or RZ on control qubit commutes with CNOT
    if (qSingle == control) {
      if (single->gateName == "Z") {
        return true;
      }
      if (single->gateName == "RZ" && std::isfinite(single->parameterValue)) {
        return true;
      }
    }

    // 2. X on target qubit commutes with CNOT
    if (qSingle == target) {
      if (single->gateName == "X" && single->parameterValue == 0.0) {
        return true;
      }
    }

    return false;
  }

  // Identical CNOTs commute
  if (first.gateName == "CNOT" && second.gateName == "CNOT" &&
      first.qubitIndices == second.qubitIndices) {
    return true;
  }

  // Conservative default: non-commuting
  return false;
}

CircuitDag::CircuitDag(std::size_t qubitCount)
    : qubitCount_(qubitCount) {}

std::size_t CircuitDag::addNode(const Instruction &instruction) {
  DagNode node;
  node.id = nodes_.size();
  node.instruction = instruction;
  node.layer = 0;
  nodes_.push_back(node);
  return node.id;
}

bool CircuitDag::hasEdge(std::size_t fromId, std::size_t toId) const {
  if (fromId >= nodes_.size() || toId >= nodes_.size()) {
    return false;
  }
  const auto &succs = nodes_[fromId].successors;
  return std::find(succs.begin(), succs.end(), toId) != succs.end();
}

bool CircuitDag::addEdge(std::size_t fromId, std::size_t toId) {
  if (fromId >= nodes_.size() || toId >= nodes_.size() || fromId == toId) {
    return false;
  }

  if (hasEdge(fromId, toId)) {
    return true;
  }

  // Temporarily add edge
  nodes_[fromId].successors.push_back(toId);
  nodes_[toId].predecessors.push_back(fromId);

  // Check if edge introduced a cycle
  if (hasCycle()) {
    nodes_[fromId].successors.pop_back();
    nodes_[toId].predecessors.pop_back();
    return false;
  }

  computeLayers();
  return true;
}

bool CircuitDag::hasCycle() const {
  std::vector<std::size_t> order;
  return !topologicalSort(order);
}

bool CircuitDag::topologicalSort(std::vector<std::size_t> &order) const {
  order.clear();
  if (nodes_.empty()) {
    return true;
  }

  std::vector<std::size_t> inDegree(nodes_.size(), 0);
  for (const auto &node : nodes_) {
    for (std::size_t succ : node.successors) {
      if (succ < nodes_.size()) {
        inDegree[succ]++;
      }
    }
  }

  std::queue<std::size_t> q;
  for (std::size_t i = 0; i < nodes_.size(); ++i) {
    if (inDegree[i] == 0) {
      q.push(i);
    }
  }

  while (!q.empty()) {
    std::size_t curr = q.front();
    q.pop();
    order.push_back(curr);

    for (std::size_t succ : nodes_[curr].successors) {
      if (succ < nodes_.size()) {
        if (--inDegree[succ] == 0) {
          q.push(succ);
        }
      }
    }
  }

  return order.size() == nodes_.size();
}

void CircuitDag::computeLayers() {
  std::vector<std::size_t> order;
  if (!topologicalSort(order)) {
    return;
  }

  for (std::size_t id : order) {
    DagNode &curr = nodes_[id];
    int maxPredLayer = -1;
    for (std::size_t pred : curr.predecessors) {
      if (nodes_[pred].layer > maxPredLayer) {
        maxPredLayer = nodes_[pred].layer;
      }
    }
    curr.layer = maxPredLayer + 1;
  }
}

CircuitDag CircuitDag::fromInstructions(
    const std::vector<Instruction> &instructions,
    std::size_t qubitCount
) {
  CircuitDag dag(qubitCount);
  std::unordered_map<int, std::size_t> lastNodeOnWire;

  for (const auto &inst : instructions) {
    std::size_t nodeId = dag.addNode(inst);

    if (inst.qubitIndices.empty()) {
      // Connect to immediately preceding node if present to preserve global sequencing
      if (nodeId > 0) {
        dag.nodes_[nodeId - 1].successors.push_back(nodeId);
        dag.nodes_[nodeId].predecessors.push_back(nodeId - 1);
      }
    } else {
      for (int q : inst.qubitIndices) {
        auto it = lastNodeOnWire.find(q);
        if (it != lastNodeOnWire.end()) {
          std::size_t predId = it->second;
          if (predId != nodeId && !dag.hasEdge(predId, nodeId)) {
            dag.nodes_[predId].successors.push_back(nodeId);
            dag.nodes_[nodeId].predecessors.push_back(predId);
          }
        }
        lastNodeOnWire[q] = nodeId;
      }
    }
  }

  dag.computeLayers();
  return dag;
}

std::size_t CircuitDag::instructionCount() const {
  return nodes_.size();
}

std::size_t CircuitDag::gateCount(InstructionType type) const {
  std::size_t count = 0;
  for (const auto &n : nodes_) {
    if (n.instruction.type == type) {
      count++;
    }
  }
  return count;
}

std::size_t CircuitDag::twoQubitGateCount() const {
  std::size_t count = 0;
  for (const auto &n : nodes_) {
    if (n.instruction.type == InstructionType::ApplyGate &&
        n.instruction.qubitIndices.size() == 2) {
      count++;
    }
  }
  return count;
}

std::size_t CircuitDag::criticalPathDepth() const {
  if (nodes_.empty()) {
    return 0;
  }
  int maxL = 0;
  for (const auto &n : nodes_) {
    if (n.layer > maxL) {
      maxL = n.layer;
    }
  }
  return static_cast<std::size_t>(maxL + 1);
}

std::vector<std::vector<Instruction>> CircuitDag::layerSchedule() const {
  std::size_t depth = criticalPathDepth();
  if (depth == 0) {
    return {};
  }

  std::vector<std::vector<Instruction>> layers(depth);
  std::vector<std::size_t> order;
  if (!topologicalSort(order)) {
    return {};
  }

  for (std::size_t id : order) {
    const DagNode &n = nodes_[id];
    if (n.layer >= 0 && static_cast<std::size_t>(n.layer) < depth) {
      layers[n.layer].push_back(n.instruction);
    }
  }

  return layers;
}

std::vector<Instruction> CircuitDag::toInstructions() const {
  std::vector<std::size_t> order;
  if (!topologicalSort(order)) {
    return {};
  }

  std::vector<Instruction> result;
  result.reserve(order.size());
  for (std::size_t id : order) {
    result.push_back(nodes_[id].instruction);
  }
  return result;
}

CircuitDag CircuitDag::optimize() const {
  std::vector<Instruction> current = toInstructions();
  bool progress = true;

  while (progress) {
    progress = false;

    for (std::size_t i = 0; i < current.size() && !progress; ++i) {
      const Instruction &first = current[i];
      if (!isCandidateForCancellationOrMerge(first, qubitCount_)) {
        continue;
      }

      for (std::size_t j = i + 1; j < current.size(); ++j) {
        const Instruction &second = current[j];

        // Check if intervening gates all commute with first
        bool allInterveningCommute = true;
        for (std::size_t k = i + 1; k < j; ++k) {
          if (!canCommute(current[k], first, qubitCount_)) {
            allInterveningCommute = false;
            break;
          }
        }
        if (!allInterveningCommute) {
          // Blocked by non-commuting gate or barrier
          break;
        }

        // Gates can meet! Check for self-inverse cancellation
        if (canCancel(first, second, qubitCount_)) {
          current.erase(current.begin() + j);
          current.erase(current.begin() + i);
          progress = true;
          break;
        }

        // Check for RZ merging
        double mergedAngle = 0.0;
        if (canMergeRz(first, second, qubitCount_, mergedAngle)) {
          if (std::abs(mergedAngle) < 1e-12) {
            current.erase(current.begin() + j);
            current.erase(current.begin() + i);
          } else {
            current[j].parameterValue = mergedAngle;
            current.erase(current.begin() + i);
          }
          progress = true;
          break;
        }
      }
    }
  }

  return CircuitDag::fromInstructions(current, qubitCount_);
}

std::vector<Instruction> optimizeWithDag(
    const std::vector<Instruction> &instructions,
    std::size_t qubitCount
) {
  CircuitDag dag = CircuitDag::fromInstructions(instructions, qubitCount);
  CircuitDag optDag = dag.optimize();
  return optDag.toInstructions();
}

} // namespace compiler
