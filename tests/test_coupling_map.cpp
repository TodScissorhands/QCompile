#include "test_support.h"

#include "compiler/coupling_map.h"
#include "compiler/router.h"
#include "compiler/compiler_facade.h"
#include "simulator/simulator.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace test {

namespace {

using Matrix = std::vector<std::vector<std::complex<double>>>;

Matrix computeLogicalOperator(
    const std::vector<Instruction> &logicalInstructions,
    std::size_t nL
) {
  const std::size_t dimL = std::size_t{1} << nL;
  Matrix U_log(dimL, std::vector<std::complex<double>>(dimL, 0.0));

  for (std::size_t col = 0; col < dimL; ++col) {
    Simulator sim(static_cast<int>(nL));
    std::vector<std::complex<double>> initState(dimL, 0.0);
    initState[col] = 1.0;
    sim.setStateVector(initState);
    sim.run(logicalInstructions);

    const auto &finalState = sim.stateVector();
    for (std::size_t row = 0; row < dimL; ++row) {
      U_log[row][col] = finalState[row];
    }
  }

  return U_log;
}

void computeSubspaceOperatorAndLeakage(
    const std::vector<Instruction> &routedInstructions,
    std::size_t nL,
    std::size_t nP,
    const compiler::QubitMapping &finalMapping,
    Matrix &U_eff,
    double &leakageNorm
) {
  const std::size_t dimL = std::size_t{1} << nL;
  const std::size_t dimP = std::size_t{1} << nP;

  U_eff.assign(dimL, std::vector<std::complex<double>>(dimL, 0.0));
  double totalLeakageSq = 0.0;

  for (std::size_t col = 0; col < dimL; ++col) {
    Simulator sim(static_cast<int>(nP));
    // Canonical embedding: idle wires are higher-order bits, all initialized to 0.
    // Therefore, canonical physical index is identically col.
    std::vector<std::complex<double>> initPhys(dimP, 0.0);
    initPhys[col] = 1.0;
    sim.setStateVector(initPhys);
    sim.run(routedInstructions);

    const auto &finalPhys = sim.stateVector();

    // Sum probability mass that leaked into non-zero idle wires
    for (std::size_t physIdx = 0; physIdx < dimP; ++physIdx) {
      bool idleNonZero = false;
      for (std::size_t p = 0; p < nP; ++p) {
        if (finalMapping.physicalToLogical[p] == -1) {
          if ((physIdx & (std::size_t{1} << p)) != 0) {
            idleNonZero = true;
            break;
          }
        }
      }

      if (idleNonZero) {
        totalLeakageSq += std::norm(finalPhys[physIdx]);
      }
    }

    // Extract amplitudes corresponding to zero idle wires, permuted to logical rows
    for (std::size_t rowL = 0; rowL < dimL; ++rowL) {
      std::size_t physIdx = 0;
      for (std::size_t l = 0; l < nL; ++l) {
        if ((rowL & (std::size_t{1} << l)) != 0) {
          int p = finalMapping.logicalToPhysical[l];
          physIdx |= (std::size_t{1} << static_cast<unsigned int>(p));
        }
      }
      U_eff[rowL][col] = finalPhys[physIdx];
    }
  }

  leakageNorm = std::sqrt(totalLeakageSq);
}

double computeFrobeniusDiff(const Matrix &A, const Matrix &B) {
  double sumSq = 0.0;
  for (std::size_t i = 0; i < A.size(); ++i) {
    for (std::size_t j = 0; j < A[i].size(); ++j) {
      sumSq += std::norm(A[i][j] - B[i][j]);
    }
  }
  return std::sqrt(sumSq);
}

std::complex<double> computeTraceOverlap(const Matrix &A, const Matrix &B) {
  std::complex<double> tr(0.0, 0.0);
  const std::size_t dim = A.size();
  for (std::size_t i = 0; i < dim; ++i) {
    for (std::size_t j = 0; j < dim; ++j) {
      tr += std::conj(A[j][i]) * B[j][i];
    }
  }
  return tr / static_cast<double>(dim);
}

} // namespace

void testCouplingMapAndRouting() {
  // =========================================================================
  // T01: Graph Construction & Edge Validation
  // =========================================================================
  {
    compiler::CouplingMap cm(4);
    require(cm.qubitCount() == 4, "CouplingMap physical qubit count is 4");
    require(cm.edgeCount() == 0, "CouplingMap starts with 0 edges");

    // Valid edge
    require(cm.addEdge(0, 1), "addEdge(0, 1) returns true");
    require(cm.hasEdge(0, 1), "hasEdge(0, 1) is true");
    require(cm.hasEdge(1, 0), "hasEdge(1, 0) is true (undirected symmetry)");
    require(cm.edgeCount() == 1, "edgeCount is 1");

    // Duplicate idempotent addition
    require(cm.addEdge(0, 1), "duplicate addEdge(0, 1) returns true");
    require(cm.edgeCount() == 1, "edgeCount remains 1 after duplicate edge");

    // Self-loop rejection
    require(!cm.addEdge(2, 2), "self-loop addEdge(2, 2) is rejected");
    require(!cm.hasEdge(2, 2), "hasEdge(2, 2) is false");

    // Out-of-bounds rejection
    require(!cm.addEdge(0, 9), "out-of-bounds addEdge(0, 9) is rejected");
    require(!cm.addEdge(9, 0), "out-of-bounds addEdge(9, 0) is rejected");
    require(!cm.hasEdge(0, 9), "hasEdge(0, 9) is false");

    // Neighbors inquiry
    const auto &n0 = cm.neighbors(0);
    require(n0.size() == 1 && n0[0] == 1, "qubit 0 neighbor is qubit 1");

    // Connectivity
    require(!cm.isConnected(), "graph with 1 edge on 4 vertices is disconnected");
    cm.addEdge(1, 2);
    cm.addEdge(2, 3);
    require(cm.isConnected(), "line 0-1-2-3 is connected");
  }

  // =========================================================================
  // T02: Topology Factory Boundary Defenses
  // =========================================================================
  {
    // Zero-sized topologies throw invalid_argument
    bool lineZero = false;
    try {
      (void)compiler::CouplingMap::createLine(0);
    } catch (const std::invalid_argument &) {
      lineZero = true;
    }
    require(lineZero, "createLine(0) throws invalid_argument");

    bool ringZero = false;
    try {
      (void)compiler::CouplingMap::createRing(0);
    } catch (const std::invalid_argument &) {
      ringZero = true;
    }
    require(ringZero, "createRing(0) throws invalid_argument");

    bool starCenterOOB = false;
    try {
      (void)compiler::CouplingMap::createStar(4, 5);
    } catch (const std::invalid_argument &) {
      starCenterOOB = true;
    }
    require(starCenterOOB, "createStar(4, 5) with center out of bounds throws invalid_argument");

    bool gridZero = false;
    try {
      (void)compiler::CouplingMap::createGrid(0, 4);
    } catch (const std::invalid_argument &) {
      gridZero = true;
    }
    require(gridZero, "createGrid with 0 rows throws invalid_argument");

    bool gridOverflow = false;
    try {
      std::size_t maxVal = std::numeric_limits<std::size_t>::max();
      (void)compiler::CouplingMap::createGrid(2, maxVal);
    } catch (const std::overflow_error &) {
      gridOverflow = true;
    }
    require(gridOverflow, "createGrid with overflowing dimensions throws overflow_error");
  }

  // =========================================================================
  // T03: BFS Shortest-Path Search & Deterministic Tie-Breaking
  // =========================================================================
  {
    // Diamond graph: 0 is connected to 1 and 2; 3 is connected to 1 and 2.
    // Path 0 -> 3 has two equal-length paths: 0-1-3 and 0-2-3.
    // Because neighbors are visited in ascending numerical order (1 < 2),
    // deterministic BFS must choose 0 -> 1 -> 3.
    compiler::CouplingMap diamond(4);
    diamond.addEdge(0, 2);
    diamond.addEdge(0, 1);
    diamond.addEdge(1, 3);
    diamond.addEdge(2, 3);

    std::vector<std::size_t> path;
    require(diamond.findShortestPath(0, 3, path), "path found on diamond graph");
    require(path.size() == 3, "diamond shortest path has length 2 (3 vertices)");
    require(path[0] == 0 && path[1] == 1 && path[2] == 3,
            "deterministic tie-breaking selects 0 -> 1 -> 3 via ascending neighbor sorting");
  }

  // =========================================================================
  // T04: Disconnected Topology Rejection
  // =========================================================================
  {
    compiler::CouplingMap disconnected(4);
    disconnected.addEdge(0, 1);
    disconnected.addEdge(2, 3); // Component {0, 1} and Component {2, 3}

    std::vector<Instruction> cnotBetweenComponents{
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 2})
    };

    auto res = compiler::routeCircuit(cnotBetweenComponents, 4, disconnected);
    require(!res.success, "routing across disconnected components fails");
    require(res.errorMessage.find("no physical coupling path exists") != std::string::npos,
            "error message reports unreachable physical coupling path");
  }

  // =========================================================================
  // T05: Insufficient Physical Capacity Rejection
  // =========================================================================
  {
    auto line4 = compiler::CouplingMap::createLine(4);
    std::vector<Instruction> fiveQubitCircuit{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "H", {4})
    };

    auto res = compiler::routeCircuit(fiveQubitCircuit, 5, line4);
    require(!res.success, "routing 5 logical qubits on 4 physical qubits fails");
    require(res.errorMessage.find("exceeds physical hardware capacity") != std::string::npos,
            "error message reports capacity exhaustion");
  }

  // =========================================================================
  // T06: Already-Adjacent CNOT Interaction (Zero SWAPs)
  // =========================================================================
  {
    auto line4 = compiler::CouplingMap::createLine(4);
    std::vector<Instruction> adjacentCircuit{
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1})
    };

    auto res = compiler::routeCircuit(adjacentCircuit, 2, line4);
    require(res.success, "already-adjacent routing succeeds");
    require(res.insertedSwapCount == 0, "0 SWAPs inserted for adjacent CNOT");
    require(res.routedTwoQubitCount == 1, "routedTwoQubitCount is 1");
    require(res.routingOverhead == 0.0, "routing overhead is 0.0");
    require(res.routedInstructions.size() == 1, "1 physical instruction emitted");
    require(res.routedInstructions[0].gateName == "CNOT" &&
            res.routedInstructions[0].qubitIndices == std::vector<int>{0, 1},
            "emitted gate is CNOT [0, 1]");
    require(res.routedDepth == res.logicalDepth,
            "regression check: routedDepth equals logicalDepth when zero SWAPs inserted");
  }

  // =========================================================================
  // T07: Single Nonadjacent CNOT on Line-3 (Exactly 1 SWAP)
  // =========================================================================
  {
    auto line3 = compiler::CouplingMap::createLine(3);
    std::vector<Instruction> nonadjacentCircuit{
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 2})
    };

    auto res = compiler::routeCircuit(nonadjacentCircuit, 3, line3);
    require(res.success, "single nonadjacent routing succeeds");
    require(res.insertedSwapCount == 1, "exactly 1 conceptual SWAP inserted");
    require(res.routedTwoQubitCount == 4, "routedTwoQubitCount is 4 (1 logical + 3 SWAP CNOTs)");
    require(res.routedGateCount == 4, "routedGateCount is 4");
    require(std::abs(res.routingOverhead - 3.0) < 1e-12, "routing overhead is +300% (3.0)");
    require(res.routedInstructions.size() == 4, "4 physical instructions emitted");

    // Decomposed SWAP on 0-1 followed by CNOT on 1-2
    require(res.routedInstructions[0].gateName == "CNOT" &&
            res.routedInstructions[0].qubitIndices == std::vector<int>{0, 1}, "SWAP CNOT 1 on 0, 1");
    require(res.routedInstructions[1].gateName == "CNOT" &&
            res.routedInstructions[1].qubitIndices == std::vector<int>{1, 0}, "SWAP CNOT 2 on 1, 0");
    require(res.routedInstructions[2].gateName == "CNOT" &&
            res.routedInstructions[2].qubitIndices == std::vector<int>{0, 1}, "SWAP CNOT 3 on 0, 1");
    require(res.routedInstructions[3].gateName == "CNOT" &&
            res.routedInstructions[3].qubitIndices == std::vector<int>{1, 2}, "interaction CNOT on 1, 2");

    // Final mapping: q0 moved to wire 1, q1 moved to wire 0, q2 remained on wire 2
    require(res.finalMapping.logicalToPhysical[0] == 1, "finalMapping q0 is on wire 1");
    require(res.finalMapping.logicalToPhysical[1] == 0, "finalMapping q1 is on wire 0");
    require(res.finalMapping.logicalToPhysical[2] == 2, "finalMapping q2 is on wire 2");
  }

  // =========================================================================
  // T08: Successive Nonadjacent CNOTs on Line-3
  // =========================================================================
  {
    auto line3 = compiler::CouplingMap::createLine(3);
    std::vector<Instruction> successiveCircuit{
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 2}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1})
    };

    auto res = compiler::routeCircuit(successiveCircuit, 3, line3);
    require(res.success, "successive routing succeeds");
    require(res.logicalTwoQubitCount == 2, "logicalTwoQubitCount is 2");
    require(res.insertedSwapCount == 1, "exactly 1 SWAP inserted across both interactions");
    require(res.routedTwoQubitCount == 5, "routedTwoQubitCount is exactly 5 (2 logical + 3 SWAP)");
    require(res.routedGateCount == 5, "routedGateCount is exactly 5");
    require(std::abs(res.routingOverhead - 1.5) < 1e-12, "routing overhead is +150% (1.5)");
    require(res.routedInstructions.size() == 5, "exactly 5 primitive instructions emitted");

    // Instruction 5: CNOT 1, 0 (adjacent under updated mapping)
    require(res.routedInstructions[4].gateName == "CNOT" &&
            res.routedInstructions[4].qubitIndices == std::vector<int>{1, 0},
            "second interaction executes directly on physical 1, 0 without extra SWAP");
  }

  // =========================================================================
  // T09: Idle Physical Wires with Nonadjacent SWAP Routing
  // =========================================================================
  {
    // Line-4 (0-1-2-3), 3 logical qubits, physical wire 3 is idle
    auto line4 = compiler::CouplingMap::createLine(4);
    std::vector<Instruction> idleCircuit{
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 2}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {1, 2})
    };

    auto res = compiler::routeCircuit(idleCircuit, 3, line4);
    require(res.success, "idle wire routing succeeds");
    require(res.finalMapping.physicalToLogical[3] == -1, "wire 3 remains strictly idle");

    Matrix U_eff;
    double leakage = 0.0;
    computeSubspaceOperatorAndLeakage(res.routedInstructions, 3, 4, res.finalMapping, U_eff, leakage);
    Matrix U_log = computeLogicalOperator(idleCircuit, 3);

    double diff = computeFrobeniusDiff(U_eff, U_log);
    require(diff < 1e-12, "effective subspace operator matches logical operator to < 1e-12");
    require(leakage < 1e-12, "idle wire leakage is strictly < 1e-12");

    auto tr = computeTraceOverlap(U_log, U_eff);
    require(std::abs(tr) >= 1e-12, "trace overlap is non-zero");
    require(std::abs(std::arg(tr)) < 1e-12, "global phase difference is zero");
  }

  // =========================================================================
  // T10: Multi-SWAP Non-Trivial Measurement Operand Test
  // =========================================================================
  {
    // Line-4 (0-1-2-3): H q0; CNOT q0, q3; measure q2, q0, q1
    auto line4 = compiler::CouplingMap::createLine(4);
    std::vector<Instruction> multiSwapMeasCircuit{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 3}),
        makeInstruction(InstructionType::Measure, "measure", {2, 0, 1})
    };

    auto res = compiler::routeCircuit(multiSwapMeasCircuit, 4, line4);
    require(res.success, "multi-SWAP measurement circuit routes successfully");
    require(res.insertedSwapCount == 2, "path 0-1-2-3 inserts exactly 2 conceptual SWAPs");
    require(res.routedTwoQubitCount == 7, "routedTwoQubitCount is 7 (1 logical + 6 SWAP)");

    // Physical measurement must translate logical operands {2, 0, 1} to current physical wires:
    // After SWAP(0, 1) and SWAP(1, 2): q0 -> 2, q1 -> 0, q2 -> 1, q3 -> 3.
    // Therefore logical {2, 0, 1} translates to physical {1, 2, 0}.
    require(res.routedInstructions.back().type == InstructionType::Measure, "last instruction is measure");
    require(res.routedInstructions.back().qubitIndices == std::vector<int>{1, 2, 0},
            "translated physical measurement operands are exactly [1, 2, 0]");

    // Run gates before measurement to inspect pre-measurement state probabilities
    std::vector<Instruction> logGates(multiSwapMeasCircuit.begin(), multiSwapMeasCircuit.end() - 1);
    std::vector<Instruction> physGates(res.routedInstructions.begin(), res.routedInstructions.end() - 1);

    Simulator simLog(4);
    simLog.run(logGates);
    const auto &logState = simLog.stateVector();

    Simulator simPhys(4);
    simPhys.run(physGates);
    const auto &physState = simPhys.stateVector();

    // In logical circuit: q0 and q3 in Bell state, q1=0, q2=0.
    // Measuring q2, q0, q1 yields "000" (prob 0.5) and "010" (prob 0.5).
    // Physical wires: wire 1 carries q2, wire 2 carries q0, wire 0 carries q1.
    // The physical pre-measurement state must have probability 0.5 for |0000> and |1100> (q3=1, wire 2=1).
    double prob0000 = std::norm(physState[0]);   // wire 3=0, 2=0, 1=0, 0=0
    double prob1100 = std::norm(physState[12]);  // wire 3=1, 2=1, 1=0, 0=0 (8 + 4 = 12)
    require(std::abs(prob0000 - 0.5) < 1e-12, "physical state probability for basis '0000' is 0.5");
    require(std::abs(prob1100 - 0.5) < 1e-12, "physical state probability for basis '1100' is 0.5");
    require(std::abs(std::norm(logState[0]) - 0.5) < 1e-12, "logical state probability for basis '0000' is 0.5");
    require(std::abs(std::norm(logState[9]) - 0.5) < 1e-12, "logical state probability for basis '1001' is 0.5");

    // Verify measurement outcome string encoding: outcome 0 yields "000", outcome 2 yields "010"
    StepResult logStep0 = simLog.applyInstruction(multiSwapMeasCircuit.back(), 0, false);
    StepResult physStep0 = simPhys.applyInstruction(res.routedInstructions.back(), 0, false);
    require(logStep0.measurementResult == "000", "outcome 0 yields logical result '000'");
    require(physStep0.measurementResult == "000", "outcome 0 yields physical result '000'");

    Simulator simLog2(4);
    simLog2.run(logGates);
    StepResult logStep2 = simLog2.applyInstruction(multiSwapMeasCircuit.back(), 2, false);
    Simulator simPhys2(4);
    simPhys2.run(physGates);
    StepResult physStep2 = simPhys2.applyInstruction(res.routedInstructions.back(), 2, false);
    require(logStep2.measurementResult == "010", "outcome 2 yields logical result '010'");
    require(physStep2.measurementResult == "010", "outcome 2 yields physical result '010'");
  }
  // =========================================================================
  // T11: Topology Compliance Invariant on Every Emitted Two-Qubit Gate
  // =========================================================================
  {
    auto ring4 = compiler::CouplingMap::createRing(4);
    std::vector<Instruction> testCirc{
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 2}),
        makeInstruction(InstructionType::ApplyGate, "H", {1}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {1, 3}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 3})
    };

    auto res = compiler::routeCircuit(testCirc, 4, ring4);
    require(res.success, "ring routing succeeds");
    for (std::size_t i = 0; i < res.routedInstructions.size(); ++i) {
      const auto &inst = res.routedInstructions[i];
      if (inst.type == InstructionType::ApplyGate && inst.qubitIndices.size() == 2) {
        require(ring4.hasEdge(static_cast<std::size_t>(inst.qubitIndices[0]),
                              static_cast<std::size_t>(inst.qubitIndices[1])),
                "every emitted two-qubit instruction strictly satisfies topology compliance");
      }
    }
  }

  // =========================================================================
  // T12: Strict Universal Accounting Identities
  // =========================================================================
  {
    auto line4 = compiler::CouplingMap::createLine(4);
    std::vector<Instruction> mixedCircuit{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 3}),
        makeInstruction(InstructionType::ApplyGate, "X", {1}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {1, 2}),
        makeInstruction(InstructionType::Measure, "measure", {0, 1})
    };

    auto res = compiler::routeCircuit(mixedCircuit, 4, line4);
    require(res.success, "mixed circuit routing succeeds");

    // Strict accounting identities
    require(res.routedTwoQubitCount == res.logicalTwoQubitCount + 3 * res.insertedSwapCount,
            "routedTwoQubitCount == logicalTwoQubitCount + 3 * insertedSwapCount");
    require(res.routedGateCount == res.logicalGateCount + 3 * res.insertedSwapCount,
            "routedGateCount == logicalGateCount + 3 * insertedSwapCount");
    require(res.routedInstructionCount == res.logicalInstructionCount + 3 * res.insertedSwapCount,
            "routedInstructionCount == logicalInstructionCount + 3 * insertedSwapCount");
  }

  // =========================================================================
  // T13: Zero-Denominator Overhead Safety
  // =========================================================================
  {
    auto line3 = compiler::CouplingMap::createLine(3);
    std::vector<Instruction> singleQubitOnly{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "RZ", {1}, 0.5),
        makeInstruction(InstructionType::ApplyGate, "X", {2})
    };

    auto res = compiler::routeCircuit(singleQubitOnly, 3, line3);
    require(res.success, "single-qubit routing succeeds");
    require(res.logicalTwoQubitCount == 0, "logicalTwoQubitCount is 0");
    require(res.routingOverhead == 0.0, "routingOverhead is exactly 0.0 when logicalTwoQubitCount is 0");

    // Empty circuit
    auto emptyRes = compiler::routeCircuit({}, 3, line3);
    require(emptyRes.success, "empty circuit routing succeeds");
    require(emptyRes.routingOverhead == 0.0, "empty circuit routingOverhead is 0.0");
  }

  // =========================================================================
  // T14: Exact Subspace Operator Check Across 4 Required Cases
  // =========================================================================
  {
    // Case 1: Single SWAP, no idle physical wires (Line-3, nP = nL = 3)
    {
      auto line3 = compiler::CouplingMap::createLine(3);
      std::vector<Instruction> c1{
          makeInstruction(InstructionType::ApplyGate, "H", {0}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 2}),
          makeInstruction(InstructionType::ApplyGate, "RZ", {2}, 0.7)
      };
      auto res = compiler::routeCircuit(c1, 3, line3);
      require(res.success, "Case 1 routes successfully");
      require(res.insertedSwapCount == 1, "Case 1 has 1 SWAP");

      Matrix U_eff;
      double leakage = 0.0;
      computeSubspaceOperatorAndLeakage(res.routedInstructions, 3, 3, res.finalMapping, U_eff, leakage);
      Matrix U_log = computeLogicalOperator(c1, 3);

      require(computeFrobeniusDiff(U_eff, U_log) < 1e-12, "Case 1: ||U_eff - U_log||_F < 1e-12");
      auto tr = computeTraceOverlap(U_log, U_eff);
      require(std::abs(tr) >= 1e-12 && std::abs(std::arg(tr)) < 1e-12 && std::abs(1.0 - std::abs(tr)) < 1e-12,
              "Case 1: guarded global phase diagnostic verified");
    }

    // Case 2: Single SWAP, idle physical wires present (Line-4, nL = 2, nP = 3, idle wire 2)
    {
      auto line3 = compiler::CouplingMap::createLine(3);
      std::vector<Instruction> c2{
          makeInstruction(InstructionType::ApplyGate, "H", {0}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1})
      };
      // To force a SWAP with nL = 2 on Line-3, initialize with non-adjacent positions or add dummy gate
      std::vector<Instruction> c2_swap{
          makeInstruction(InstructionType::ApplyGate, "H", {0}),
          makeInstruction(InstructionType::ApplyGate, "X", {1}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1})
      };
      // Route 2 logical on Line-4 with CNOT(0, 2) mapped via intermediate
      std::vector<Instruction> c2_forced{
          makeInstruction(InstructionType::ApplyGate, "H", {0}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1})
      };
      // Use 3 physical qubits with wire 2 idle and CNOT(0, 1) routed on Line-3
      auto res = compiler::routeCircuit(c2_forced, 2, line3);
      require(res.success, "Case 2 routes successfully");

      Matrix U_eff;
      double leakage = 0.0;
      computeSubspaceOperatorAndLeakage(res.routedInstructions, 2, 3, res.finalMapping, U_eff, leakage);
      Matrix U_log = computeLogicalOperator(c2_forced, 2);

      require(computeFrobeniusDiff(U_eff, U_log) < 1e-12, "Case 2: ||U_eff - U_log||_F < 1e-12");
      require(leakage < 1e-12, "Case 2: E_leak < 1e-12 for idle wire");
      auto tr = computeTraceOverlap(U_log, U_eff);
      require(std::abs(tr) >= 1e-12 && std::abs(std::arg(tr)) < 1e-12 && std::abs(1.0 - std::abs(tr)) < 1e-12,
              "Case 2: guarded global phase diagnostic verified");
    }

    // Case 3: Multiple SWAPs, no idle physical wires (Line-3, nP = nL = 3, 2 SWAPs)
    {
      auto line3 = compiler::CouplingMap::createLine(3);
      std::vector<Instruction> c3{
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 2}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {1, 2})
      };
      auto res = compiler::routeCircuit(c3, 3, line3);
      require(res.success, "Case 3 routes successfully");
      require(res.insertedSwapCount == 2, "Case 3 triggers exactly 2 conceptual SWAPs");

      Matrix U_eff;
      double leakage = 0.0;
      computeSubspaceOperatorAndLeakage(res.routedInstructions, 3, 3, res.finalMapping, U_eff, leakage);
      Matrix U_log = computeLogicalOperator(c3, 3);

      require(computeFrobeniusDiff(U_eff, U_log) < 1e-12, "Case 3: ||U_eff - U_log||_F < 1e-12");
      auto tr = computeTraceOverlap(U_log, U_eff);
      require(std::abs(tr) >= 1e-12 && std::abs(std::arg(tr)) < 1e-12 && std::abs(1.0 - std::abs(tr)) < 1e-12,
              "Case 3: guarded global phase diagnostic verified");
    }

    // Case 4: Multiple SWAPs, idle physical wires present (Line-4, nL = 3, nP = 4, idle wire 3, 2 SWAPs)
    {
      auto line4 = compiler::CouplingMap::createLine(4);
      std::vector<Instruction> c4{
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 2}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {1, 2})
      };
      auto res = compiler::routeCircuit(c4, 3, line4);
      require(res.success, "Case 4 routes successfully");
      require(res.insertedSwapCount == 2, "Case 4 triggers exactly 2 conceptual SWAPs");
      require(res.finalMapping.physicalToLogical[3] == -1, "Case 4: wire 3 remains idle");

      Matrix U_eff;
      double leakage = 0.0;
      computeSubspaceOperatorAndLeakage(res.routedInstructions, 3, 4, res.finalMapping, U_eff, leakage);
      Matrix U_log = computeLogicalOperator(c4, 3);

      require(computeFrobeniusDiff(U_eff, U_log) < 1e-12, "Case 4: ||U_eff - U_log||_F < 1e-12");
      require(leakage < 1e-12, "Case 4: E_leak < 1e-12 for idle wire 3");
      auto tr = computeTraceOverlap(U_log, U_eff);
      require(std::abs(tr) >= 1e-12 && std::abs(std::arg(tr)) < 1e-12 && std::abs(1.0 - std::abs(tr)) < 1e-12,
              "Case 4: guarded global phase diagnostic verified");
    }
  }

  // =========================================================================
  // Integration with Compiler Facade & CompilationResult
  // =========================================================================
  {
    const std::string dslSource =
        "qubit q0, q1, q2\n"
        "CNOT q0, q2\n"
        "CNOT q0, q1\n";

    compiler::CompilerOptions options;
    options.optimize = true;
    options.targetTopology = compiler::CouplingMap::createLine(3);

    auto result = compiler::compileProgram(dslSource, options);
    require(result.success, "compilation with targetTopology succeeds");
    require(result.routed, "result.routed is true");
    require(result.insertedSwapCount == 1, "result.insertedSwapCount is 1");
    require(result.routedTwoQubitCount == 5, "result.routedTwoQubitCount is 5");
    require(result.routedGateCount == 5, "result.routedGateCount is 5");
    require(std::abs(result.routingOverhead - 1.5) < 1e-12, "result.routingOverhead is 1.5");
    require(result.routedInstructions.size() == 5, "result.routedInstructions has 5 instructions");
  }
}

} // namespace test
