#include "test_support.h"

#include "compiler/circuit_dag.h"
#include "compiler/compiler_facade.h"
#include "simulator/simulator.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

namespace test {

void testCircuitDag() {
  // =========================================================================
  // AC1: DAG Construction, Dependency Tracking, Topological Order & Cycle Rejection
  // =========================================================================
  {
    compiler::CircuitDag dag(3);
    require(dag.empty(), "new DAG is empty");
    require(dag.nodeCount() == 0, "node count is 0");
    require(dag.criticalPathDepth() == 0, "empty DAG has depth 0");
    require(!dag.hasCycle(), "empty DAG has no cycles");

    std::vector<std::size_t> order;
    require(dag.topologicalSort(order), "topological sort succeeds on empty DAG");
    require(order.empty(), "order is empty");

    // Add nodes manually
    Instruction h0 = makeInstruction(InstructionType::ApplyGate, "H", {0});
    Instruction x1 = makeInstruction(InstructionType::ApplyGate, "X", {1});
    Instruction cnot01 = makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1});

    std::size_t id0 = dag.addNode(h0);
    std::size_t id1 = dag.addNode(x1);
    std::size_t id2 = dag.addNode(cnot01);

    require(id0 == 0 && id1 == 1 && id2 == 2, "sequential node IDs assigned");
    require(dag.nodeCount() == 3, "node count is 3");

    // Add valid causal edges: id0 -> id2 (q0), id1 -> id2 (q1)
    require(dag.addEdge(id0, id2), "edge id0 -> id2 added successfully");
    require(dag.addEdge(id1, id2), "edge id1 -> id2 added successfully");
    require(dag.hasEdge(id0, id2), "hasEdge confirms id0 -> id2");
    require(dag.hasEdge(id1, id2), "hasEdge confirms id1 -> id2");
    require(!dag.hasEdge(id2, id0), "hasEdge confirms no back-edge");

    // Idempotent duplicate edge addition
    require(dag.addEdge(id0, id2), "duplicate edge addition returns true idempotently");

    // Cycle rejection: adding id2 -> id0 would form a cycle (id0 -> id2 -> id0)
    require(!dag.addEdge(id2, id0), "cycle-inducing edge id2 -> id0 is rejected");
    require(!dag.hasEdge(id2, id0), "graph remains acyclic without back-edge");
    require(!dag.hasCycle(), "DAG hasCycle is false after cycle rejection");

    // Self-loop rejection
    require(!dag.addEdge(id0, id0), "self-loop edge is rejected");

    // Out-of-bounds node edge rejection
    require(!dag.addEdge(id0, 999), "edge to non-existent node rejected");

    // Topological order validation
    require(dag.topologicalSort(order), "topologicalSort succeeds on valid DAG");
    require(order.size() == 3, "topological order contains all 3 nodes");
    // In order, id0 and id1 must both appear before id2
    auto pos0 = std::find(order.begin(), order.end(), id0);
    auto pos1 = std::find(order.begin(), order.end(), id1);
    auto pos2 = std::find(order.begin(), order.end(), id2);
    require(pos0 < pos2, "id0 appears before id2 in topological sort");
    require(pos1 < pos2, "id1 appears before id2 in topological sort");
  }

  // =========================================================================
  // AC2: Parallel ASAP Layer Assignment & Critical-Path Depth Verification
  // =========================================================================
  {
    // 1. Single gate
    {
      std::vector<Instruction> single{makeInstruction(InstructionType::ApplyGate, "H", {0})};
      auto dag = compiler::CircuitDag::fromInstructions(single, 1);
      require(dag.instructionCount() == 1, "single gate count is 1");
      require(dag.criticalPathDepth() == 1, "single gate depth is 1");
      auto layers = dag.layerSchedule();
      require(layers.size() == 1, "layer schedule has 1 layer");
      require(layers[0].size() == 1, "layer 0 has 1 gate");
    }

    // 2. N parallel gates on disjoint qubits -> Depth 1
    {
      std::vector<Instruction> parallelGates{
          makeInstruction(InstructionType::ApplyGate, "H", {0}),
          makeInstruction(InstructionType::ApplyGate, "X", {1}),
          makeInstruction(InstructionType::ApplyGate, "Z", {2}),
          makeInstruction(InstructionType::ApplyGate, "RZ", {3}, 0.5)
      };
      auto dag = compiler::CircuitDag::fromInstructions(parallelGates, 4);
      require(dag.instructionCount() == 4, "parallel instruction count is 4");
      require(dag.twoQubitGateCount() == 0, "parallel two-qubit count is 0");
      require(dag.criticalPathDepth() == 1, "4 parallel gates on disjoint qubits have depth 1");
      auto layers = dag.layerSchedule();
      require(layers.size() == 1, "layer schedule has exactly 1 layer");
      require(layers[0].size() == 4, "layer 0 contains all 4 parallel gates");
    }

    // 3. N serial gates on the same qubit -> Depth N
    {
      std::vector<Instruction> serialGates{
          makeInstruction(InstructionType::ApplyGate, "H", {0}),
          makeInstruction(InstructionType::ApplyGate, "X", {0}),
          makeInstruction(InstructionType::ApplyGate, "Z", {0}),
          makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 1.2)
      };
      auto dag = compiler::CircuitDag::fromInstructions(serialGates, 1);
      require(dag.instructionCount() == 4, "serial instruction count is 4");
      require(dag.criticalPathDepth() == 4, "4 serial gates on same qubit have depth 4");
      auto layers = dag.layerSchedule();
      require(layers.size() == 4, "layer schedule has 4 distinct layers");
      for (std::size_t i = 0; i < 4; ++i) {
        require(layers[i].size() == 1, "each layer has exactly 1 gate in serial circuit");
      }
    }

    // 4. Standard Bell circuit: H(q0), CNOT(q0, q1) -> Depth 2
    {
      std::vector<Instruction> bell{
          makeInstruction(InstructionType::ApplyGate, "H", {0}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1})
      };
      auto dag = compiler::CircuitDag::fromInstructions(bell, 2);
      require(dag.instructionCount() == 2, "Bell instruction count is 2");
      require(dag.twoQubitGateCount() == 1, "Bell two-qubit count is 1");
      require(dag.criticalPathDepth() == 2, "Bell circuit has depth 2");
      auto layers = dag.layerSchedule();
      require(layers.size() == 2, "Bell has 2 layers");
      require(layers[0].size() == 1 && layers[0][0].gateName == "H", "layer 0 is H");
      require(layers[1].size() == 1 && layers[1][0].gateName == "CNOT", "layer 1 is CNOT");
    }

    // 5. Invariant proof: In every layer of layerSchedule, instructions have disjoint qubits
    {
      std::vector<Instruction> complexCircuit{
          makeInstruction(InstructionType::ApplyGate, "H", {0}),
          makeInstruction(InstructionType::ApplyGate, "H", {1}),
          makeInstruction(InstructionType::ApplyGate, "H", {2}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
          makeInstruction(InstructionType::ApplyGate, "RZ", {2}, 0.7),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {1, 2}),
          makeInstruction(InstructionType::ApplyGate, "X", {0})
      };
      auto dag = compiler::CircuitDag::fromInstructions(complexCircuit, 3);
      auto layers = dag.layerSchedule();
      for (std::size_t layerIdx = 0; layerIdx < layers.size(); ++layerIdx) {
        const auto &layer = layers[layerIdx];
        std::vector<int> touchedQubits;
        for (const auto &inst : layer) {
          for (int q : inst.qubitIndices) {
            auto it = std::find(touchedQubits.begin(), touchedQubits.end(), q);
            require(it == touchedQubits.end(),
                    "layer invariant holds: no two gates in the same layer share a qubit");
            touchedQubits.push_back(q);
          }
        }
      }
    }
  }

  // =========================================================================
  // AC3: Three-Gate Disjoint Cancellation Example (H q0; X q1; H q0 -> X q1)
  // =========================================================================
  {
    std::vector<Instruction> circuit{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "X", {1}),
        makeInstruction(InstructionType::ApplyGate, "H", {0})
    };

    auto rawDag = compiler::CircuitDag::fromInstructions(circuit, 2);
    require(rawDag.instructionCount() == 3, "raw instruction count is 3");
    // Layer 0 has H(0) and X(1); Layer 1 has H(0) -> Raw parallel depth is 2
    require(rawDag.criticalPathDepth() == 2, "raw parallel circuit depth is 2 (not 3)");

    auto optDag = rawDag.optimize();
    require(optDag.instructionCount() == 1, "optimized instruction count is 1");
    require(optDag.criticalPathDepth() == 1, "optimized parallel circuit depth is 1");

    auto optInstructions = optDag.toInstructions();
    require(optInstructions.size() == 1, "optimized sequence has 1 instruction");
    require(optInstructions[0].gateName == "X" && optInstructions[0].qubitIndices == std::vector<int>{1},
            "remaining gate is X on q1");

    // Simulator equivalence
    Simulator simRaw(2);
    Simulator simOpt(2);
    simRaw.run(circuit);
    simOpt.run(optInstructions);
    requireSameState(simRaw.stateVector(), simOpt.stateVector(),
                     "statevectors match exactly for H(0)-X(1)-H(0) -> X(1)");
  }

  // =========================================================================
  // AC4: Commutation Through CNOT Control and Target
  // =========================================================================
  {
    // 1. Control phase commutation: RZ(theta) q0; CNOT q0, q1; RZ(phi) q0
    // Commutes through control and merges to CNOT q0, q1; RZ(theta + phi) q0
    {
      std::vector<Instruction> controlCircuit{
          makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 0.6),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
          makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 0.8)
      };
      auto rawDag = compiler::CircuitDag::fromInstructions(controlCircuit, 2);
      auto optDag = rawDag.optimize();
      auto optInsts = optDag.toInstructions();
      require(optInsts.size() == 2, "RZ on control commutes and merges: 3 -> 2 instructions");

      Simulator simRaw(2);
      Simulator simOpt(2);
      simRaw.run(controlCircuit);
      simOpt.run(optInsts);
      requireSameState(simRaw.stateVector(), simOpt.stateVector(),
                       "control RZ commutation preserves exact statevector");
    }

    // 2. Control Z commutation: Z q0; CNOT q0, q1; Z q0 cancels to CNOT q0, q1
    {
      std::vector<Instruction> zControlCircuit{
          makeInstruction(InstructionType::ApplyGate, "Z", {0}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
          makeInstruction(InstructionType::ApplyGate, "Z", {0})
      };
      auto optDag = compiler::CircuitDag::fromInstructions(zControlCircuit, 2).optimize();
      auto optInsts = optDag.toInstructions();
      require(optInsts.size() == 1, "Z on control commutes through CNOT and cancels: 3 -> 1");
      require(optInsts[0].gateName == "CNOT", "remaining gate is CNOT");

      Simulator simRaw(2);
      Simulator simOpt(2);
      simRaw.run(zControlCircuit);
      simOpt.run(optInsts);
      requireSameState(simRaw.stateVector(), simOpt.stateVector(),
                       "control Z cancellation preserves exact statevector");
    }

    // 3. Target bit-flip commutation: X q1; CNOT q0, q1; X q1 cancels to CNOT q0, q1
    {
      std::vector<Instruction> xTargetCircuit{
          makeInstruction(InstructionType::ApplyGate, "X", {1}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
          makeInstruction(InstructionType::ApplyGate, "X", {1})
      };
      auto optDag = compiler::CircuitDag::fromInstructions(xTargetCircuit, 2).optimize();
      auto optInsts = optDag.toInstructions();
      require(optInsts.size() == 1, "X on target commutes through CNOT and cancels: 3 -> 1");
      require(optInsts[0].gateName == "CNOT", "remaining gate is CNOT");

      Simulator simRaw(2);
      Simulator simOpt(2);
      simRaw.run(xTargetCircuit);
      simOpt.run(optInsts);
      requireSameState(simRaw.stateVector(), simOpt.stateVector(),
                       "target X cancellation preserves exact statevector");
    }

    // 4. Non-commuting gate combinations MUST NOT commute or cancel:
    // 4a. X on control qubit (X_0 CNOT_01 != CNOT_01 X_0)
    {
      std::vector<Instruction> xControlCircuit{
          makeInstruction(InstructionType::ApplyGate, "X", {0}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
          makeInstruction(InstructionType::ApplyGate, "X", {0})
      };
      auto optDag = compiler::CircuitDag::fromInstructions(xControlCircuit, 2).optimize();
      require(optDag.instructionCount() == 3, "X on control qubit does not commute through CNOT");
    }

    // 4b. Z on target qubit (Z_1 CNOT_01 != CNOT_01 Z_1)
    {
      std::vector<Instruction> zTargetCircuit{
          makeInstruction(InstructionType::ApplyGate, "Z", {1}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
          makeInstruction(InstructionType::ApplyGate, "Z", {1})
      };
      auto optDag = compiler::CircuitDag::fromInstructions(zTargetCircuit, 2).optimize();
      require(optDag.instructionCount() == 3, "Z on target qubit does not commute through CNOT");
    }

    // 4c. H on control qubit (H_0 CNOT_01 != CNOT_01 H_0)
    {
      std::vector<Instruction> hControlCircuit{
          makeInstruction(InstructionType::ApplyGate, "H", {0}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
          makeInstruction(InstructionType::ApplyGate, "H", {0})
      };
      auto optDag = compiler::CircuitDag::fromInstructions(hControlCircuit, 2).optimize();
      require(optDag.instructionCount() == 3, "H on control qubit does not commute through CNOT");
    }

    // 4d. H on target qubit (H_1 CNOT_01 != CNOT_01 H_1)
    {
      std::vector<Instruction> hTargetCircuit{
          makeInstruction(InstructionType::ApplyGate, "H", {1}),
          makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
          makeInstruction(InstructionType::ApplyGate, "H", {1})
      };
      auto optDag = compiler::CircuitDag::fromInstructions(hTargetCircuit, 2).optimize();
      require(optDag.instructionCount() == 3, "H on target qubit does not commute through CNOT");
    }
  }

  // =========================================================================
  // AC5: Semantic Measurement Invariant
  // =========================================================================
  {
    // Gates cannot commute across measurement on measured qubits
    Instruction h0 = makeInstruction(InstructionType::ApplyGate, "H", {0});
    Instruction meas0 = makeInstruction(InstructionType::Measure, "measure", {0});
    std::vector<Instruction> measCircuit{h0, meas0, h0};

    auto optDag = compiler::CircuitDag::fromInstructions(measCircuit, 1).optimize();
    require(optDag.instructionCount() == 3, "H does not cancel across measurement on same qubit");
    auto optInsts = optDag.toInstructions();
    require(sameInstructions(measCircuit, optInsts),
            "instructions around measurement barrier are preserved in original order");

    // Multi-qubit measurement acts as barrier on all measured qubits
    Instruction measAll = makeInstruction(InstructionType::Measure, "measure", {0, 1});
    std::vector<Instruction> multiMeasCircuit{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "X", {1}),
        measAll,
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "X", {1})
    };
    auto optMulti = compiler::CircuitDag::fromInstructions(multiMeasCircuit, 2).optimize();
    require(optMulti.instructionCount() == 5,
            "multi-qubit measurement strictly isolates pre- and post-measurement operations");
  }

  // =========================================================================
  // AC6: Exact Unitary Equivalence & Statevector Invariance
  // =========================================================================
  {
    // Test on all 4 computational basis inputs for a 2-qubit circuit with commuting gates:
    // RZ(0.5) q0; X q1; CNOT q0, q1; X q1; RZ(-0.5) q0
    // - X(q1) commutes through target of CNOT(q0, q1) and cancels: X*X = I
    // - RZ(0.5) q0 commutes through control of CNOT(q0, q1) and merges with RZ(-0.5) q0 to cancel: I
    // Net result: CNOT(q0, q1)
    std::vector<Instruction> testEquiv{
        makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 0.5),
        makeInstruction(InstructionType::ApplyGate, "X", {1}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
        makeInstruction(InstructionType::ApplyGate, "X", {1}),
        makeInstruction(InstructionType::ApplyGate, "RZ", {0}, -0.5)
    };

    auto optDag = compiler::CircuitDag::fromInstructions(testEquiv, 2).optimize();
    auto optInsts = optDag.toInstructions();
    require(optInsts.size() == 1, "testEquiv simplifies completely to CNOT");
    require(optInsts[0].gateName == "CNOT", "simplified gate is CNOT");
    // Verify exact unitary matrix equivalence: U_opt |j> == U_orig |j> for all j in {0,1,2,3}
    for (int basis = 0; basis < 4; ++basis) {
      // Prepare basis state by applying bit flips
      std::vector<Instruction> prep;
      if (basis & 1) prep.push_back(makeInstruction(InstructionType::ApplyGate, "X", {0}));
      if (basis & 2) prep.push_back(makeInstruction(InstructionType::ApplyGate, "X", {1}));

      std::vector<Instruction> runOrig = prep;
      runOrig.insert(runOrig.end(), testEquiv.begin(), testEquiv.end());

      std::vector<Instruction> runOpt = prep;
      runOpt.insert(runOpt.end(), optInsts.begin(), optInsts.end());

      Simulator simOrig(2);
      Simulator simOpt(2);
      simOrig.run(runOrig);
      simOpt.run(runOpt);

      requireSameState(simOrig.stateVector(), simOpt.stateVector(),
                       "exact unitary equivalence verified on basis state " + std::to_string(basis));
    }
  }

  // =========================================================================
  // AC7: Comprehensive Instruction Set Coverage & Idempotence
  // =========================================================================
  {
    // Circuit with all supported gate types: H, X, Z, RZ, CNOT
    std::vector<Instruction> fullSet{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 0.5),
        makeInstruction(InstructionType::ApplyGate, "Z", {0}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
        makeInstruction(InstructionType::ApplyGate, "X", {1}),
        makeInstruction(InstructionType::ApplyGate, "H", {2})
    };

    auto dag1 = compiler::CircuitDag::fromInstructions(fullSet, 3);
    auto dagOpt1 = dag1.optimize();
    auto dagOpt2 = dagOpt1.optimize();

    require(sameInstructions(dagOpt1.toInstructions(), dagOpt2.toInstructions()),
            "CircuitDag optimization is strictly idempotent");

    // Out-of-range instruction preservation
    Instruction outOfRange = makeInstruction(InstructionType::ApplyGate, "H", {99});
    auto dagOOR = compiler::CircuitDag::fromInstructions({outOfRange, outOfRange}, 2).optimize();
    require(dagOOR.instructionCount() == 2, "out-of-range instructions are preserved");
  }

  // =========================================================================
  // AC8: Compiler Facade & CompilationResult Integration
  // =========================================================================
  {
    const std::string dslSource =
        "qubit q0, q1\n"
        "H q0\n"
        "X q1\n"
        "H q0\n";

    // 1. Optimized compilation
    compiler::CompilerOptions optOptions;
    optOptions.optimize = true;
    compiler::CompilationResult optRes = compiler::compileProgram(dslSource, optOptions);

    require(optRes.success, "compilation succeeds");
    require(optRes.rawGateCount == 3, "rawGateCount is 3");
    require(optRes.rawDepth == 2, "rawDepth is 2 (parallel ASAP layers)");
    require(optRes.optimizedGateCount == 1, "optimizedGateCount is 1");
    require(optRes.optimizedDepth == 1, "optimizedDepth is 1");
    require(optRes.rawTwoQubitCount == 0, "rawTwoQubitCount is 0");
    require(optRes.optimizedTwoQubitCount == 0, "optimizedTwoQubitCount is 0");

    // 2. Unoptimized compilation toggle
    compiler::CompilerOptions unoptOptions;
    unoptOptions.optimize = false;
    compiler::CompilationResult unoptRes = compiler::compileProgram(dslSource, unoptOptions);

    require(unoptRes.success, "unoptimized compilation succeeds");
    require(unoptRes.instructions.size() == 3, "unoptimized instruction count is 3");
    require(unoptRes.optimizedDepth == 2, "unoptimized depth equals raw depth 2");
    require(unoptRes.optimizedGateCount == 3, "unoptimized gate count is 3");
  }
}

} // namespace test
