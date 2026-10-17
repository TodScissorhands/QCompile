#include "test_support.h"

#include <algorithm>
#include <vector>

namespace test {

void testInstructionOptimization() {
  Program declaredCircuit;
  const SemanticResult declaredSemantics =
      check("qubit q0, q1\n", declaredCircuit);
  require(declaredSemantics.errors.empty(),
          "optimizer test circuit has validated qubit declarations");
  const std::size_t qubitCount = declaredSemantics.qubitNames.size();
  const auto optimize = [qubitCount](
                           const std::vector<Instruction> &instructions) {
    return optimizeInstructions(instructions, qubitCount);
  };

  const Instruction h = makeInstruction(InstructionType::ApplyGate, "H", {0});
  const Instruction x = makeInstruction(InstructionType::ApplyGate, "X", {0});
  const Instruction z = makeInstruction(InstructionType::ApplyGate, "Z", {0});
  const Instruction cnot =
      makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1});
  for (const Instruction &gate : {h, x, z, cnot}) {
    require(optimize({gate, gate}).empty(),
            gate.gateName + " identical adjacent pair cancels");
  }
  const int lastValidQubit = static_cast<int>(qubitCount) - 1;
  const Instruction boundaryH =
      makeInstruction(InstructionType::ApplyGate, "H", {lastValidQubit});
  require(optimize({boundaryH, boundaryH}).empty(),
          "gate at the last declared qubit remains eligible");
  const int firstOutOfRangeQubit = static_cast<int>(qubitCount);
  const Instruction outOfRangeH =
      makeInstruction(InstructionType::ApplyGate, "H", {firstOutOfRangeQubit});
  require(sameInstructions(optimize({outOfRangeH, outOfRangeH}),
                           {outOfRangeH, outOfRangeH}),
          "identical H gates at a positive out-of-range index are preserved");
  const Instruction outOfRangeCnot = makeInstruction(
      InstructionType::ApplyGate, "CNOT", {0, firstOutOfRangeQubit});
  require(sameInstructions(optimize({outOfRangeCnot, outOfRangeCnot}),
                           {outOfRangeCnot, outOfRangeCnot}),
          "CNOT with a positive out-of-range operand is preserved");
  require(sameInstructions(optimize({h, h, h}), {h}),
          "repeated H-H-H cancellation leaves H");
  require(sameInstructions(
              optimize({h, makeInstruction(
                                           InstructionType::ApplyGate, "H", {1})}),
              {h, makeInstruction(InstructionType::ApplyGate, "H", {1})}),
          "single-qubit gates on different qubits do not cancel");
  const Instruction reversedCnot =
      makeInstruction(InstructionType::ApplyGate, "CNOT", {1, 0});
  require(sameInstructions(optimize({cnot, reversedCnot}),
                           {cnot, reversedCnot}),
          "CNOT with reversed operands does not cancel");
  require(sameInstructions(optimize({h, x, h}), {h, x, h}),
          "non-adjacent gates do not cancel across an intervening gate");

  const Instruction rz =
      makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 0.75);

  // Phase 8: Consecutive RZ rotation merging regression tests
  // 1. RZ(1.0), RZ(2.0) on the same qubit merging to RZ(3.0)
  const Instruction rz1 = makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 1.0);
  const Instruction rz2 = makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 2.0);
  const Instruction rz3 = makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 3.0);
  require(sameInstructions(optimize({rz1, rz2}), {rz3}),
          "consecutive RZ(1.0), RZ(2.0) on same qubit merge to RZ(3.0)");

  // 2. RZ(1.5), RZ(-1.5) cancelling to no rotation
  const Instruction rzPos = makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 1.5);
  const Instruction rzNeg = makeInstruction(InstructionType::ApplyGate, "RZ", {0}, -1.5);
  require(optimize({rzPos, rzNeg}).empty(),
          "consecutive RZ(1.5), RZ(-1.5) cancel completely");

  // 3. Rotations on different qubits remaining separate
  const Instruction rzQ0 = makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 1.0);
  const Instruction rzQ1 = makeInstruction(InstructionType::ApplyGate, "RZ", {1}, 2.0);
  require(sameInstructions(optimize({rzQ0, rzQ1}), {rzQ0, rzQ1}),
          "RZ on different qubits remain separate");

  // 4. Non-adjacent rotations remaining unmerged
  require(sameInstructions(optimize({rz1, h, rz2}), {rz1, h, rz2}),
          "non-adjacent RZ rotations across intervening gate remain unmerged");

  // 5. Measurement acting as a strict ordering barrier
  const Instruction measBarrier = makeInstruction(InstructionType::Measure, "measure", {0});
  require(sameInstructions(optimize({rz1, measBarrier, rz2}), {rz1, measBarrier, rz2}),
          "measurement acts as strict barrier preventing RZ merging");

  // 6. Z and RZ not being incorrectly merged or cancelled
  require(sameInstructions(optimize({z, rz1}), {z, rz1}),
          "Z followed by RZ is not merged or cancelled");
  require(sameInstructions(optimize({rz1, z}), {rz1, z}),
          "RZ followed by Z is not merged or cancelled");

  // 7. Non-finite parameters preserved without producing invalid instructions
  const double infVal = std::numeric_limits<double>::infinity();
  const double nanVal = std::numeric_limits<double>::quiet_NaN();
  const Instruction rzInf = makeInstruction(InstructionType::ApplyGate, "RZ", {0}, infVal);
  const Instruction rzNan = makeInstruction(InstructionType::ApplyGate, "RZ", {0}, nanVal);
  require(sameInstructions(optimize({rzInf, rz1}), {rzInf, rz1}),
          "RZ with inf parameter is preserved without merging");
  require(sameInstructions(optimize({rz1, rzNan}), {rz1, rzNan}),
          "RZ with nan parameter is preserved without merging");

  // Rotations near non-zero multiples of 2pi are preserved (not eliminated)
  const double twoPi = 2.0 * std::acos(-1.0);
  const Instruction rzTwoPi = makeInstruction(InstructionType::ApplyGate, "RZ", {0}, twoPi);
  require(sameInstructions(optimize({rzTwoPi}), {rzTwoPi}),
          "RZ near 2pi is preserved and not eliminated");

  // 8. State-vector equivalence between optimized and unoptimized circuits
  {
    const std::vector<Instruction> rzCircuit{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 0.4),
        makeInstruction(InstructionType::ApplyGate, "RZ", {0}, 0.8),
        makeInstruction(InstructionType::ApplyGate, "RZ", {0}, -0.2),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
        makeInstruction(InstructionType::ApplyGate, "RZ", {1}, 1.2),
        makeInstruction(InstructionType::ApplyGate, "RZ", {1}, -1.2)
    };
    const std::vector<Instruction> optimizedRzCircuit = optimize(rzCircuit);
    require(optimizedRzCircuit.size() == 3, "circuit reduced to H, merged RZ(1.0), and CNOT");
    Simulator simUnopt(2);
    Simulator simOpt(2);
    simUnopt.run(rzCircuit);
    simOpt.run(optimizedRzCircuit);
    requireSameState(simUnopt.stateVector(), simOpt.stateVector(),
                     "state vectors match between unoptimized and RZ-merged circuits");
  }
  const Instruction measurement =
      makeInstruction(InstructionType::Measure, "measure", {1, 0});
  const Instruction secondMeasurement =
      makeInstruction(InstructionType::Measure, "measure", {0});
  const std::vector<Instruction> aroundMeasurement{
      h, h, measurement, x, x, secondMeasurement};
  const std::vector<Instruction> optimizedMeasurement =
      optimize(aroundMeasurement);
  require(sameInstructions(optimizedMeasurement,
                           {measurement, secondMeasurement}),
          "cancellations preserve measurement barriers, order, and qubit lists");
  require(std::count_if(optimizedMeasurement.begin(), optimizedMeasurement.end(),
                        [](const Instruction &instruction) {
                          return instruction.type == InstructionType::Measure;
                        }) == 2,
          "measurement count is preserved");

  const Instruction unknown =
      makeInstruction(InstructionType::ApplyGate, "UNKNOWN", {0});
  const Instruction malformedH =
      makeInstruction(InstructionType::ApplyGate, "H", {0, 1});
  const Instruction malformedCnot =
      makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 0});
  const Instruction malformedParameter =
      makeInstruction(InstructionType::ApplyGate, "H", {0}, 1.0);
  const Instruction malformedIndex =
      makeInstruction(InstructionType::ApplyGate, "X", {-1});
  for (const Instruction &invalid : {unknown, malformedH, malformedCnot,
                                     malformedParameter, malformedIndex}) {
    require(sameInstructions(optimize({invalid, invalid}),
                             {invalid, invalid}),
            "unknown or malformed instructions are preserved");
  }

  const std::vector<Instruction> fixedPoint{
      h, x, makeInstruction(InstructionType::ApplyGate, "Z", {1}), cnot, rz,
      measurement, h, unknown, malformedH};
  const std::vector<Instruction> optimized =
      optimize(fixedPoint);
  require(sameInstructions(optimize(optimized), optimized),
          "optimization pass is idempotent");
  for (std::size_t index = 1; index < optimized.size(); ++index) {
    require(!optimize(
                 {optimized[index - 1], optimized[index]}).empty(),
            "no eligible adjacent pair remains");
  }

  const std::vector<std::vector<Instruction>> contexts{
      {makeInstruction(InstructionType::ApplyGate, "X", {1}), h, h, cnot, cnot},
      {makeInstruction(InstructionType::ApplyGate, "H", {0}), x, x, rz,
       rz, makeInstruction(InstructionType::ApplyGate, "H", {1}), z, z}};
  for (const std::vector<Instruction> &circuit : contexts) {
    const std::vector<Instruction> reduced = optimize(circuit);
    Simulator originalSimulator(2);
    Simulator optimizedSimulator(2);
    originalSimulator.run(circuit);
    optimizedSimulator.run(reduced);
    requireSameState(originalSimulator.stateVector(),
                     optimizedSimulator.stateVector(),
                     "optimized and original final state vectors agree");
  }
}

} // namespace test
