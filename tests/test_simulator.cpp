#include "test_support.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

namespace test {

void testCodeGeneration() {
  Program program;
  SemanticResult result = check(
      "qubit second, first\nlet theta = 2 + 3 * 4\nRZ(theta) first\n"
      "CNOT second, first\nmeasure first, second\n",
      program);
  require(result.errors.empty(), "valid code-generation input accepted");
  std::vector<Instruction> instructions = generateInstructions(program, result);
  require(instructions.size() == 3, "three executable instructions generated");
  require(instructions[0].type == InstructionType::ApplyGate &&
              instructions[0].gateName == "RZ" &&
              instructions[0].qubitIndices == std::vector<int>{1} &&
              std::abs(instructions[0].parameterValue - 14.0) < 1e-12,
          "RZ instruction has evaluated parameter and declaration-order operand");
  require(instructions[1].gateName == "CNOT" &&
              instructions[1].qubitIndices == std::vector<int>({0, 1}),
          "CNOT operands use declaration-order indices");
  require(instructions[2].type == InstructionType::Measure &&
              instructions[2].qubitIndices == std::vector<int>({1, 0}),
          "measurement instruction and operand order preserved");
}

void testSimulator() {
  std::string initial = run("qubit q0\nmeasure q0\n");
  require(initial.find("0: 1.000000") != std::string::npos,
          "initial state has all probability at zero");
  std::string x = run("qubit q0\nX q0\nmeasure q0\n");
  require(x.find("1: 1.000000") != std::string::npos, "X flips initial basis state");
  std::string h = run("qubit q0\nH q0\nmeasure q0\n");
  require(h.find("0: 0.500000") != std::string::npos &&
              h.find("1: 0.500000") != std::string::npos,
          "H produces equal basis probabilities");
  std::string bell = run("qubit q0, q1\nH q0\nCNOT q0, q1\nmeasure q0, q1\n");
  require(bell.find("00: 0.500000") != std::string::npos &&
              bell.find("11: 0.500000") != std::string::npos &&
              bell.find("01:") == std::string::npos &&
              bell.find("10:") == std::string::npos,
          "CNOT creates expected two-qubit correlated state");
  std::string rz = run("qubit q0\nH q0\nRZ(1.2) q0\nmeasure q0\n");
  require(rz.find("0: 0.500000") != std::string::npos &&
              rz.find("1: 0.500000") != std::string::npos,
          "RZ preserves probabilities and normalization");
  require(run("qubit q0\nH q0\nmeasure q0\n") ==
              run("qubit q0\nH q0\nmeasure q0\n"),
          "fixed-seed measurement output is repeatable");

  // Phase 5 C6: Full-register basis-state bitstrings q(n-1)...q0 (q0 least significant, rightmost)
  {
    Simulator sim(2);
    require(sim.bitString(0) == "00", "state |0> is '00'");
    require(sim.bitString(1) == "01", "state |1> (q0=1, q1=0) is '01'");
    require(sim.bitString(2) == "10", "state |2> (q0=0, q1=1) is '10'");
    require(sim.bitString(3) == "11", "state |3> (q0=1, q1=1) is '11'");

    Simulator sim3(3);
    require(sim3.bitString(1) == "001", "q0=1 is '001'");
    require(sim3.bitString(2) == "010", "q1=1 is '010'");
    require(sim3.bitString(4) == "100", "q2=1 is '100'");
    require(sim3.bitString(3) == "011", "q0=1, q1=1 is '011'");
  }

  // Phase 5 C7: Operand-order measurement result mapping
  {
    // Prepare state |q2 q1 q0> = |0 1 1> (state index 3)
    // q0=1 (X q0), q1=1 (X q1), q2=0
    std::string out02 = run("qubit q0, q1, q2\nX q0\nX q1\nmeasure q0, q2\n");
    require(out02.find("Result: 10") != std::string::npos,
            "measure q0, q2 on |011> yields '10' (q0=1, q2=0)");

    std::string out20 = run("qubit q0, q1, q2\nX q0\nX q1\nmeasure q2, q0\n");
    require(out20.find("Result: 01") != std::string::npos,
            "measure q2, q0 on |011> yields '01' (q2=0, q0=1)");
  }

  // Phase 5 C8: Pre-measurement probability labeling
  {
    std::string out = run("qubit q0\nH q0\nmeasure q0\n");
    require(out.find("Pre-measurement probabilities:\n") != std::string::npos,
            "output explicitly labels pre-measurement probabilities");
  }

  // Phase 5 C2: Repeat measurement invariant without intervening gates
  {
    std::string out = run("qubit q0\nH q0\nmeasure q0\nmeasure q0\n");
    std::size_t firstResPos = out.find("Result: ");
    require(firstResPos != std::string::npos, "first measurement result found");
    std::size_t secondResPos = out.find("Result: ", firstResPos + 8);
    require(secondResPos != std::string::npos, "second measurement result found");
    std::string firstResult = out.substr(firstResPos + 8, 1);
    std::string secondResult = out.substr(secondResPos + 8, 1);
    require(firstResult == secondResult,
            "repeated measurement of q0 yields identical outcome (state collapsed)");

    // Inspect second pre-measurement probabilities: collapsed to observed outcome with probability 1.0
    std::size_t secondPrePos = out.find("Pre-measurement probabilities:\n", firstResPos);
    require(secondPrePos != std::string::npos, "second pre-measurement section found");
    std::string secondProbs = out.substr(secondPrePos, secondResPos - secondPrePos);
    require(secondProbs.find(firstResult + ": 1.000000") != std::string::npos,
            "second measurement observes collapsed state with probability 1.0");
  }

  // Phase 5 C3: Intervening gate acts on collapsed state
  {
    // Prepare |+>, measure q0 -> collapsed to outcome B in {0, 1}.
    // Then apply X q0 -> flips outcome to (1 - B).
    // Second measure q0 must return flipped outcome with probability 1.0.
    std::string out = run("qubit q0\nH q0\nmeasure q0\nX q0\nmeasure q0\n");
    std::size_t r1Pos = out.find("Result: ");
    std::size_t r2Pos = out.find("Result: ", r1Pos + 8);
    require(r1Pos != std::string::npos && r2Pos != std::string::npos, "both results found");
    char r1 = out[r1Pos + 8];
    char r2 = out[r2Pos + 8];
    require((r1 == '0' && r2 == '1') || (r1 == '1' && r2 == '0'),
            "X gate flips collapsed state outcome");
  }

  // Phase 5 C4: Partial measurement and correlation preservation on entangled Bell state
  {
    // Bell state (|00> + |11>)/sqrt(2). Measure only q0.
    // If q0 collapses to 0, state collapses to |00> (q1 is 0 with probability 1.0).
    // If q0 collapses to 1, state collapses to |11> (q1 is 1 with probability 1.0).
    // Then measure q1: its outcome must match q0's outcome!
    std::string out = run("qubit q0, q1\nH q0\nCNOT q0, q1\nmeasure q0\nmeasure q1\n");
    std::size_t r1Pos = out.find("Result: ");
    std::size_t r2Pos = out.find("Result: ", r1Pos + 8);
    require(r1Pos != std::string::npos && r2Pos != std::string::npos, "both results found");
    char r0 = out[r1Pos + 8];
    char r1 = out[r2Pos + 8];
    require(r0 == r1, "measuring q0 collapses entangled partner q1 to matching outcome");
  }

  // Phase 5 C5: Joint multi-qubit collapse
  {
    // In 3-qubit state (|000> + |111>)/sqrt(2), measure q0, q1 simultaneously.
    // Outcome must be either 00 or 11; state vector must collapse completely to |000> or |111>.
    std::string out = run("qubit q0, q1, q2\nH q0\nCNOT q0, q1\nCNOT q0, q2\nmeasure q0, q1\nmeasure q2\n");
    std::size_t r1Pos = out.find("Result: ");
    std::size_t r2Pos = out.find("Result: ", r1Pos + 8);
    require(r1Pos != std::string::npos && r2Pos != std::string::npos, "both results found");
    std::string jointResult = out.substr(r1Pos + 8, 2);
    char q2Result = out[r2Pos + 8];
    require((jointResult == "00" && q2Result == '0') ||
            (jointResult == "11" && q2Result == '1'),
            "joint multi-qubit measurement collapses 3-qubit GHZ state consistently");
  }

  // Phase 5 strengthening: Deterministic partial measurement preserving nontrivial superposition on unmeasured subsystem
  const double invSqrt2 = 1.0 / std::sqrt(2.0);

  // Case A: prepare q1 = |0> and q0 = |+>, measure q1.
  // In basis |q1 q0>: state is |0>_1 (x) |+>_0 = (|00> + |01>)/sqrt(2).
  // Outcome on q1 must be 0 with probability 1.0.
  // Surviving amplitudes: state[0] (|00>) == 1/sqrt(2), state[1] (|01>) == 1/sqrt(2),
  //                      state[2] (|10>) == 0.0,       state[3] (|11>) == 0.0.
  {
    Simulator simA(2);
    std::vector<Instruction> instrsA;
    instrsA.push_back(makeInstruction(InstructionType::ApplyGate, "H", {0})); // q0 = |+>, q1 = |0>
    instrsA.push_back(makeInstruction(InstructionType::Measure, "measure", {1})); // measure q1

    std::ostringstream dummyA;
    std::streambuf *origBufA = std::cout.rdbuf(dummyA.rdbuf());
    simA.run(instrsA);
    std::cout.rdbuf(origBufA);

    const auto &vecA = simA.stateVector();
    require(vecA.size() == 4, "Case A: 2-qubit statevector has 4 amplitudes");

    double normA = 0.0;
    for (const auto &amp : vecA) {
      normA += std::norm(amp);
    }
    require(std::abs(normA - 1.0) < 1e-12, "Case A: projected state vector has unit norm");

    require(std::abs(vecA[0] - Complex(invSqrt2, 0.0)) < 1e-12, "Case A: |00> amplitude is 1/sqrt(2)");
    require(std::abs(vecA[1] - Complex(invSqrt2, 0.0)) < 1e-12, "Case A: |01> amplitude is 1/sqrt(2)");
    require(std::norm(vecA[2]) < 1e-12, "Case A: |10> amplitude is 0.0 outside measured subspace");
    require(std::norm(vecA[3]) < 1e-12, "Case A: |11> amplitude is 0.0 outside measured subspace");
  }

  // Case B: prepare q1 = |1> (via X q1) and q0 = |+>, measure q1.
  // In basis |q1 q0>: state is |1>_1 (x) |+>_0 = (|10> + |11>)/sqrt(2).
  // Outcome on q1 must be 1 with probability 1.0.
  // Surviving amplitudes: state[0] (|00>) == 0.0,       state[1] (|01>) == 0.0,
  //                      state[2] (|10>) == 1/sqrt(2), state[3] (|11>) == 1/sqrt(2).
  {
    Simulator simB(2);
    std::vector<Instruction> instrsB;
    instrsB.push_back(makeInstruction(InstructionType::ApplyGate, "H", {0})); // q0 = |+>
    instrsB.push_back(makeInstruction(InstructionType::ApplyGate, "X", {1})); // q1 = |1>
    instrsB.push_back(makeInstruction(InstructionType::Measure, "measure", {1})); // measure q1

    std::ostringstream dummyB;
    std::streambuf *origBufB = std::cout.rdbuf(dummyB.rdbuf());
    simB.run(instrsB);
    std::cout.rdbuf(origBufB);

    const auto &vecB = simB.stateVector();
    require(vecB.size() == 4, "Case B: 2-qubit statevector has 4 amplitudes");

    double normB = 0.0;
    for (const auto &amp : vecB) {
      normB += std::norm(amp);
    }
    require(std::abs(normB - 1.0) < 1e-12, "Case B: projected state vector has unit norm");

    require(std::norm(vecB[0]) < 1e-12, "Case B: |00> amplitude is 0.0 outside measured subspace");
    require(std::norm(vecB[1]) < 1e-12, "Case B: |01> amplitude is 0.0 outside measured subspace");
    require(std::abs(vecB[2] - Complex(invSqrt2, 0.0)) < 1e-12, "Case B: |10> amplitude is 1/sqrt(2)");
    require(std::abs(vecB[3] - Complex(invSqrt2, 0.0)) < 1e-12, "Case B: |11> amplitude is 1/sqrt(2)");
  }
}

void testCapacity() {
  for (int count : {-1, 25}) {
    bool rejected = false;
    try {
      Simulator simulator(count);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    require(rejected, "negative and over-limit qubit counts rejected");
  }

  const auto rejectsInstruction = [](const Instruction &instruction) {
    try {
      Simulator simulator(2);
      simulator.run({instruction});
    } catch (const std::invalid_argument &) {
      return true;
    }
    return false;
  };
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "X", {-1}, 0.0}),
          "negative single-qubit index rejected");
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "X", {2}, 0.0}),
          "index equal to qubit count rejected");
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "CNOT", {-1, 1}, 0.0}),
          "negative CNOT control index rejected");
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "CNOT", {0, 2}, 0.0}),
          "out-of-range CNOT target index rejected");
  require(rejectsInstruction(
              {InstructionType::Measure, "measure", {0, -1}, 0.0}),
          "invalid measurement index rejected");
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "CNOT", {0}, 0.0}),
          "missing CNOT operand rejected");
  require(rejectsInstruction(
              {InstructionType::ApplyGate, "UNKNOWN", {0}, 0.0}),
          "unknown gate instruction rejected");
  Simulator smallCircuit(2);
  std::vector<Instruction> instructions;
  instructions.push_back({InstructionType::ApplyGate, "X", {1}, 0.0});
  instructions.push_back({InstructionType::Measure, "measure", {0, 1}, 0.0});
  std::ostringstream output;
  std::streambuf *original = std::cout.rdbuf(output.rdbuf());
  smallCircuit.run(instructions);
  std::cout.rdbuf(original);
  require(output.str().find("10: 1.000000") != std::string::npos,
          "supported small circuit executes normally with q1q0 basis order");
}

} // namespace test
