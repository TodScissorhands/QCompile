#include "test_support.h"
#include "chemistry/h2_energy.h"
#include "simulator/simulation_session.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace test {

void testSimulationSession() {
  // 1. Empty circuit behavior
  {
    SimulationSession session({}, 2);
    require(session.totalSteps() == 0, "empty circuit totalSteps is 0");
    require(session.currentStep() == 0, "empty circuit currentStep is 0");
    require(session.isComplete(), "empty circuit is immediately complete");
    require(!session.stepForward(), "stepForward() returns false on empty circuit");
    require(!session.stepBackward(), "stepBackward() returns false on empty circuit");
    require(session.qubitCount() == 2, "qubitCount is 2");

    // Initial state must be |00> with probability 1.0
    const auto &state = session.currentState();
    require(state.size() == 4, "2-qubit state dimension is 4");
    require(std::abs(state[0] - Complex(1.0, 0.0)) < 1e-12, "initial state is |00>");
    require(std::abs(state[1]) < 1e-12 && std::abs(state[2]) < 1e-12 && std::abs(state[3]) < 1e-12,
            "non-zero basis states have zero amplitude");
    require(std::abs(session.stateNorm() - 1.0) < 1e-12, "initial state norm is 1.0");

    const auto probs = session.basisProbabilities();
    require(probs.size() == 4 && std::abs(probs[0] - 1.0) < 1e-12,
            "basis probabilities have index 0 at 1.0");

    // Seek to 0 is valid on empty circuit
    session.seek(0);
    require(session.currentStep() == 0, "seek(0) succeeds on empty circuit");

    // Seek to 1 must throw out_of_range
    bool seekRejected = false;
    try {
      session.seek(1);
    } catch (const std::out_of_range &) {
      seekRejected = true;
    }
    require(seekRejected, "seek past totalSteps throws std::out_of_range");

    // Reset works on empty circuit
    session.reset();
    require(session.currentStep() == 0, "reset on empty circuit leaves currentStep 0");
  }

  // 2. Single-qubit gates and stepwise forward/backward reversibility
  {
    std::vector<Instruction> circuit{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "X", {1}),
        makeInstruction(InstructionType::ApplyGate, "Z", {0})
    };

    SimulationSession session(circuit, 2);
    require(session.totalSteps() == 3, "total steps is 3");
    require(session.currentStep() == 0, "initial step is 0");
    require(!session.isComplete(), "not complete at step 0");

    // Step 1: H on q0 -> (|00> + |01>) / sqrt(2)
    require(session.stepForward(), "stepForward 1 succeeds");
    require(session.currentStep() == 1, "current step is 1");
    const double invSqrt2 = 1.0 / std::sqrt(2.0);
    require(std::abs(session.currentState()[0] - Complex(invSqrt2, 0.0)) < 1e-12 &&
            std::abs(session.currentState()[1] - Complex(invSqrt2, 0.0)) < 1e-12 &&
            std::abs(session.currentState()[2]) < 1e-12 &&
            std::abs(session.currentState()[3]) < 1e-12,
            "step 1 state matches (|00> + |01>)/sqrt(2)");
    require(std::abs(session.stateNorm() - 1.0) < 1e-12, "norm conserved at step 1");

    // Step 2: X on q1 -> (|10> + |11>) / sqrt(2) (indices 2 and 3)
    require(session.stepForward(), "stepForward 2 succeeds");
    require(session.currentStep() == 2, "current step is 2");
    require(std::abs(session.currentState()[2] - Complex(invSqrt2, 0.0)) < 1e-12 &&
            std::abs(session.currentState()[3] - Complex(invSqrt2, 0.0)) < 1e-12 &&
            std::abs(session.currentState()[0]) < 1e-12 &&
            std::abs(session.currentState()[1]) < 1e-12,
            "step 2 state matches (|10> + |11>)/sqrt(2)");

    // Step backward to step 1
    require(session.stepBackward(), "stepBackward 1 succeeds");
    require(session.currentStep() == 1, "step restored to 1");
    require(std::abs(session.currentState()[0] - Complex(invSqrt2, 0.0)) < 1e-12 &&
            std::abs(session.currentState()[1] - Complex(invSqrt2, 0.0)) < 1e-12,
            "step 1 state exactly restored after backward step");

    // Step backward to step 0
    require(session.stepBackward(), "stepBackward 2 succeeds");
    require(session.currentStep() == 0, "step restored to 0");
    require(std::abs(session.currentState()[0] - Complex(1.0, 0.0)) < 1e-12,
            "initial |00> restored after rewinding to step 0");
    require(!session.stepBackward(), "stepBackward at 0 returns false");

    // Step all the way to completion (step 3)
    require(session.stepForward(), "stepForward 1 re-executed");
    require(session.stepForward(), "stepForward 2 re-executed");
    require(session.stepForward(), "stepForward 3 executed");
    require(session.isComplete(), "isComplete() true after step 3");
    require(!session.stepForward(), "stepForward returns false when complete");

    // Step 3 applies Z to q0 on (|10> + |11>)/sqrt(2):
    // q0 is bit 0. State |10> has q0=0 (phase +1), |11> has q0=1 (phase -1).
    require(std::abs(session.currentState()[2] - Complex(invSqrt2, 0.0)) < 1e-12 &&
            std::abs(session.currentState()[3] - Complex(-invSqrt2, 0.0)) < 1e-12,
            "Z gate applies -1 phase to |11>");
  }

  // 3. Bell-state preparation and seek reversibility
  {
    std::vector<Instruction> bellCircuit{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1})
    };

    SimulationSession session(bellCircuit, 2);

    // Seek directly to step 2 (complete Bell state)
    session.seek(2);
    require(session.currentStep() == 2, "seek(2) advances to step 2");
    require(session.isComplete(), "session complete after seek(2)");

    const auto probs = session.basisProbabilities();
    require(std::abs(probs[0] - 0.5) < 1e-12, "Bell state P(|00>) = 0.5");
    require(std::abs(probs[3] - 0.5) < 1e-12, "Bell state P(|11>) = 0.5");
    require(probs[1] < 1e-12 && probs[2] < 1e-12, "P(|01>) = P(|10>) = 0");

    // Seek back to 1
    session.seek(1);
    require(session.currentStep() == 1, "seek(1) rewinds to step 1");
    require(std::abs(session.basisProbabilities()[0] - 0.5) < 1e-12 &&
            std::abs(session.basisProbabilities()[1] - 0.5) < 1e-12,
            "step 1 has equal superposition on q0");

    // Seek to 0
    session.seek(0);
    require(session.currentStep() == 0, "seek(0) rewinds to initial state");
    require(std::abs(session.basisProbabilities()[0] - 1.0) < 1e-12, "P(|00>) = 1.0 at step 0");

    // Seek directly back to 2
    session.seek(2);
    require(session.currentStep() == 2, "seek(2) restores completed Bell state");
    require(std::abs(session.basisProbabilities()[0] - 0.5) < 1e-12 &&
            std::abs(session.basisProbabilities()[3] - 0.5) < 1e-12,
            "Bell state probabilities intact after seeking");
  }

  // 4. Measurement collapse and deterministic history preservation
  {
    std::vector<Instruction> measCircuit{
        makeInstruction(InstructionType::ApplyGate, "H", {0}),
        makeInstruction(InstructionType::ApplyGate, "CNOT", {0, 1}),
        makeInstruction(InstructionType::Measure, "measure", {0})
    };

    SimulationSession session(measCircuit, 2);

    // Step through H and CNOT to create Bell state
    session.stepForward(); // step 1: H
    session.stepForward(); // step 2: CNOT
    require(!session.lastMeasurementResult().has_value(),
            "no measurement result prior to measure instruction");

    // Step 3: Measure q0
    session.stepForward(); // step 3: measure
    require(session.isComplete(), "circuit complete after measure");
    require(session.lastMeasurementResult().has_value(),
            "measurement result available after measurement step");

    const std::string recordedResult = *session.lastMeasurementResult();
    require(recordedResult == "0" || recordedResult == "1",
            "measured single qubit outcome is '0' or '1'");

    // Inspect collapsed statevector
    const auto collapsedState = session.currentState();
    if (recordedResult == "0") {
      require(std::abs(collapsedState[0] - Complex(1.0, 0.0)) < 1e-12,
              "outcome 0 projects Bell state to |00>");
      require(std::abs(collapsedState[3]) < 1e-12, "|11> amplitude is 0 when q0 collapsed to 0");
    } else {
      require(std::abs(collapsedState[3] - Complex(1.0, 0.0)) < 1e-12,
              "outcome 1 projects Bell state to |11>");
      require(std::abs(collapsedState[0]) < 1e-12, "|00> amplitude is 0 when q0 collapsed to 1");
    }

    // Step backward to step 2 (prior to measurement)
    session.stepBackward();
    require(session.currentStep() == 2, "rewound to step 2");
    require(std::abs(session.basisProbabilities()[0] - 0.5) < 1e-12 &&
            std::abs(session.basisProbabilities()[3] - 0.5) < 1e-12,
            "unmeasured superposition restored before measurement");

    // Step forward back to step 3: MUST reproduce exact recorded measurement!
    session.stepForward();
    require(session.currentStep() == 3, "advanced back to step 3");
    require(session.lastMeasurementResult().has_value(), "measurement result present");
    require(*session.lastMeasurementResult() == recordedResult,
            "measurement outcome preserved stably through backward/forward stepping");

    // Seek back to 0, then seek to 3: outcome MUST remain identical
    session.seek(0);
    session.seek(3);
    require(*session.lastMeasurementResult() == recordedResult,
            "measurement outcome preserved across seek(0) -> seek(3)");

    // Reset without clearHistory retains the outcome
    session.reset(/*clearHistory=*/false);
    require(session.currentStep() == 0, "reset rewinds to step 0");
    session.stepForward();
    session.stepForward();
    session.stepForward();
    require(*session.lastMeasurementResult() == recordedResult,
            "measurement outcome preserved across reset(false)");

    // Reset with clearHistory = true discards snapshots beyond step 0
    session.reset(/*clearHistory=*/true);
    require(session.currentStep() == 0, "reset(true) rewinds to step 0");
    require(!session.isComplete(), "session not complete at step 0");
    require(!session.lastMeasurementResult().has_value(), "no measurement result at step 0");
    require(!session.measurementResultAt(1).has_value(), "step 1 measurement discarded");
    require(!session.measurementResultAt(2).has_value(), "step 2 measurement discarded");
    require(!session.measurementResultAt(3).has_value(), "step 3 measurement discarded");

    // Subsequent execution creates a fresh history without assuming a specific outcome
    session.stepForward(); // step 1: H
    session.stepForward(); // step 2: CNOT
    session.stepForward(); // step 3: measure
    require(session.currentStep() == 3, "advanced through fresh run to step 3");
    require(session.isComplete(), "fresh run complete at step 3");
    require(session.lastMeasurementResult().has_value(), "fresh measurement result obtained");

    const std::string freshResult = *session.lastMeasurementResult();
    require(freshResult == "0" || freshResult == "1", "fresh measurement result is '0' or '1'");

    const auto freshState = session.currentState();
    if (freshResult == "0") {
      require(std::abs(freshState[0] - Complex(1.0, 0.0)) < 1e-12, "fresh run 0 projects to |00>");
      require(std::abs(freshState[3]) < 1e-12, "|11> amplitude is 0 when fresh outcome is 0");
    } else {
      require(std::abs(freshState[3] - Complex(1.0, 0.0)) < 1e-12, "fresh run 1 projects to |11>");
      require(std::abs(freshState[0]) < 1e-12, "|00> amplitude is 0 when fresh outcome is 1");
    }

    // Verify this fresh history is preserved across subsequent backward/forward steps
    session.stepBackward();
    session.stepForward();
    require(*session.lastMeasurementResult() == freshResult,
            "fresh history is preserved across subsequent backward/forward");
  }

  // 5. Invalid inputs and resource limits
  {
    // Out of range qubit counts
    bool negativeQubitsRejected = false;
    try {
      SimulationSession session({}, -1);
    } catch (const std::invalid_argument &) {
      negativeQubitsRejected = true;
    }
    require(negativeQubitsRejected, "negative qubit count rejected");

    bool excessiveQubitsRejected = false;
    try {
      SimulationSession session({}, 25);
    } catch (const std::invalid_argument &) {
      excessiveQubitsRejected = true;
    }
    require(excessiveQubitsRejected, "qubit count > 24 rejected");

    // Instruction with out of range qubit index rejected at construction
    bool invalidQubitRejected = false;
    try {
      SimulationSession session({makeInstruction(InstructionType::ApplyGate, "H", {3})}, 2);
    } catch (const std::invalid_argument &) {
      invalidQubitRejected = true;
    }
    require(invalidQubitRejected, "instruction with out-of-range qubit index rejected");

    // Out of range seek
    SimulationSession validSession({makeInstruction(InstructionType::ApplyGate, "H", {0})}, 1);
    bool outOfRangeSeekRejected = false;
    try {
      validSession.seek(2);
    } catch (const std::out_of_range &) {
      outOfRangeSeekRejected = true;
    }
    require(outOfRangeSeekRejected, "seek beyond total steps throws std::out_of_range");

    // 24-qubit boundary memory accounting:
    // A single 24-qubit statevector is 256 MiB (2^24 * 16 bytes).
    // The session retains 1 active simulator buffer + step 0 snapshot = 512 MiB at construction.
    // Advancing to step 1 would require a 3rd buffer (768 MiB > 512 MiB limit),
    // so stepForward() must reject before allocating and preserve step 0.
    {
      SimulationSession session24({makeInstruction(InstructionType::ApplyGate, "X", {0})}, 24);
      require(session24.currentStep() == 0, "24-qubit session constructs successfully at step 0");
      bool stepRejected = false;
      try {
        session24.stepForward();
      } catch (const std::length_error &) {
        stepRejected = true;
      }
      require(stepRejected, "stepForward rejected when total state buffers exceed 512 MiB");
      require(session24.currentStep() == 0, "session remains at step 0 after memory limit rejection");
    }
  }

  // 6. Integration with H2 trial instructions
  {
    const auto h2Instructions = h2TrialInstructions(0.5);
    require(h2Instructions.size() == 8, "H2 trial sequence has 8 instructions");

    SimulationSession session(h2Instructions, 2);
    require(session.totalSteps() == 8, "session totalSteps is 8");

    // Step through all 8 gates and check norm conservation at every step
    for (std::size_t step = 0; step < 8; ++step) {
      require(session.stepForward(), "stepForward succeeds on H2 gate");
      require(std::abs(session.stateNorm() - 1.0) < 1e-12, "norm strictly conserved on H2 gate");
    }
    require(session.isComplete(), "H2 session complete after 8 steps");

    // Compare completed state with standard Simulator::run
    Simulator directSimulator(2);
    directSimulator.run(h2Instructions);
    requireSameState(session.currentState(), directSimulator.stateVector(),
                     "SimulationSession final state matches direct Simulator run");
  }
}

} // namespace test
