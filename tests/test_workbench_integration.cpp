#include "test_support.h"
#include "compiler/compiler_facade.h"
#include "simulator/simulation_session.h"

#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace test {

void testWorkbenchIntegration() {
  // 1. Preset 0: Bell State compilation and interactive stepping synchronization
  {
    const std::string bellDsl =
        "# Quantum Bell State Preparation\n"
        "qubit q0, q1\n"
        "H q0\n"
        "CNOT q0, q1\n"
        "measure q0, q1\n";

    compiler::CompilerOptions options;
    options.optimize = true;
    const compiler::CompilationResult compilation = compiler::compileProgram(bellDsl, options);

    require(compilation.success, "Bell state preset compiles successfully");
    require(compilation.qubitCount == 2, "Bell state has 2 declared qubits");
    require(compilation.instructions.size() == 3, "Bell state has 3 instructions (H, CNOT, measure)");
    require(compilation.diagnostics.empty(), "zero diagnostics on valid Bell state");

    SimulationSession session(compilation.instructions, static_cast<int>(compilation.qubitCount));
    require(session.totalSteps() == 3, "total session steps is 3");
    require(session.currentStep() == 0, "session starts at step 0");
    require(!session.isComplete(), "session not complete at step 0");

    // Initial state |00> with probability 1.0
    require(std::abs(session.currentState()[0] - Complex(1.0, 0.0)) < 1e-12, "initial state is |00>");
    require(std::abs(session.stateNorm() - 1.0) < 1e-12, "initial state norm is 1.0");

    // Step 1: H on q0 -> (|00> + |01>) / sqrt(2)
    require(session.stepForward(), "stepForward 1 succeeds");
    require(session.currentStep() == 1, "step is 1 after H");
    const double invSqrt2 = 1.0 / std::sqrt(2.0);
    require(std::abs(session.currentState()[0] - Complex(invSqrt2, 0.0)) < 1e-12 &&
            std::abs(session.currentState()[1] - Complex(invSqrt2, 0.0)) < 1e-12,
            "step 1 state matches (|00> + |01>)/sqrt(2)");
    require(session.basisProbabilities()[0] > 0.499999 && session.basisProbabilities()[1] > 0.499999,
            "step 1 basis probabilities are 0.5 for |00> and |01>");

    // Step 2: CNOT q0, q1 -> Entangled Bell state (|00> + |11>) / sqrt(2)
    require(session.stepForward(), "stepForward 2 succeeds");
    require(session.currentStep() == 2, "step is 2 after CNOT");
    require(std::abs(session.currentState()[0] - Complex(invSqrt2, 0.0)) < 1e-12 &&
            std::abs(session.currentState()[3] - Complex(invSqrt2, 0.0)) < 1e-12,
            "step 2 state matches Bell state (|00> + |11>)/sqrt(2)");
    require(std::abs(session.currentState()[1]) < 1e-12 && std::abs(session.currentState()[2]) < 1e-12,
            "orthogonal basis states |01> and |10> have zero amplitude");

    // Step 3: Projective measurement of q0 and q1
    require(session.stepForward(), "stepForward 3 succeeds");
    require(session.currentStep() == 3, "step is 3 after measure");
    require(session.isComplete(), "session is complete at step 3");
    require(session.lastMeasurementResult().has_value(), "measurement result recorded");

    const std::string outcome = *session.lastMeasurementResult();
    require(outcome == "00" || outcome == "11", "Bell measurement outcome is '00' or '11'");

    // Collapsed state matches outcome
    if (outcome == "00") {
      require(std::abs(session.currentState()[0] - Complex(1.0, 0.0)) < 1e-12,
              "outcome '00' collapses register to |00>");
      require(std::abs(session.currentState()[3]) < 1e-12, "|11> amplitude is 0");
    } else {
      require(std::abs(session.currentState()[3] - Complex(1.0, 0.0)) < 1e-12,
              "outcome '11' collapses register to |11>");
      require(std::abs(session.currentState()[0]) < 1e-12, "|00> amplitude is 0");
    }

    // Step backward restores superposition
    require(session.stepBackward(), "stepBackward restores step 2");
    require(session.currentStep() == 2, "step is 2 after stepBackward");
    require(std::abs(session.basisProbabilities()[0] - 0.5) < 1e-12 &&
            std::abs(session.basisProbabilities()[3] - 0.5) < 1e-12,
            "superposition restored on backward step");

    // Step forward restores identical outcome
    require(session.stepForward(), "stepForward re-advances to step 3");
    require(*session.lastMeasurementResult() == outcome,
            "measurement outcome stably preserved across stepping");

    // Reset rewinds to initial state
    session.reset();
    require(session.currentStep() == 0, "currentStep is 0 after reset");
    require(std::abs(session.currentState()[0] - Complex(1.0, 0.0)) < 1e-12,
            "initial state restored after reset");
  }

  // 2. Preset 1: Superposition and RZ Phase Rotation (Analytical derivation check)
  {
    const std::string phaseDsl =
        "qubit q0\n"
        "H q0\n"
        "let theta = 3.14159265358979323846 / 4.0\n"
        "RZ(theta) q0\n"
        "H q0\n";

    const compiler::CompilationResult compilation = compiler::compileProgram(phaseDsl);
    require(compilation.success, "phase rotation program compiles");
    require(compilation.instructions.size() == 3, "3 instructions generated");

    SimulationSession session(compilation.instructions, 1);
    session.seek(3);
    require(session.isComplete(), "session completed 3 steps");

    // Analytical calculation for H * RZ(pi/4) * H |0>:
    // P(|0>) = cos^2(pi/8) = (1 + 1/sqrt(2))/2 ~= 0.85355339059
    // P(|1>) = sin^2(pi/8) = (1 - 1/sqrt(2))/2 ~= 0.14644660940
    const double expectedP0 = (1.0 + 1.0 / std::sqrt(2.0)) / 2.0;
    const double expectedP1 = (1.0 - 1.0 / std::sqrt(2.0)) / 2.0;

    const auto probs = session.basisProbabilities();
    require(std::abs(probs[0] - expectedP0) < 1e-10, "P(|0>) matches exact analytical cos^2(pi/8)");
    require(std::abs(probs[1] - expectedP1) < 1e-10, "P(|1>) matches exact analytical sin^2(pi/8)");
    require(std::abs(session.stateNorm() - 1.0) < 1e-12, "state norm is 1.0");
  }

  // 3. Preset 2: 3-Qubit GHZ State Entanglement
  {
    const std::string ghzDsl =
        "qubit q0, q1, q2\n"
        "H q0\n"
        "CNOT q0, q1\n"
        "CNOT q1, q2\n"
        "measure q0, q1, q2\n";

    const compiler::CompilationResult compilation = compiler::compileProgram(ghzDsl);
    require(compilation.success, "GHZ circuit compiles");
    require(compilation.qubitCount == 3, "3 qubits in GHZ circuit");

    SimulationSession session(compilation.instructions, 3);
    session.seek(3); // prior to measurement

    // State at step 3: (|000> + |111>) / sqrt(2) -> indices 0 and 7
    const auto probs = session.basisProbabilities();
    require(std::abs(probs[0] - 0.5) < 1e-12, "GHZ state P(|000>) = 0.5");
    require(std::abs(probs[7] - 0.5) < 1e-12, "GHZ state P(|111>) = 0.5");
    for (int i = 1; i < 7; ++i) {
      require(probs[static_cast<std::size_t>(i)] < 1e-12, "non-GHZ states have zero probability");
    }

    // Step 4: measurement collapses to '000' or '111'
    session.stepForward();
    require(session.isComplete(), "session complete after GHZ measurement");
    const std::string ghzOutcome = *session.lastMeasurementResult();
    require(ghzOutcome == "000" || ghzOutcome == "111",
            "GHZ joint measurement outcome is strictly '000' or '111'");
  }

  // 4. Invalid source recovery
  {
    const std::string invalidDsl =
        "qubit q0\n"
        "H missing_qubit\n"
        "measure q0, q0\n"
        "let bad = 1.0 / 0.0\n";

    const compiler::CompilationResult badResult = compiler::compileProgram(invalidDsl);
    require(!badResult.success, "invalid DSL rejected with success = false");
    require(badResult.diagnostics.size() >= 3, "all 3 diagnostics reported");
    require(badResult.diagnostics[0].line == 2, "missing qubit diagnostic on line 2");
    require(badResult.diagnostics[1].line == 3, "duplicate measure diagnostic on line 3");
    require(badResult.diagnostics[2].line == 4, "division by zero diagnostic on line 4");

    // Clean recovery on subsequent compilation of valid source
    const std::string recoveryDsl = "qubit q0\nH q0\n";
    const compiler::CompilationResult goodResult = compiler::compileProgram(recoveryDsl);
    require(goodResult.success, "compiler recovers cleanly from error state");
    require(goodResult.diagnostics.empty(), "clean compilation has zero diagnostics");
    require(goodResult.qubitCount == 1, "1 qubit declared");
    require(goodResult.instructions.size() == 1, "1 instruction emitted");
  }

  // 5. Raw vs Optimized circuit comparison
  {
    const std::string redundantDsl =
        "qubit q0\n"
        "H q0\n"
        "H q0\n"
        "X q0\n";

    compiler::CompilerOptions unopt;
    unopt.optimize = false;
    const auto unoptResult = compiler::compileProgram(redundantDsl, unopt);
    require(unoptResult.rawInstructions.size() == 3, "raw instruction count is 3");
    require(unoptResult.instructions.size() == 3, "unoptimized instruction count is 3");

    compiler::CompilerOptions opt;
    opt.optimize = true;
    const auto optResult = compiler::compileProgram(redundantDsl, opt);
    require(optResult.rawInstructions.size() == 3, "raw instructions preserved in opt mode");
    require(optResult.instructions.size() == 1, "optimized instructions reduced adjacent H-H");
    require(optResult.instructions[0].gateName == "X", "remaining instruction is X");

    // Verify both produce identical state |1>
    SimulationSession unoptSession(unoptResult.instructions, 1);
    unoptSession.seek(3);
    SimulationSession optSession(optResult.instructions, 1);
    optSession.seek(1);

    requireSameState(unoptSession.currentState(), optSession.currentState(),
                     "unoptimized and optimized circuits yield identical final statevector");
    require(std::abs(optSession.currentState()[1] - Complex(1.0, 0.0)) < 1e-12,
            "final state is |1>");
  }
}

} // namespace test
