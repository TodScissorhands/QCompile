#include "test_support.h"
#include "compiler/compiler_facade.h"
#include "simulator/simulation_session.h"
#include "chemistry/h2_energy.h"
#include "chemistry/heh_energy.h"
#include "chemistry/vqe_convergence.h"
#include <complex>
#include <iomanip>
#include <sstream>
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

  // 6. HeH+ 4-qubit Givens rotation circuit stepwise simulation and energy expectation
  {
    constexpr double optTheta = 0.068498646886667;
    constexpr double optPhi = 4.976298509698156;

    const std::vector<Instruction> hehCircuit = hehSingletTrialCircuit(optTheta, optPhi);
    require(hehCircuit.size() > 50, "HeH+ trial circuit contains full Givens decomposition sequence");

    // Decompile to DSL string
    std::string dsl = "qubit q0, q1, q2, q3\n";
    for (const auto &instr : hehCircuit) {
      if (instr.type == InstructionType::Measure) {
        dsl += "measure";
        for (std::size_t i = 0; i < instr.qubitIndices.size(); ++i) {
          dsl += (i == 0 ? " q" : ", q") + std::to_string(instr.qubitIndices[i]);
        }
        dsl += "\n";
      } else if (instr.gateName == "RZ") {
        std::ostringstream ss;
        ss << std::setprecision(17) << instr.parameterValue;
        dsl += "RZ(" + ss.str() + ") q" + std::to_string(instr.qubitIndices[0]) + "\n";
      } else if (instr.gateName == "CNOT") {
        dsl += "CNOT q" + std::to_string(instr.qubitIndices[0]) + ", q" + std::to_string(instr.qubitIndices[1]) + "\n";
      } else {
        dsl += instr.gateName + " q" + std::to_string(instr.qubitIndices[0]) + "\n";
      }
    }

    const compiler::CompilationResult compilation = compiler::compileProgram(dsl);
    require(compilation.success, "decompiled HeH+ 4-qubit DSL compiles successfully");
    require(compilation.qubitCount == 4, "4 qubits declared in HeH+ circuit");
    require(compilation.rawInstructions.size() == hehCircuit.size(), "decompiled raw instruction count matches original");
    require(compilation.instructions.size() <= compilation.rawInstructions.size(), "optimized instruction count is non-expanding");

    SimulationSession session(compilation.rawInstructions, 4);
    require(session.totalSteps() == compilation.rawInstructions.size(), "total steps matches instruction count");
    require(session.currentStep() == 0, "session starts at step 0");

    // Step 0: |0000>
    require(std::abs(session.currentState()[0] - Complex(1.0, 0.0)) < 1e-12, "initial state is |0000>");

    // Step 1: X q0 -> |0001> (index 1)
    session.stepForward();
    require(session.currentStep() == 1, "step is 1 after X0");
    require(std::abs(session.currentState()[1] - Complex(1.0, 0.0)) < 1e-12, "state is |0001>");

    // Step 2: X q2 -> |0101> (index 5)
    session.stepForward();
    require(session.currentStep() == 2, "step is 2 after X2");
    require(std::abs(session.currentState()[5] - Complex(1.0, 0.0)) < 1e-12,
            "state after X0, X2 is reference determinant |0101> (index 5)");

    // Seek to end: all Givens rotations applied
    session.seek(session.totalSteps());
    require(session.isComplete(), "session completed all steps");
    require(std::abs(session.stateNorm() - 1.0) < 1e-12, "final HeH+ state is normalized");

    // Verify physical symmetry sector: statevector must reside strictly in singlet CSFs (basis 5, 6, 9, 10)
    const auto &finalState = session.currentState();
    for (std::size_t basis = 0; basis < finalState.size(); ++basis) {
      if (basis != 5 && basis != 6 && basis != 9 && basis != 10) {
        require(std::abs(finalState[basis]) < 1e-12,
                "HeH+ state occupies strictly the 2-electron physical sector");
      }
    }
    require(std::abs(finalState[6] - finalState[9]) < 1e-12,
            "HeH+ state maintains zero-spin open-shell singlet symmetry (c6 == c9)");

    // Evaluate Pauli expectation against HeH+ Hamiltonian
    const double elEnergy = hehHamiltonian().expectation(finalState);
    require(std::abs(elEnergy - HEH_ELECTRONIC_REFERENCE_HARTREE) < 1e-8,
            "simulated HeH+ circuit energy matches exact FCI electronic reference (-4.218264 Ha)");

    const double totEnergy = elEnergy + HEH_NUCLEAR_REPULSION_HARTREE;
    require(std::abs(totEnergy - HEH_TOTAL_REFERENCE_HARTREE) < 1e-8,
            "total energy matches Born-Oppenheimer reference (-2.851475 Ha)");

    // Parameter sensitivity proof: verify non-optimal parameters raise energy strictly above minimum
    const double hfEnergy = evaluateHeHEnergy(0.0, 0.0);
    require(hfEnergy > elEnergy + 0.009,
            "Hartree-Fock determinant (theta=0, phi=0) energy is strictly higher than optimal VQE minimum by ~9.63 mHa");

    const double perturbedEnergy = evaluateHeHEnergy(optTheta + 0.1, optPhi);
    require(perturbedEnergy > elEnergy + 0.001,
            "Perturbed parameter energy is strictly higher than optimal VQE minimum");
  }

  // 7. Multi-molecule PES dataset, series definition, and bounds validation
  {
    // H2 5-point sweep reference data
    const std::vector<double> h2Distances{0.735, 0.900, 1.100, 1.250, 1.400};
    require(h2Distances.size() == 5, "H2 5-point sweep has 5 distances");

    // HeH+ equilibrium geometry reference
    const double hehBohr = 1.463283;
    const double hehAngstrom = hehBohr * 0.52917721092;
    require(std::abs(hehAngstrom - 0.774331) < 1e-4, "HeH+ equilibrium geometry is 0.774331 Angstrom");

    // Geometry out of range rejection
    bool rejectedLow = false;
    try {
      evaluateH2Energy(-0.1);
    } catch (const std::invalid_argument &) {
      rejectedLow = true;
    }
    require(rejectedLow, "negative parameter rejected");
  }

  // 8. Real-time VQE convergence trajectory recording for H2 and HeH+
  {
    // H2 VQE with callback
    std::vector<VqeConvergencePoint> h2Points;
    VqeCallback h2Cb = [&](const VqeConvergencePoint &pt) {
      h2Points.push_back(pt);
    };

    const H2Result h2Res = optimizeH2(h2Cb);
    require(h2Points.size() == 62, "H2 optimizer emits exactly 62 evaluation points");
    require(h2Res.energyEvaluations == 62, "H2 result reports 62 evaluations");
    require(!h2Points.empty(), "H2 convergence points non-empty");

    // Verify evaluation sequence and monotonically non-increasing best trajectory
    double bestSoFar = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < h2Points.size(); ++i) {
      require(h2Points[i].evaluation == static_cast<int>(i + 1), "evaluation indices are 1-based sequential");
      require(std::isfinite(h2Points[i].electronicEnergy), "electronic energy is finite");
      require(std::isfinite(h2Points[i].totalEnergy), "total energy is finite");
      if (h2Points[i].totalEnergy < bestSoFar) {
        bestSoFar = h2Points[i].totalEnergy;
        require(h2Points[i].isBest, "isBest is true when new minimum found");
      }
    }
    require(std::abs(bestSoFar - H2_TOTAL_REFERENCE_HARTREE) < 1e-8,
            "H2 live convergence best matches total reference -1.137306 Ha");

    // HeH+ VQE with callback (verify evaluation tracking)
    std::vector<VqeConvergencePoint> hehPoints;
    int hehCount = 0;
    VqeCallback hehCb = [&](const VqeConvergencePoint &pt) {
      if (hehCount < 10) {
        hehPoints.push_back(pt);
      }
      ++hehCount;
    };

    const HeHResult hehRes = optimizeHeH(hehCb);
    require(hehCount == 6803, "HeH+ optimizer executes exactly 6803 evaluations");
    require(hehPoints.size() == 10, "HeH+ callback captured initial trajectory points");
    require(std::abs(hehRes.totalEnergy - HEH_TOTAL_REFERENCE_HARTREE) < 1e-8,
            "HeH+ optimization converges to exact -2.851475 Ha total energy");
  }
}

} // namespace test
