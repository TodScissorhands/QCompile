#pragma once

#include "compiler/code_generator.h"
#include "simulator/simulator.h"

#include <complex>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

class SimulationSession {
public:
  // Constructs a simulation session for the given circuit and qubit count.
  // Validates qubitCount (0 <= qubitCount <= MAX_SIMULATOR_QUBITS) and initial allocation safety.
  // Initializes state to |0...0> at step 0.
  SimulationSession(std::vector<Instruction> instructions, int qubitCount);

  // Step definition:
  // - currentStep() counts the number of instructions executed so far (range: [0, totalSteps()]).
  // - When currentStep() == 0, 0 instructions have been executed (initial |0...0> state).
  // - When currentStep() == k, the first k instructions (indices 0 to k-1) have been executed.
  // - When currentStep() < totalSteps(), currentStep() is also the 0-based index of the next instruction.
  // - isComplete() is true when currentStep() == totalSteps().

  // Advances simulation by executing the next instruction.
  // If the next instruction has been executed previously in this session, the recorded
  // state and measurement outcome are restored without re-sampling.
  // Returns true if an instruction was executed; false if already complete.
  bool stepForward();

  // Rewinds simulation by one instruction, restoring the exact state prior to that instruction.
  // Returns true if rewound; false if already at step 0.
  bool stepBackward();

  // Resets execution to step 0.
  // By default (clearHistory = false), recorded snapshots and measurement outcomes are
  // retained so subsequent stepForward() calls reproduce the exact historical outcomes.
  // If clearHistory is true, all snapshots after step 0 are discarded.
  void reset(bool clearHistory = false);

  // Seeks to an arbitrary step in [0, totalSteps()].
  // Throws std::out_of_range if targetStep > totalSteps().
  void seek(std::size_t targetStep);

  // Current execution progress: number of instructions already executed.
  std::size_t currentStep() const;

  // Total number of instructions in the circuit.
  std::size_t totalSteps() const;

  // Read-only reference to the current statevector (dimension: 2^qubitCount).
  const std::vector<std::complex<double>> &currentState() const;

  // Basis state probabilities |c_i|^2 for all 2^qubitCount computational basis states.
  std::vector<double> basisProbabilities() const;

  // Euclidean norm (L2) of the current statevector.
  double stateNorm() const;

  // True if all instructions have been executed (currentStep() == totalSteps()).
  bool isComplete() const;

  // Qubit count and instruction access.
  int qubitCount() const;
  const std::vector<Instruction> &instructions() const;

  // Measurement outcome queries.
  // measurementResultAt(step) returns the recorded outcome bitstring if instruction at (step - 1)
  // was a measurement; returns std::nullopt otherwise.
  std::optional<std::string> measurementResultAt(std::size_t step) const;
  std::optional<std::string> lastMeasurementResult() const;

private:
  struct Snapshot {
    std::vector<std::complex<double>> state;
    std::optional<std::size_t> measurementOutcome;
    std::string measurementResultString;
  };

  int qubitCount_;
  std::vector<Instruction> instructions_;
  std::size_t currentStep_;
  Simulator simulator_;
  std::vector<Snapshot> snapshots_;

  void executeStep(std::size_t stepIndex);
};
