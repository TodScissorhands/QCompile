#include "simulator/simulation_session.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

// Configured statevector payload memory limit across the active simulator buffer
// and all retained step snapshots.
//
// Memory accounting definition and boundaries:
// The 512 MiB bound applies strictly to statevector payload allocations:
//   total_statevector_bytes = total_statevector_buffers * stateBytes
// where stateBytes = 2^qubitCount * sizeof(std::complex<double>).
// At any execution point, total_statevector_buffers = snapshots_.size() + 1
// (the retained historical snapshots plus the active simulator_.state buffer).
//
// This limit intentionally excludes fixed process metadata, heap allocator alignment
// padding, and small object control blocks (e.g. std::string result text), as those
// platform-dependent heap costs cannot be portably measured. It guarantees a deterministic
// and predictable resource bound for primary statevector payload memory.
constexpr std::size_t kMaxSessionStateMemoryBytes = 512 * 1024 * 1024; // 512 MiB limit

std::size_t statevectorDimension(int qubitCount) {
  if (qubitCount < 0 || qubitCount > MAX_SIMULATOR_QUBITS) {
    throw std::invalid_argument(
        "Simulator supports between 0 and " +
        std::to_string(MAX_SIMULATOR_QUBITS) + " qubits"
    );
  }
  return std::size_t{1} << static_cast<unsigned int>(qubitCount);
}

// Overflow-safe check whether allocating requiredBuffers of stateBytes exceeds maxLimitBytes.
bool wouldExceedStateMemoryLimit(
    std::size_t requiredBuffers,
    std::size_t stateBytes,
    std::size_t maxLimitBytes
) {
  if (stateBytes == 0 || requiredBuffers == 0) {
    return false;
  }
  const std::size_t maxAllowedBuffers = maxLimitBytes / stateBytes;
  return requiredBuffers > maxAllowedBuffers;
}

} // namespace

SimulationSession::SimulationSession(
    std::vector<Instruction> instructions,
    int qubitCount
)
    : qubitCount_(qubitCount),
      instructions_(std::move(instructions)),
      currentStep_(0),
      simulator_(qubitCount) {
  const std::size_t dimension = statevectorDimension(qubitCount_);
  const std::size_t stateBytes = dimension * sizeof(std::complex<double>);
  // Initial session holds 2 buffers: the active simulator state plus the step 0 snapshot.
  constexpr std::size_t kInitialBuffers = 2;
  if (wouldExceedStateMemoryLimit(kInitialBuffers, stateBytes, kMaxSessionStateMemoryBytes)) {
    throw std::length_error(
        "Initial session statevector memory exceeds configured limit"
    );
  }

  // Pre-validate instruction qubit indices to fail early on malformed circuits
  for (std::size_t i = 0; i < instructions_.size(); ++i) {
    const auto &instr = instructions_[i];
    for (int qubit : instr.qubitIndices) {
      if (qubit < 0 || qubit >= qubitCount_) {
        throw std::invalid_argument(
            "Instruction at index " + std::to_string(i) +
            " references out-of-range qubit " + std::to_string(qubit)
        );
      }
    }
  }

  // Snapshot at step 0: initial all-zero computational basis state
  snapshots_.push_back(Snapshot{
      simulator_.stateVector(),
      std::nullopt,
      ""
  });
}

bool SimulationSession::stepForward() {
  if (isComplete()) {
    return false;
  }

  const std::size_t nextStep = currentStep_ + 1;
  if (nextStep < snapshots_.size()) {
    // Restore already recorded snapshot from history without re-sampling
    currentStep_ = nextStep;
    simulator_.setStateVector(snapshots_[currentStep_].state);
    return true;
  }

  // Check memory limit before computing and storing next state.
  // The session will retain snapshots_.size() existing snapshots, plus 1 new snapshot,
  // plus 1 active simulator statevector buffer: total = snapshots_.size() + 2.
  const std::size_t dimension = statevectorDimension(qubitCount_);
  const std::size_t stateBytes = dimension * sizeof(std::complex<double>);
  if (snapshots_.size() > std::numeric_limits<std::size_t>::max() - 2) {
    throw std::length_error("Snapshot count exceeds maximum index limit");
  }
  const std::size_t requiredBuffers = snapshots_.size() + 2;
  if (wouldExceedStateMemoryLimit(requiredBuffers, stateBytes, kMaxSessionStateMemoryBytes)) {
    throw std::length_error(
        "SimulationSession statevector memory exceeds maximum limit"
    );
  }

  const Instruction &instruction = instructions_[currentStep_];
  StepResult stepResult = simulator_.applyInstruction(
      instruction,
      std::nullopt,
      /*printOutput=*/false
  );

  snapshots_.push_back(Snapshot{
      simulator_.stateVector(),
      stepResult.measurementOutcome,
      stepResult.measurementResult
  });

  currentStep_ = nextStep;
  return true;
}

bool SimulationSession::stepBackward() {
  if (currentStep_ == 0) {
    return false;
  }

  currentStep_ -= 1;
  simulator_.setStateVector(snapshots_[currentStep_].state);
  return true;
}

void SimulationSession::reset(bool clearHistory) {
  currentStep_ = 0;
  simulator_.setStateVector(snapshots_[0].state);

  if (clearHistory && snapshots_.size() > 1) {
    snapshots_.resize(1);
  }
}

void SimulationSession::seek(std::size_t targetStep) {
  if (targetStep > totalSteps()) {
    throw std::out_of_range(
        "Seek step " + std::to_string(targetStep) +
        " exceeds total instructions " + std::to_string(totalSteps())
    );
  }

  while (currentStep_ < targetStep) {
    stepForward();
  }
  while (currentStep_ > targetStep) {
    stepBackward();
  }
}

std::size_t SimulationSession::currentStep() const {
  return currentStep_;
}

std::size_t SimulationSession::totalSteps() const {
  return instructions_.size();
}

const std::vector<std::complex<double>> &SimulationSession::currentState() const {
  return snapshots_[currentStep_].state;
}

std::vector<double> SimulationSession::basisProbabilities() const {
  const auto &state = currentState();
  std::vector<double> probabilities(state.size(), 0.0);
  for (std::size_t i = 0; i < state.size(); ++i) {
    probabilities[i] = std::norm(state[i]);
  }
  return probabilities;
}

double SimulationSession::stateNorm() const {
  const auto &state = currentState();
  double sum = 0.0;
  for (const auto &amplitude : state) {
    sum += std::norm(amplitude);
  }
  return std::sqrt(sum);
}

bool SimulationSession::isComplete() const {
  return currentStep_ == instructions_.size();
}

int SimulationSession::qubitCount() const {
  return qubitCount_;
}

const std::vector<Instruction> &SimulationSession::instructions() const {
  return instructions_;
}

std::optional<std::string> SimulationSession::measurementResultAt(
    std::size_t step
) const {
  if (step == 0 || step > currentStep_ || step >= snapshots_.size()) {
    return std::nullopt;
  }
  if (snapshots_[step].measurementOutcome.has_value()) {
    return snapshots_[step].measurementResultString;
  }
  return std::nullopt;
}

std::optional<std::string> SimulationSession::lastMeasurementResult() const {
  if (currentStep_ == 0) {
    return std::nullopt;
  }
  return measurementResultAt(currentStep_);
}
