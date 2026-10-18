#pragma once

#include "chemistry/vqe_convergence.h"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace compiler {

enum class OptimizerStatus {
  Converged,
  IterationLimitExceeded,
  EvaluationLimitExceeded,
  DegenerateSimplex,
  InvalidInput
};

inline const char *optimizerStatusToString(OptimizerStatus status) {
  switch (status) {
  case OptimizerStatus::Converged: return "Converged";
  case OptimizerStatus::IterationLimitExceeded: return "IterationLimitExceeded";
  case OptimizerStatus::EvaluationLimitExceeded: return "EvaluationLimitExceeded";
  case OptimizerStatus::DegenerateSimplex: return "DegenerateSimplex";
  case OptimizerStatus::InvalidInput: return "InvalidInput";
  default: return "Unknown";
  }
}

struct OptimizerResult {
  bool success = false;
  OptimizerStatus status = OptimizerStatus::InvalidInput;
  std::string statusMessage;
  std::vector<double> optimalParameters;
  double optimalValue = 0.0;
  int iterations = 0;
  int evaluations = 0;
};

struct OptimizerOptions {
  double tolerance = 1e-8;          // Simplex function value spread threshold
  double parameterTolerance = 1e-8; // Simplex coordinate diameter threshold
  int maxIterations = 500;         // Maximum simplex transformations
  int maxEvaluations = 2000;       // Hard evaluation ceiling
  double initialStep = 0.2;        // Initial simplex edge step size
  double alpha = 1.0;              // Reflection coefficient (> 0)
  double gamma = 2.0;              // Expansion coefficient (> alpha)
  double rho = 0.5;                // Contraction coefficient (0 < rho < 1)
  double sigma = 0.5;              // Shrink coefficient (0 < sigma < 1)
};

using MultiVqeCallback = std::function<void(
    int evaluation,
    const std::vector<double> &parameters,
    double value,
    bool isBest
)>;

class NelderMeadOptimizer {
public:
  static OptimizerResult minimize(
      std::function<double(const std::vector<double> &)> objective,
      const std::vector<double> &initialGuess,
      const std::vector<double> &lowerBounds = {},
      const std::vector<double> &upperBounds = {},
      const OptimizerOptions &options = {},
      MultiVqeCallback callback = nullptr
  );
};

// Adapter to forward MultiVqeCallback events to existing legacy VqeCallback
MultiVqeCallback makeVqeCallbackAdapter(VqeCallback callback, double nuclearRepulsion = 0.0);

} // namespace compiler
