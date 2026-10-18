#include "chemistry/nelder_mead.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace compiler {

namespace {

struct ObjectiveEvaluation {
  double value = 0.0;
  bool finite = false;
};

struct SimplexVertex {
  std::vector<double> point;
  double value = 0.0;
  bool finite = false;
};

bool isBetter(bool aFinite, double aValue, bool bFinite, double bValue) {
  if (aFinite != bFinite) return aFinite;
  return aValue < bValue;
}

} // namespace

OptimizerResult NelderMeadOptimizer::minimize(
    std::function<double(const std::vector<double> &)> objective,
    const std::vector<double> &initialGuess,
    const std::vector<double> &lowerBounds,
    const std::vector<double> &upperBounds,
    const OptimizerOptions &options,
    MultiVqeCallback callback
) {
  OptimizerResult result;

  // 1. Validation of input parameters
  if (!objective) {
    result.status = OptimizerStatus::InvalidInput;
    result.statusMessage = "Objective function must not be empty";
    return result;
  }

  if (initialGuess.empty()) {
    result.status = OptimizerStatus::InvalidInput;
    result.statusMessage = "Initial guess must have positive dimension";
    return result;
  }

  for (double coordinate : initialGuess) {
    if (!std::isfinite(coordinate)) {
      result.status = OptimizerStatus::InvalidInput;
      result.statusMessage = "Initial guess coordinates must be finite";
      return result;
    }
  }

  const std::size_t n = initialGuess.size();
  const bool hasLower = !lowerBounds.empty();
  const bool hasUpper = !upperBounds.empty();

  if ((hasLower && lowerBounds.size() != n) || (hasUpper && upperBounds.size() != n)) {
    result.status = OptimizerStatus::InvalidInput;
    result.statusMessage = "Bounds dimensions must match initial guess dimension";
    return result;
  }

  if (hasLower && hasUpper) {
    for (std::size_t i = 0; i < n; ++i) {
      if (lowerBounds[i] > upperBounds[i]) {
        result.status = OptimizerStatus::InvalidInput;
        result.statusMessage = "Lower bound exceeds upper bound at coordinate " + std::to_string(i);
        return result;
      }
    }
  }

  if (options.tolerance <= 0.0 || options.parameterTolerance <= 0.0 ||
      options.maxIterations <= 0 || options.maxEvaluations <= 0 ||
      options.initialStep <= 0.0 || options.alpha <= 0.0 ||
      options.gamma <= options.alpha || options.rho <= 0.0 ||
      options.rho >= 1.0 || options.sigma <= 0.0 || options.sigma >= 1.0) {
    result.status = OptimizerStatus::InvalidInput;
    result.statusMessage = "Invalid optimizer options configuration";
    return result;
  }

  // 2. Bound projection / penalty helper
  const bool hasBounds = hasLower && hasUpper;
  auto clampPoint = [&](const std::vector<double> &p) {
    if (!hasBounds) return p;
    std::vector<double> clamped = p;
    for (std::size_t d = 0; d < n; ++d) {
      if (clamped[d] < lowerBounds[d]) clamped[d] = lowerBounds[d];
      if (clamped[d] > upperBounds[d]) clamped[d] = upperBounds[d];
    }
    return clamped;
  };

  int evalCount = 0;
  bool observedFiniteValue = false;
  bool hasBestValueSeen = false;
  bool bestValueFinite = false;
  double bestValueSeen = std::numeric_limits<double>::infinity();

  auto evalFunction = [&](const std::vector<double> &rawPoint) -> ObjectiveEvaluation {
    if (evalCount >= options.maxEvaluations) {
      return {std::numeric_limits<double>::infinity(), false};
    }

    std::vector<double> clampedPoint = clampPoint(rawPoint);
    double val = objective(clampedPoint);
    ++evalCount;

    const bool objectiveFinite = std::isfinite(val);
    if (objectiveFinite) {
      observedFiniteValue = true;
    } else {
      val = 1e12; // callback penalty for a non-finite objective result
    }

    // Boundary penalty if rawPoint is outside bounds
    if (hasBounds) {
      double penalty = 0.0;
      for (std::size_t d = 0; d < n; ++d) {
        if (rawPoint[d] < lowerBounds[d]) {
          const double diff = lowerBounds[d] - rawPoint[d];
          penalty += 1e4 * diff * diff + 100.0 * diff;
        } else if (rawPoint[d] > upperBounds[d]) {
          const double diff = rawPoint[d] - upperBounds[d];
          penalty += 1e4 * diff * diff + 100.0 * diff;
        }
      }
      val += penalty;
    }

    const bool isBest = !hasBestValueSeen ||
        isBetter(objectiveFinite, val, bestValueFinite, bestValueSeen);
    if (isBest) {
      hasBestValueSeen = true;
      bestValueFinite = objectiveFinite;
      bestValueSeen = val;
    }
    if (callback) {
      callback(evalCount, clampedPoint, val, isBest);
    }
    return {val, objectiveFinite};
  };

  // 3. Construct initial simplex of n + 1 vertices
  std::vector<SimplexVertex> simplex(n + 1);
  std::vector<double> x0 = clampPoint(initialGuess);
  simplex[0].point = x0;
  const ObjectiveEvaluation initialEvaluation = evalFunction(x0);
  simplex[0].value = initialEvaluation.value;
  simplex[0].finite = initialEvaluation.finite;
  for (std::size_t i = 0; i < n; ++i) {
    std::vector<double> xi = x0;
    if (hasBounds) {
      if (x0[i] + options.initialStep <= upperBounds[i]) {
        xi[i] = x0[i] + options.initialStep;
      } else if (x0[i] - options.initialStep >= lowerBounds[i]) {
        xi[i] = x0[i] - options.initialStep;
      } else {
        xi[i] = 0.5 * (lowerBounds[i] + upperBounds[i]);
      }
      // Guarantee non-degenerate simplex
      if (std::abs(xi[i] - x0[i]) < 1e-14) {
        xi[i] = lowerBounds[i] + 0.05 * (upperBounds[i] - lowerBounds[i]);
      }
    } else {
      xi[i] = x0[i] + options.initialStep;
    }
    simplex[i + 1].point = xi;
    const ObjectiveEvaluation evaluation = evalFunction(xi);
    simplex[i + 1].value = evaluation.value;
    simplex[i + 1].finite = evaluation.finite;
  }

  // 4. Downhill Simplex Transformations Loop
  int iter = 0;
  OptimizerStatus finalStatus = OptimizerStatus::IterationLimitExceeded;

  for (iter = 0; iter < options.maxIterations; ++iter) {
    // Rank finite objective evaluations first, then by value.
    std::sort(simplex.begin(), simplex.end(), [](const SimplexVertex &a, const SimplexVertex &b) {
      return isBetter(a.finite, a.value, b.finite, b.value);
    });

    // Convergence check 1: Function value standard deviation across simplex
    double meanVal = 0.0;
    for (const auto &v : simplex) meanVal += v.value;
    meanVal /= static_cast<double>(n + 1);

    double stddevSq = 0.0;
    for (const auto &v : simplex) {
      const double diff = v.value - meanVal;
      stddevSq += diff * diff;
    }
    const double fSpread = std::sqrt(stddevSq / static_cast<double>(n + 1));

    // Convergence check 2: Simplex coordinate diameter (infinity norm)
    double coordDiameter = 0.0;
    for (std::size_t i = 1; i <= n; ++i) {
      for (std::size_t d = 0; d < n; ++d) {
        coordDiameter = std::max(coordDiameter, std::abs(simplex[i].point[d] - simplex[0].point[d]));
      }
    }

    if (fSpread < options.tolerance && coordDiameter < options.parameterTolerance) {
      finalStatus = OptimizerStatus::Converged;
      break;
    }

    if (evalCount >= options.maxEvaluations) {
      finalStatus = OptimizerStatus::EvaluationLimitExceeded;
      break;
    }

    // Compute centroid of the best n vertices (excluding worst at index n)
    std::vector<double> xCentroid(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t d = 0; d < n; ++d) {
        xCentroid[d] += simplex[i].point[d];
      }
    }
    for (std::size_t d = 0; d < n; ++d) {
      xCentroid[d] /= static_cast<double>(n);
    }

    // Transformation 1: Reflection
    // xr = xCentroid + alpha * (xCentroid - worst)
    std::vector<double> xr(n);
    for (std::size_t d = 0; d < n; ++d) {
      xr[d] = xCentroid[d] + options.alpha * (xCentroid[d] - simplex[n].point[d]);
    }
    const ObjectiveEvaluation fxr = evalFunction(xr);

    if (!isBetter(fxr.finite, fxr.value, simplex[0].finite, simplex[0].value) &&
        isBetter(fxr.finite, fxr.value, simplex[n - 1].finite, simplex[n - 1].value)) {
      simplex[n].point = xr;
      simplex[n].value = fxr.value;
      simplex[n].finite = fxr.finite;
      continue;
    }

    // Transformation 2: Expansion
    if (isBetter(fxr.finite, fxr.value, simplex[0].finite, simplex[0].value)) {
      // xe = xCentroid + gamma * (xr - xCentroid)
      std::vector<double> xe(n);
      for (std::size_t d = 0; d < n; ++d) {
        xe[d] = xCentroid[d] + options.gamma * (xr[d] - xCentroid[d]);
      }
      const ObjectiveEvaluation fxe = evalFunction(xe);
      const bool expansionBetter = isBetter(fxe.finite, fxe.value, fxr.finite, fxr.value);
      const ObjectiveEvaluation &selected = expansionBetter ? fxe : fxr;
      simplex[n].point = expansionBetter ? xe : xr;
      simplex[n].value = selected.value;
      simplex[n].finite = selected.finite;
      continue;
    }

    // Transformation 3: Contraction
    if (isBetter(fxr.finite, fxr.value, simplex[n].finite, simplex[n].value)) {
      // Outside contraction: xc = xCentroid + rho * (xr - xCentroid)
      std::vector<double> xc(n);
      for (std::size_t d = 0; d < n; ++d) {
        xc[d] = xCentroid[d] + options.rho * (xr[d] - xCentroid[d]);
      }
      const ObjectiveEvaluation fxc = evalFunction(xc);
      if (!isBetter(fxr.finite, fxr.value, fxc.finite, fxc.value)) {
        simplex[n].point = xc;
        simplex[n].value = fxc.value;
        simplex[n].finite = fxc.finite;
        continue;
      }
    } else {
      // Inside contraction: xc = xCentroid + rho * (worst - xCentroid)
      std::vector<double> xc(n);
      for (std::size_t d = 0; d < n; ++d) {
        xc[d] = xCentroid[d] + options.rho * (simplex[n].point[d] - xCentroid[d]);
      }
      const ObjectiveEvaluation fxc = evalFunction(xc);
      if (isBetter(fxc.finite, fxc.value, simplex[n].finite, simplex[n].value)) {
        simplex[n].point = xc;
        simplex[n].value = fxc.value;
        simplex[n].finite = fxc.finite;
        continue;
      }
    }

    // Transformation 4: Shrink towards best vertex v0
    for (std::size_t i = 1; i <= n; ++i) {
      for (std::size_t d = 0; d < n; ++d) {
        simplex[i].point[d] = simplex[0].point[d] + options.sigma * (simplex[i].point[d] - simplex[0].point[d]);
      }
      const ObjectiveEvaluation evaluation = evalFunction(simplex[i].point);
      simplex[i].value = evaluation.value;
      simplex[i].finite = evaluation.finite;
    }
  }

  // Final sorting keeps every finite objective ahead of non-finite evaluations.
  std::sort(simplex.begin(), simplex.end(), [](const SimplexVertex &a, const SimplexVertex &b) {
    return isBetter(a.finite, a.value, b.finite, b.value);
  });

  result.optimalParameters = hasBounds ? clampPoint(simplex[0].point) : simplex[0].point;
  result.optimalValue = simplex[0].value;
  result.iterations = iter;
  result.evaluations = evalCount;
  if (!observedFiniteValue || !simplex[0].finite) {
    result.status = OptimizerStatus::DegenerateSimplex;
    result.success = false;
    result.statusMessage = observedFiniteValue
        ? "No finite objective solution remains in the final simplex"
        : "All objective function evaluations were non-finite";
  } else {
    result.status = finalStatus;
    result.success = (finalStatus == OptimizerStatus::Converged);
    result.statusMessage = optimizerStatusToString(finalStatus);
  }

  return result;
}

MultiVqeCallback makeVqeCallbackAdapter(VqeCallback callback, double nuclearRepulsion) {
  if (!callback) return nullptr;
  return [callback, nuclearRepulsion](
      int evaluation,
      const std::vector<double> &parameters,
      double value,
      bool isBest
  ) {
    VqeConvergencePoint pt;
    pt.evaluation = evaluation;
    pt.parameter0 = !parameters.empty() ? parameters[0] : 0.0;
    pt.parameter1 = (parameters.size() > 1) ? parameters[1] : 0.0;
    pt.electronicEnergy = value;
    pt.totalEnergy = value + nuclearRepulsion;
    pt.isBest = isBest;
    callback(pt);
  };
}

} // namespace compiler
