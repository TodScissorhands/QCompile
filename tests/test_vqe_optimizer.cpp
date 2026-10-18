#include "test_support.h"

#include "chemistry/h2_energy.h"
#include "chemistry/heh_energy.h"
#include "chemistry/molecular_integral_engine.h"
#include "chemistry/nelder_mead.h"
#include "chemistry/vqe_ansatz.h"
#include "chemistry/vqe_convergence.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace test {

void testVqeOptimizer() {
  // =========================================================================
  // 1. Analytical Test Functions: 2D & 4D Sphere, 2D Rosenbrock
  // =========================================================================
  {
    // 1a. 2D Sphere: f(x, y) = x^2 + y^2, minimum at (0, 0), f = 0
    auto resSphere2D = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) {
          return p[0] * p[0] + p[1] * p[1];
        },
        {1.5, -2.0}
    );
    require(resSphere2D.success, "2D Sphere optimization converges");
    require(resSphere2D.status == compiler::OptimizerStatus::Converged, "Status is Converged");
    require(std::abs(resSphere2D.optimalParameters[0]) < 1e-4, "Sphere 2D optimal x is ~0");
    require(std::abs(resSphere2D.optimalParameters[1]) < 1e-4, "Sphere 2D optimal y is ~0");
    require(resSphere2D.optimalValue < 1e-8, "Sphere 2D optimal value is < 1e-8");

    // 1b. 4D Sphere: f(x) = sum x_i^2, minimum at origin
    auto resSphere4D = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) {
          double s = 0.0;
          for (double x : p) s += x * x;
          return s;
        },
        {1.0, -1.0, 0.5, -0.5}
    );
    require(resSphere4D.success, "4D Sphere optimization converges");
    require(resSphere4D.optimalValue < 1e-7, "Sphere 4D optimal value is < 1e-7");

    // 1c. 2D Rosenbrock: f(x, y) = 100*(y - x^2)^2 + (1 - x)^2, minimum at (1, 1), f = 0
    compiler::OptimizerOptions rosenbrockOpt;
    rosenbrockOpt.maxIterations = 1000;
    rosenbrockOpt.maxEvaluations = 3000;
    rosenbrockOpt.tolerance = 1e-9;
    auto resRosenbrock = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) {
          const double x = p[0];
          const double y = p[1];
          const double d1 = y - x * x;
          const double d2 = 1.0 - x;
          return 100.0 * d1 * d1 + d2 * d2;
        },
        {-1.2, 1.0},
        {},
        {},
        rosenbrockOpt
    );
    require(resRosenbrock.success, "2D Rosenbrock optimization converges");
    require(std::abs(resRosenbrock.optimalParameters[0] - 1.0) < 1e-3, "Rosenbrock x is ~1.0");
    require(std::abs(resRosenbrock.optimalParameters[1] - 1.0) < 1e-3, "Rosenbrock y is ~1.0");
    require(resRosenbrock.optimalValue < 1e-6, "Rosenbrock optimal value is < 1e-6");
  }

  // =========================================================================
  // 2. Input Validation, Bounds Enforcement, Limits & Edge Cases
  // =========================================================================
  {
    // Empty objective
    auto resNull = compiler::NelderMeadOptimizer::minimize(nullptr, {1.0});
    require(!resNull.success, "Null objective fails");
    require(resNull.status == compiler::OptimizerStatus::InvalidInput, "Status is InvalidInput");

    // Empty initial guess
    auto resEmpty = compiler::NelderMeadOptimizer::minimize([](const std::vector<double> &) { return 0.0; }, {});
    require(!resEmpty.success, "Empty initial guess fails");
    require(resEmpty.status == compiler::OptimizerStatus::InvalidInput, "Status is InvalidInput");


    // Non-finite initial guess coordinates are rejected before objective evaluation.
    const std::vector<double> invalidCoordinates{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()
    };
    for (double coordinate : invalidCoordinates) {
      int objectiveEvaluations = 0;
      const auto resInvalidGuess = compiler::NelderMeadOptimizer::minimize(
          [&](const std::vector<double> &) {
            ++objectiveEvaluations;
            return 0.0;
          },
          {coordinate}
      );
      require(!resInvalidGuess.success, "Non-finite initial guess fails");
      require(resInvalidGuess.status == compiler::OptimizerStatus::InvalidInput,
              "Non-finite initial guess reports InvalidInput");
      require(resInvalidGuess.evaluations == 0 && objectiveEvaluations == 0,
              "Non-finite initial guess is rejected before objective evaluation");
      require(resInvalidGuess.status != compiler::OptimizerStatus::Converged &&
                  resInvalidGuess.optimalParameters.empty(),
              "Non-finite initial guess returns no converged solution");
      require(resInvalidGuess.statusMessage == "Initial guess coordinates must be finite",
              "Non-finite initial guess reports clear status message");
    }
    // Mismatched lower bounds dimension
    auto resBoundMismatch = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) { return p[0]; },
        {1.0, 2.0},
        {0.0}
    );
    require(!resBoundMismatch.success, "Mismatched lower bounds fails");
    require(resBoundMismatch.status == compiler::OptimizerStatus::InvalidInput, "Status is InvalidInput");

    // Inverted bounds (lower > upper)
    auto resInvertedBounds = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) { return p[0]; },
        {1.0},
        {2.0},
        {1.0}
    );
    require(!resInvertedBounds.success, "Inverted bounds fail");
    require(resInvertedBounds.status == compiler::OptimizerStatus::InvalidInput, "Status is InvalidInput");

    // Invalid options (e.g. negative tolerance, rho > 1)
    compiler::OptimizerOptions badOptions;
    badOptions.tolerance = -1.0;
    auto resBadOpt = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) { return p[0]; },
        {1.0},
        {},
        {},
        badOptions
    );
    require(!resBadOpt.success, "Negative tolerance fails");
    require(resBadOpt.status == compiler::OptimizerStatus::InvalidInput, "Status is InvalidInput");

    // Box bounds enforcement: ensure optimizer strictly stays within bounds
    // Objective: f(x) = (x - 5)^2 on [-2, 2]. Unconstrained minimum is at 5, but bounded minimum is at 2.0!
    auto resBounded = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) {
          const double x = p[0];
          return (x - 5.0) * (x - 5.0);
        },
        {0.0},
        {-2.0},
        {2.0}
    );
    require(resBounded.success, "Bounded quadratic optimization converges");
    require(std::abs(resBounded.optimalParameters[0] - 2.0) < 1e-6, "Bounded minimum is clamped at boundary 2.0");
    require(std::abs(resBounded.optimalValue - 9.0) < 1e-6, "Optimal value is (2 - 5)^2 = 9.0");

    // Non-finite objective recovery: function returns NaN outside safe region
    auto resNanRecovery = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) {
          if (p[0] < 0.0) return std::numeric_limits<double>::quiet_NaN();
          return (p[0] - 1.0) * (p[0] - 1.0);
        },
        {0.5}
    );
    require(resNanRecovery.success, "Optimizer recovers from non-finite objective region");
    require(std::abs(resNanRecovery.optimalParameters[0] - 1.0) < 1e-4, "Optimal x is 1.0");

    // Non-finite evaluations cannot outrank a finite objective above the
    // historical 1e12 penalty.
    auto resFiniteAbovePenalty = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) {
          if (p[0] != 0.0) return std::numeric_limits<double>::quiet_NaN();
          return 1e13;
        },
        {0.0}
    );
    require(!resFiniteAbovePenalty.success ||
                (std::isfinite(resFiniteAbovePenalty.optimalValue) &&
                 resFiniteAbovePenalty.optimalValue == 1e13 &&
                 resFiniteAbovePenalty.optimalParameters.size() == 1 &&
                 resFiniteAbovePenalty.optimalParameters[0] == 0.0),
            "A successful result from the exact regression has a finite objective");
    require(resFiniteAbovePenalty.optimalParameters.size() == 1 &&
                resFiniteAbovePenalty.optimalParameters[0] == 0.0,
            "Finite 1e13 point remains ahead of non-finite objective points");
    require(resFiniteAbovePenalty.optimalValue == 1e13,
            "Finite 1e13 objective value is preserved");

    auto resFiniteAtPenalty = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &) { return 1e12; },
        {0.0}
    );
    require(resFiniteAtPenalty.success &&
                resFiniteAtPenalty.optimalValue == 1e12,
            "Finite objective equal to 1e12 remains valid");

    // Mixed valid/non-finite evaluations still converge when the finite
    // objective is above 1e12.
    auto resMixedLargeFinite = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) {
          if (p[0] < 0.2) return std::numeric_limits<double>::quiet_NaN();
          return 1e13 + (p[0] - 0.5) * (p[0] - 0.5);
        },
        {0.0}
    );
    require(resMixedLargeFinite.success,
            "Mixed non-finite and large finite objective converges");
    require(std::isfinite(resMixedLargeFinite.optimalValue) &&
                resMixedLargeFinite.optimalValue >= 1e13,
            "Mixed objective returns a finite value above 1e12");
    require(resMixedLargeFinite.optimalParameters[0] >= 0.2,
            "Mixed objective result is in its finite domain");
    require(std::abs(resMixedLargeFinite.optimalParameters[0] - 0.5) < 1e-4,
            "Mixed objective converges to its finite-domain minimum");

    // All NaN evaluations failure
    auto resAllNaN = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &) { return std::numeric_limits<double>::quiet_NaN(); },
        {1.0}
    );
    require(!resAllNaN.success, "All-NaN objective fails");
    require(resAllNaN.status == compiler::OptimizerStatus::DegenerateSimplex, "All-NaN reports DegenerateSimplex");
    require(resAllNaN.status != compiler::OptimizerStatus::Converged, "All-NaN does not report Converged");
    require(resAllNaN.statusMessage == "All objective function evaluations were non-finite",
            "All-NaN reports descriptive failure message");

    // All +Infinity evaluations failure
    auto resAllPosInf = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &) { return std::numeric_limits<double>::infinity(); },
        {1.0}
    );
    require(!resAllPosInf.success, "All +Inf objective fails");
    require(resAllPosInf.status == compiler::OptimizerStatus::DegenerateSimplex, "All +Inf reports DegenerateSimplex");
    require(resAllPosInf.status != compiler::OptimizerStatus::Converged, "All +Inf does not report Converged");
    require(resAllPosInf.statusMessage == "All objective function evaluations were non-finite",
            "All +Inf reports descriptive failure message");

    // All -Infinity evaluations failure
    auto resAllNegInf = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &) { return -std::numeric_limits<double>::infinity(); },
        {1.0}
    );
    require(!resAllNegInf.success, "All -Inf objective fails");
    require(resAllNegInf.status == compiler::OptimizerStatus::DegenerateSimplex, "All -Inf reports DegenerateSimplex");
    require(resAllNegInf.status != compiler::OptimizerStatus::Converged, "All -Inf does not report Converged");
    require(resAllNegInf.statusMessage == "All objective function evaluations were non-finite",
            "All -Inf reports descriptive failure message");

    // Evaluation limit exceeded
    compiler::OptimizerOptions evalLimitOptions;
    evalLimitOptions.maxEvaluations = 10;
    auto resEvalLimit = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) { return p[0] * p[0]; },
        {1.0},
        {},
        {},
        evalLimitOptions
    );
    require(!resEvalLimit.success, "Evaluation limit exceeded correctly fails");
    require(resEvalLimit.status == compiler::OptimizerStatus::EvaluationLimitExceeded,
            "Status is EvaluationLimitExceeded");

    // Iteration limit exceeded
    compiler::OptimizerOptions iterLimitOptions;
    iterLimitOptions.maxIterations = 3;
    auto resIterLimit = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) { return p[0] * p[0]; },
        {1.0},
        {},
        {},
        iterLimitOptions
    );
    require(!resIterLimit.success, "Iteration limit exceeded correctly fails");
    require(resIterLimit.status == compiler::OptimizerStatus::IterationLimitExceeded,
            "Status is IterationLimitExceeded");
  }

  // =========================================================================
  // 3. MultiVqeCallback & Legacy Adapter Verification
  // =========================================================================
  {
    int callbackInvocations = 0;
    int bestCount = 0;
    auto resCallback = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) { return p[0] * p[0]; },
        {2.0},
        {},
        {},
        {},
        [&](int eval, const std::vector<double> &, double, bool isBest) {
          ++callbackInvocations;
          require(eval == callbackInvocations, "Callback evaluation index increases monotonically");
          if (isBest) ++bestCount;
        }
    );
    require(callbackInvocations == resCallback.evaluations,
            "Callback invocation count matches total evaluations");
    require(bestCount > 0, "At least one 'isBest' point was recorded");

    // Legacy VqeCallback adapter check
    int legacyInvocations = 0;
    VqeCallback legacyCallback = [&](const VqeConvergencePoint &pt) {
      ++legacyInvocations;
      require(pt.totalEnergy == pt.electronicEnergy + 1.5, "Adapter applies nuclear repulsion shift to totalEnergy");
    };
    auto adapter = compiler::makeVqeCallbackAdapter(legacyCallback, 1.5);
    (void)compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) { return p[0] * p[0]; },
        {1.0},
        {},
        {},
        {},
        adapter
    );
    require(legacyInvocations > 0, "Legacy VqeCallback adapter forwards evaluations successfully");
  }

  // =========================================================================
  // 4. Molecular VQE Optimization: H2 (1D) and HeH+ (2D)
  // =========================================================================
  {
    // 4a. H2 1-Parameter Ground State Optimization
    // Reference: -1.857275029288 Ha
    const double kPi = 3.141592653589793238462643383279502884;
    auto h2Res = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) {
          return evaluateH2Energy(p[0]);
        },
        {0.1},
        {0.0},
        {kPi / 2.0}
    );
    require(h2Res.success, "Nelder-Mead H2 optimization converges");
    const double h2Min = evaluateH2Energy(h2Res.optimalParameters[0]);
    require(std::abs(h2Min - (-1.857275029288)) < 1e-8,
            "Nelder-Mead H2 ground energy matches reference within 1e-8 Ha");
    require(h2Res.evaluations < 62, "Nelder-Mead H2 converges in fewer evaluations than golden section (< 62)");

    // 4b. HeH+ 2-Parameter Ground State Optimization
    // Reference: -4.218264154490 Ha
    auto hehRes = compiler::NelderMeadOptimizer::minimize(
        [](const std::vector<double> &p) {
          return evaluateHeHEnergy(p[0], p[1]);
        },
        {0.1, 4.9},
        {0.0, 0.0},
        {kPi, 2.0 * kPi - 1e-6}
    );
    require(hehRes.success, "Nelder-Mead HeH+ optimization converges");
    const double hehMin = evaluateHeHEnergy(hehRes.optimalParameters[0], hehRes.optimalParameters[1]);
    require(std::abs(hehMin - (-4.218264154490)) < 1e-8,
            "Nelder-Mead HeH+ ground energy matches reference within 1e-8 Ha");
    require(hehRes.evaluations < 500, "Nelder-Mead HeH+ converges in under 500 evaluations (vs 6803 for grid search)");
  }

  // =========================================================================
  // 5. Parameterized Hardware-Efficient Ansatz & LiH Active Space VQE
  // =========================================================================
  {
    // 5a. Hardware-Efficient RY-CX Ansatz on H2 (2 Qubits, 4 Parameters)
    const QubitHamiltonian h2Ham = h2Hamiltonian();
    auto heaRes = compiler::NelderMeadOptimizer::minimize(
        [&](const std::vector<double> &params) {
          const auto circuit = compiler::createHardwareEfficientAnsatz(2, 1, params);
          return compiler::evaluateAnsatzEnergy(circuit, 2, h2Ham);
        },
        {0.0, 0.0, 0.0, 0.0}
    );
    require(heaRes.success, "Hardware-efficient ansatz optimization on H2 converges");
    const auto optHeaCircuit = compiler::createHardwareEfficientAnsatz(2, 1, heaRes.optimalParameters);
    const double heaEnergy = compiler::evaluateAnsatzEnergy(optHeaCircuit, 2, h2Ham);
    require(std::abs(heaEnergy - (-1.857275029288)) < 1e-6,
            "Hardware-efficient ansatz reaches H2 ground-state within 1e-6 Ha");

    // 5b. LiH Active-Space VQE (2 Qubits, 4 Parameters)
    MolecularGeometry lihGeom;
    lihGeom.atoms = {{"Li", 3, {0.0, 0.0, 0.0}}, {"H", 1, {0.0, 0.0, 1.595}}};
    lihGeom.unit = GeometryUnit::Angstrom;
    ActiveSpaceConfig actConfig{1, 2}; // 1 core, 2 active
    const auto lihHam = buildMolecularHamiltonian(lihGeom, MappingType::Parity, actConfig);

    compiler::OptimizerOptions opt;
    opt.maxIterations = 1000;
    opt.maxEvaluations = 3000;
    opt.tolerance = 1e-8;

    auto lihRes = compiler::NelderMeadOptimizer::minimize(
        [&](const std::vector<double> &p) {
          const auto circuit = compiler::createHardwareEfficientAnsatz(2, 1, p);
          return compiler::evaluateAnsatzEnergy(circuit, 2, lihHam.hamiltonian);
        },
        {0.0, 0.0, 0.0, 0.0},
        {},
        {},
        opt
    );
    require(lihRes.success, "LiH active-space VQE optimization converges");
    const auto optLihCircuit = compiler::createHardwareEfficientAnsatz(2, 1, lihRes.optimalParameters);
    const double lihOptCorrElec = compiler::evaluateAnsatzEnergy(optLihCircuit, 2, lihHam.hamiltonian);
    const double lihOptTotal = lihHam.coreEnergy + lihOptCorrElec;

    // Compare with PySCF CASCI(2,2) reference -7.8622855339 Ha
    require(std::abs(lihOptTotal - (-7.862285533865)) < 1e-6,
            "LiH active-space VQE total energy matches PySCF CASCI reference within 1e-6 Ha");
  }
}

} // namespace test
