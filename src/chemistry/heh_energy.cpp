#include "chemistry/heh_energy.h"

#include "simulator/simulator.h"
#include "chemistry/symmetric_eigensolver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <utility>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kGolden = 0.6180339887498948482;
constexpr double kParameterTolerance = 1e-8;
constexpr double kEnergyTolerance = 1e-8;
constexpr int kMaximumCoordinateIterations = 48;
constexpr int kMaximumEvaluations = 10000;
constexpr double kReferenceTolerance = 1e-10;
constexpr std::size_t kDimension = 16;

Pauli pauliFromLabel(char label) {
  switch (label) {
  case 'I': return Pauli::I;
  case 'X': return Pauli::X;
  case 'Y': return Pauli::Y;
  case 'Z': return Pauli::Z;
  default: throw std::logic_error("Invalid frozen HeH+ Pauli label");
  }
}

PauliTerm term(double coefficient, const std::string &label) {
  if (label.size() != 4) {
    throw std::logic_error("Frozen HeH+ Pauli label must have four characters");
  }
  std::vector<Pauli> operators(4);
  for (std::size_t qubit = 0; qubit < 4; ++qubit) {
    operators[qubit] = pauliFromLabel(label[3 - qubit]);
  }
  return {coefficient, std::move(operators)};
}

void appendGate(
    std::vector<Instruction> &circuit,
    const char *name,
    int qubit,
    double parameter = 0.0
) {
  circuit.push_back({InstructionType::ApplyGate, name, {qubit}, parameter});
}

void appendCnot(std::vector<Instruction> &circuit, int control, int target) {
  circuit.push_back({InstructionType::ApplyGate, "CNOT", {control, target}, 0.0});
}

void appendRY(std::vector<Instruction> &circuit, int target, double angle) {
  appendGate(circuit, "RZ", target, -kPi / 2.0);
  appendGate(circuit, "H", target);
  appendGate(circuit, "RZ", target, angle);
  appendGate(circuit, "H", target);
  appendGate(circuit, "RZ", target, kPi / 2.0);
}

void appendMultiControlledRY(
    std::vector<Instruction> &circuit,
    int target,
    const std::vector<int> &controls,
    unsigned int activePattern,
    double angle
) {
  const std::size_t count = std::size_t{1} << controls.size();
  std::vector<unsigned int> gray(count);
  for (std::size_t index = 0; index < count; ++index) {
    const unsigned int binary = static_cast<unsigned int>(index);
    gray[index] = binary ^ (binary >> 1U);
  }
  std::vector<double> patternAngles(count, 0.0);
  patternAngles[activePattern] = angle;
  for (std::size_t index = 0; index < count; ++index) {
    double rotation = 0.0;
    for (std::size_t pattern = 0; pattern < count; ++pattern) {
      unsigned int parityValue =
          gray[index] & static_cast<unsigned int>(pattern);
      unsigned int parity = 0;
      while (parityValue != 0U) {
        parity ^= parityValue & 1U;
        parityValue >>= 1U;
      }
      rotation += parity == 0U ? patternAngles[pattern] : -patternAngles[pattern];
    }
    appendRY(circuit, target, rotation / static_cast<double>(count));
    const unsigned int next = gray[(index + 1) % count];
    const unsigned int changed = gray[index] ^ next;
    unsigned int changedBits = changed;
    unsigned int controlBit = 0;
    while ((changedBits & 1U) == 0U) {
      ++controlBit;
      changedBits >>= 1U;
    }
    appendCnot(circuit, controls[controlBit], target);
  }
}

void appendGivens(
    std::vector<Instruction> &circuit,
    std::size_t first,
    std::size_t second,
    double angle
) {
  int target = -1;
  std::vector<int> changedBits;
  for (int qubit = 0; qubit < 4; ++qubit) {
    if (((first ^ second) & (std::size_t{1} << qubit)) != 0) {
      if (target < 0) target = qubit;
      changedBits.push_back(qubit);
    }
  }
  if (target < 0) {
    throw std::logic_error("Givens rotation requires distinct basis states");
  }

  std::vector<int> flippedBits;
  for (int qubit = 0; qubit < 4; ++qubit) {
    if ((first & (std::size_t{1} << qubit)) != 0) {
      appendGate(circuit, "X", qubit);
      flippedBits.push_back(qubit);
    }
  }
  std::vector<int> cxControls;
  for (int qubit : changedBits) {
    if (qubit != target) {
      appendCnot(circuit, target, qubit);
      cxControls.push_back(qubit);
    }
  }
  std::vector<int> controls;
  for (int qubit = 0; qubit < 4; ++qubit) {
    if (qubit != target) controls.push_back(qubit);
  }
  appendMultiControlledRY(circuit, target, controls, 0U, 2.0 * angle);
  for (auto iter = cxControls.rbegin(); iter != cxControls.rend(); ++iter) {
    appendCnot(circuit, target, *iter);
  }
  for (auto iter = flippedBits.rbegin(); iter != flippedBits.rend(); ++iter) {
    appendGate(circuit, "X", *iter);
  }
}

void appendSingletRotation(
    std::vector<Instruction> &circuit,
    std::size_t first,
    std::size_t second,
    double angle
) {
  // B maps |9> to (|9>+|6>)/sqrt(2); execute B†, G, B.
  appendGivens(circuit, 9, 6, -kPi / 4.0);
  appendGivens(circuit, first, second, angle);
  appendGivens(circuit, 9, 6, kPi / 4.0);
}

using Matrix16 = std::array<std::array<std::complex<double>, kDimension>, kDimension>;

Matrix16 fullHamiltonianMatrix() {
  Matrix16 matrix{};
  const QubitHamiltonian &hamiltonian = hehHamiltonian();
  for (const PauliTerm &pauliTerm : hamiltonian.terms()) {
    for (std::size_t column = 0; column < kDimension; ++column) {
      std::size_t row = column;
      std::complex<double> phase(1.0, 0.0);
      for (std::size_t qubit = 0; qubit < 4; ++qubit) {
        const std::size_t mask = std::size_t{1} << qubit;
        const bool one = (column & mask) != 0;
        switch (pauliTerm.operators[qubit]) {
        case Pauli::I: break;
        case Pauli::X: row ^= mask; break;
        case Pauli::Y:
          row ^= mask;
          phase *= one ? std::complex<double>(0.0, -1.0)
                       : std::complex<double>(0.0, 1.0);
          break;
        case Pauli::Z: if (one) phase = -phase; break;
        }
      }
      matrix[row][column] += pauliTerm.coefficient * phase;
    }
  }
  return matrix;
}

std::array<std::array<double, 3>, 3> singletMatrix() {
  const Matrix16 full = fullHamiltonianMatrix();
  const double inverseSqrtTwo = 1.0 / std::sqrt(2.0);
  const std::array<std::array<double, kDimension>, 3> basis = [&]() {
    std::array<std::array<double, kDimension>, 3> vectors{};
    vectors[0][5] = 1.0;
    vectors[1][6] = inverseSqrtTwo;
    vectors[1][9] = inverseSqrtTwo;
    vectors[2][10] = 1.0;
    return vectors;
  }();
  std::array<std::array<double, 3>, 3> projected{};
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t column = 0; column < 3; ++column) {
      std::complex<double> value(0.0, 0.0);
      for (std::size_t i = 0; i < kDimension; ++i) {
        for (std::size_t j = 0; j < kDimension; ++j) {
          value += basis[row][i] * full[i][j] * basis[column][j];
        }
      }
      if (std::abs(value.imag()) > 1e-12) {
        throw std::runtime_error("HeH+ singlet projection is not real");
      }
      projected[row][column] = value.real();
    }
  }
  return projected;
}

double wrapPhi(double phi) {
  phi = std::fmod(phi, 2.0 * kPi);
  return phi < 0.0 ? phi + 2.0 * kPi : phi;
}

} // namespace

const QubitHamiltonian &hehHamiltonian() {
  static const QubitHamiltonian hamiltonian(4, {
      term(-2.9120172734747536, "IIII"),
      term(0.7583883292075204, "IIIZ"),
      term(0.05256451032097764, "IIYY"),
      term(0.05256451032097764, "IIXX"),
      term(0.19189724869368352, "IIZI"),
      term(0.7583883292075202, "IZII"),
      term(0.05256451032097764, "YYII"),
      term(0.05256451032097764, "XXII"),
      term(0.19189724869368363, "ZIII"),
      term(0.12871109487820281, "IIZZ"),
      term(0.23577440756967472, "IZIZ"),
      term(0.04324233541007819, "YYIZ"),
      term(0.04324233541007819, "XXIZ"),
      term(0.16505925176240369, "ZIIZ"),
      term(0.04324233541007819, "IZYY"),
      term(0.04324233541007819, "IZXX"),
      term(0.03634815688420088, "YYYY"),
      term(0.03634815688420088, "XXYY"),
      term(0.03634815688420088, "YYXX"),
      term(0.03634815688420088, "XXXX"),
      term(-0.009322174907132445, "ZIYY"),
      term(-0.009322174907132445, "ZIXX"),
      term(0.16505925176240369, "IZZI"),
      term(-0.009322174907132445, "YYZI"),
      term(-0.009322174907132445, "XXZI"),
      term(0.18813101682146008, "ZIZI"),
      term(0.12871109487820281, "ZZII")
  });
  return hamiltonian;
}

std::vector<Instruction> hehSingletTrialCircuit(double theta, double phi) {
  if (!std::isfinite(theta) || !std::isfinite(phi) ||
      theta < 0.0 || theta > kPi || phi < 0.0 || phi >= 2.0 * kPi) {
    throw std::invalid_argument("HeH+ trial parameters are outside their bounds");
  }
  std::vector<Instruction> circuit;
  circuit.push_back({InstructionType::ApplyGate, "X", {0}, 0.0});
  circuit.push_back({InstructionType::ApplyGate, "X", {2}, 0.0});
  appendSingletRotation(circuit, 5, 9, theta);
  appendSingletRotation(circuit, 9, 10, phi);
  return circuit;
}

double evaluateHeHEnergy(double theta, double phi) {
  Simulator simulator(4);
  simulator.run(hehSingletTrialCircuit(theta, phi));
  return hehHamiltonian().expectation(simulator);
}

double exactHeHSingletEnergy() {
  return lowestSymmetricEigenvalue(singletMatrix());
}

HeHResult optimizeHeH() {
  const double exactEnergy = exactHeHSingletEnergy();
  if (std::abs(exactEnergy - HEH_ELECTRONIC_REFERENCE_HARTREE) >
      kReferenceTolerance) {
    throw std::runtime_error("HeH+ exact singlet energy disagrees with Psi4 FCI");
  }

  int evaluations = 0;
  auto energy = [&](double theta, double phi) {
    if (evaluations >= kMaximumEvaluations) {
      throw std::runtime_error("HeH+ optimizer exceeded its evaluation limit");
    }
    ++evaluations;
    return evaluateHeHEnergy(theta, wrapPhi(phi));
  };
  double bestTheta = 0.0;
  double bestPhi = 0.0;
  double bestEnergy = std::numeric_limits<double>::infinity();
  constexpr int thetaStarts = 3;
  constexpr int phiStarts = 4;
  const double thetaStep = kPi / 2.0;
  const double phiStep = kPi / 2.0;

  auto goldenSection = [&](double lower, double upper, double fixed, bool thetaAxis) {
    double left = upper - kGolden * (upper - lower);
    double right = lower + kGolden * (upper - lower);
    auto evaluate = [&](double value) {
      return thetaAxis ? energy(value, fixed) : energy(fixed, value);
    };
    double leftEnergy = evaluate(left);
    double rightEnergy = evaluate(right);
    int iterations = 0;
    while (upper - lower > kParameterTolerance &&
           iterations < kMaximumCoordinateIterations) {
      if (leftEnergy <= rightEnergy) {
        upper = right;
        right = left;
        rightEnergy = leftEnergy;
        left = upper - kGolden * (upper - lower);
        leftEnergy = evaluate(left);
      } else {
        lower = left;
        left = right;
        leftEnergy = rightEnergy;
        right = lower + kGolden * (upper - lower);
        rightEnergy = evaluate(right);
      }
      ++iterations;
    }
    if (upper - lower > kParameterTolerance) {
      throw std::runtime_error("HeH+ optimizer coordinate search did not converge");
    }
    const double optimum = 0.5 * (lower + upper);
    return std::pair<double, double>{optimum, evaluate(optimum)};
  };

  for (int thetaIndex = 0; thetaIndex < thetaStarts; ++thetaIndex) {
    for (int phiIndex = 0; phiIndex < phiStarts; ++phiIndex) {
      double theta = thetaIndex * thetaStep;
      double phi = phiIndex * phiStep;
      double currentEnergy = energy(theta, phi);
      for (int sweep = 0; sweep < 8; ++sweep) {
        const double previousEnergy = currentEnergy;
        auto thetaMinimum = goldenSection(
            std::max(0.0, theta - thetaStep),
            std::min(kPi, theta + thetaStep), phi, true
        );
        theta = thetaMinimum.first;
        currentEnergy = thetaMinimum.second;
        auto phiMinimum = goldenSection(
            phi - phiStep, phi + phiStep, theta, false
        );
        phi = wrapPhi(phiMinimum.first);
        currentEnergy = phiMinimum.second;
        if (std::abs(previousEnergy - currentEnergy) < 1e-13) break;
      }
      if (currentEnergy < bestEnergy) {
        bestEnergy = currentEnergy;
        bestTheta = theta;
        bestPhi = phi;
      }
    }
  }
  const double finalEnergy = energy(bestTheta, bestPhi);
  if (finalEnergy < exactEnergy - 1e-10) {
    throw std::runtime_error("HeH+ variational energy is below the exact minimum");
  }
  if (std::abs(finalEnergy - exactEnergy) > kEnergyTolerance) {
    throw std::runtime_error("HeH+ variational state failed exact-energy validation");
  }
  return {
      bestTheta,
      bestPhi,
      finalEnergy,
      exactEnergy,
      HEH_NUCLEAR_REPULSION_HARTREE,
      finalEnergy + HEH_NUCLEAR_REPULSION_HARTREE,
      evaluations
  };
}
