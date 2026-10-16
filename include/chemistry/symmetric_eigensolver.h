#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

template <std::size_t Dimension>
double lowestSymmetricEigenvalue(
    std::array<std::array<double, Dimension>, Dimension> matrix
) {
  static_assert(Dimension > 1, "Matrix dimension must be at least two");
  // Relative symmetry bound, scaled by max(1, |a_ij|, |a_ji|).
  constexpr double relativeSymmetryTolerance = 1e-12;
  for (std::size_t row = 0; row < Dimension; ++row) {
    for (std::size_t column = 0; column < Dimension; ++column) {
      if (!std::isfinite(matrix[row][column])) {
        throw std::invalid_argument(
            "Symmetric matrix contains a non-finite value"
        );
      }
      if (column > row) {
        const double scale = std::max({
            1.0,
            std::abs(matrix[row][column]),
            std::abs(matrix[column][row])
        });
        if (std::abs(matrix[row][column] - matrix[column][row]) >
            relativeSymmetryTolerance * scale) {
          throw std::invalid_argument("Matrix is not symmetric");
        }
      }
    }
  }
  constexpr double offDiagonalTolerance = 1e-14;
  constexpr int maximumRotations = 128 * static_cast<int>(Dimension * Dimension);
  bool converged = false;

  for (int rotation = 0; rotation < maximumRotations; ++rotation) {
    std::size_t p = 0;
    std::size_t q = 1;
    double largest = 0.0;
    for (std::size_t row = 0; row < Dimension; ++row) {
      for (std::size_t column = row + 1; column < Dimension; ++column) {
        if (std::abs(matrix[row][column]) > largest) {
          largest = std::abs(matrix[row][column]);
          p = row;
          q = column;
        }
      }
    }
    if (largest <= offDiagonalTolerance) {
      converged = true;
      break;
    }

    const double angle = 0.5 * std::atan2(
        2.0 * matrix[p][q], matrix[q][q] - matrix[p][p]
    );
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double app = matrix[p][p];
    const double aqq = matrix[q][q];
    const double apq = matrix[p][q];

    for (std::size_t index = 0; index < Dimension; ++index) {
      if (index == p || index == q) {
        continue;
      }
      const double aip = matrix[index][p];
      const double aiq = matrix[index][q];
      matrix[index][p] = cosine * aip - sine * aiq;
      matrix[p][index] = matrix[index][p];
      matrix[index][q] = sine * aip + cosine * aiq;
      matrix[q][index] = matrix[index][q];
    }
    matrix[p][p] = cosine * cosine * app -
                   2.0 * sine * cosine * apq + sine * sine * aqq;
    matrix[q][q] = sine * sine * app +
                   2.0 * sine * cosine * apq + cosine * cosine * aqq;
    matrix[p][q] = 0.0;
    matrix[q][p] = 0.0;
  }
  if (!converged) {
    throw std::runtime_error("Symmetric Jacobi diagonalization did not converge");
  }

  double minimum = matrix[0][0];
  for (std::size_t index = 1; index < Dimension; ++index) {
    minimum = std::min(minimum, matrix[index][index]);
  }
  return minimum;
}
