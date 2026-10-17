#include "chemistry/molecular_integral_engine.h"
#include "chemistry/symmetric_eigensolver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kAngstromToBohr = 1.0 / 0.52917721092;

const std::vector<GaussianPrimitive> kH_1s{
    {3.42525091, 0.15432897},
    {0.62391373, 0.53532814},
    {0.16885540, 0.44463454}
};

const std::vector<GaussianPrimitive> kHe_1s{
    {6.36242139, 0.15432897},
    {1.15892300, 0.53532814},
    {0.31364979, 0.44463454}
};

const std::vector<GaussianPrimitive> kLi_1s{
    {16.11957500, 0.15432897},
    {2.93620070, 0.53532814},
    {0.79465050, 0.44463454}
};

const std::vector<GaussianPrimitive> kLi_2s{
    {0.63628970, -0.09996723},
    {0.14786010, 0.39951283},
    {0.04808870, 0.70011547}
};

const std::vector<GaussianPrimitive> kLi_2p{
    {0.63628970, 0.15591627},
    {0.14786010, 0.60768372},
    {0.04808870, 0.39195739}
};

// General Boys function F_m(t) for m in [0, mMax]
std::vector<double> boysFunctionArray(int mMax, double t) {
  std::vector<double> res(static_cast<std::size_t>(mMax + 1), 0.0);
  if (t < 1e-7) {
    for (int m = 0; m <= mMax; ++m) {
      res[static_cast<std::size_t>(m)] =
          1.0 / (2.0 * m + 1.0) - t / (2.0 * m + 3.0) +
          0.5 * t * t / (2.0 * m + 5.0) - (t * t * t) / (6.0 * (2.0 * m + 7.0));
    }
    return res;
  }
  if (t > 25.0) {
    res[0] = 0.5 * std::sqrt(kPi / t) * std::erf(std::sqrt(t));
    const double expT = std::exp(-t);
    for (int m = 0; m < mMax; ++m) {
      res[static_cast<std::size_t>(m + 1)] =
          ((2.0 * m + 1.0) * res[static_cast<std::size_t>(m)] - expT) / (2.0 * t);
    }
    return res;
  }

  // Downward recurrence (Miller algorithm)
  const int M = std::max(mMax + 30, static_cast<int>(t) + 30);
  const double expT = std::exp(-t);
  std::vector<double> table(static_cast<std::size_t>(M + 1), 0.0);
  table[static_cast<std::size_t>(M)] = 0.0;
  for (int k = M - 1; k >= 0; --k) {
    table[static_cast<std::size_t>(k)] =
        (2.0 * t * table[static_cast<std::size_t>(k + 1)] + expT) / (2.0 * k + 1.0);
  }
  const double f0Exact = 0.5 * std::sqrt(kPi / t) * std::erf(std::sqrt(t));
  const double scale = f0Exact / table[0];
  for (int m = 0; m <= mMax; ++m) {
    res[static_cast<std::size_t>(m)] = table[static_cast<std::size_t>(m)] * scale;
  }
  return res;
}

double primitiveNormFactor(double alpha, int lx, int ly, int lz) {
  const int lTot = lx + ly + lz;
  const double pre = std::pow(2.0 * alpha / kPi, 0.75);
  if (lTot == 0) return pre;
  if (lTot == 1) return pre * 2.0 * std::sqrt(alpha);
  throw std::invalid_argument("Angular momentum > 1 unsupported");
}

double overlap1D(int la, int lb, double Pa, double Pb, double gamma) {
  std::vector<std::vector<double>> table(la + 1, std::vector<double>(lb + 1, 0.0));
  table[0][0] = std::sqrt(kPi / gamma);
  for (int i = 0; i <= la; ++i) {
    for (int j = 0; j <= lb; ++j) {
      if (i == 0 && j == 0) continue;
      if (j == 0) {
        table[i][0] = Pa * table[i - 1][0] +
                      ((i >= 2) ? ((i - 1.0) / (2.0 * gamma)) * table[i - 2][0] : 0.0);
      } else {
        table[i][j] = Pb * table[i][j - 1] +
                      ((i >= 1) ? (static_cast<double>(i) / (2.0 * gamma)) * table[i - 1][j - 1] : 0.0) +
                      ((j >= 2) ? ((j - 1.0) / (2.0 * gamma)) * table[i][j - 2] : 0.0);
      }
    }
  }
  return table[la][lb];
}

double primitiveOverlap(double a, const std::array<int, 3> &la, const std::array<double, 3> &A,
                        double b, const std::array<int, 3> &lb, const std::array<double, 3> &B) {
  const double gamma = a + b;
  std::array<double, 3> P{};
  double r2 = 0.0;
  for (int d = 0; d < 3; ++d) {
    P[d] = (a * A[d] + b * B[d]) / gamma;
    const double diff = A[d] - B[d];
    r2 += diff * diff;
  }
  const double K = std::exp(-a * b / gamma * r2);
  const double Ix = overlap1D(la[0], lb[0], P[0] - A[0], P[0] - B[0], gamma);
  const double Iy = overlap1D(la[1], lb[1], P[1] - A[1], P[1] - B[1], gamma);
  const double Iz = overlap1D(la[2], lb[2], P[2] - A[2], P[2] - B[2], gamma);
  const double normA = primitiveNormFactor(a, la[0], la[1], la[2]);
  const double normB = primitiveNormFactor(b, lb[0], lb[1], lb[2]);
  return normA * normB * K * Ix * Iy * Iz;
}

double primitiveKinetic(double a, const std::array<int, 3> &la, const std::array<double, 3> &A,
                        double b, const std::array<int, 3> &lb, const std::array<double, 3> &B) {
  const double gamma = a + b;
  std::array<double, 3> P{};
  double r2 = 0.0;
  for (int d = 0; d < 3; ++d) {
    P[d] = (a * A[d] + b * B[d]) / gamma;
    const double diff = A[d] - B[d];
    r2 += diff * diff;
  }
  const double K = std::exp(-a * b / gamma * r2);

  auto getD = [&](int la_w, int lb_w, double Pa_w, double Pb_w) {
    const double I0 = (lb_w >= 2) ? overlap1D(la_w, lb_w - 2, Pa_w, Pb_w, gamma) : 0.0;
    const double I1 = overlap1D(la_w, lb_w, Pa_w, Pb_w, gamma);
    const double I2 = overlap1D(la_w, lb_w + 2, Pa_w, Pb_w, gamma);
    const double t1 = lb_w * (lb_w - 1) * I0;
    const double t2 = -2.0 * b * (2 * lb_w + 1) * I1;
    const double t3 = 4.0 * b * b * I2;
    return -0.5 * (t1 + t2 + t3);
  };

  const double Dx = getD(la[0], lb[0], P[0] - A[0], P[0] - B[0]);
  const double Dy = getD(la[1], lb[1], P[1] - A[1], P[1] - B[1]);
  const double Dz = getD(la[2], lb[2], P[2] - A[2], P[2] - B[2]);

  const double Ix = overlap1D(la[0], lb[0], P[0] - A[0], P[0] - B[0], gamma);
  const double Iy = overlap1D(la[1], lb[1], P[1] - A[1], P[1] - B[1], gamma);
  const double Iz = overlap1D(la[2], lb[2], P[2] - A[2], P[2] - B[2], gamma);

  const double tUnnorm = K * (Dx * Iy * Iz + Ix * Dy * Iz + Ix * Iy * Dz);
  const double normA = primitiveNormFactor(a, la[0], la[1], la[2]);
  const double normB = primitiveNormFactor(b, lb[0], lb[1], lb[2]);
  return normA * normB * tUnnorm;
}

double primitiveNuclearSingleCenter(
    double a, const std::array<int, 3> &la, const std::array<double, 3> &A,
    double b, const std::array<int, 3> &lb, const std::array<double, 3> &B,
    const std::array<double, 3> &C
) {
  const double gamma = a + b;
  std::array<double, 3> P{};
  double rAB2 = 0.0;
  for (int d = 0; d < 3; ++d) {
    P[d] = (a * A[d] + b * B[d]) / gamma;
    const double diff = A[d] - B[d];
    rAB2 += diff * diff;
  }
  const double K_AB = std::exp(-a * b / gamma * rAB2);
  double rPC2 = 0.0;
  for (int d = 0; d < 3; ++d) {
    const double diff = P[d] - C[d];
    rPC2 += diff * diff;
  }
  const double T = gamma * rPC2;
  const auto F = boysFunctionArray(2, T);

  const int laTot = la[0] + la[1] + la[2];
  const int lbTot = lb[0] + lb[1] + lb[2];

  std::array<double, 3> v00{
      2.0 * kPi / gamma * K_AB * F[0],
      2.0 * kPi / gamma * K_AB * F[1],
      2.0 * kPi / gamma * K_AB * F[2]
  };

  double val = 0.0;
  if (laTot == 0 && lbTot == 0) {
    val = v00[0];
  } else if (laTot == 1 && lbTot == 0) {
    int idx = (la[0] == 1) ? 0 : ((la[1] == 1) ? 1 : 2);
    val = (P[idx] - A[idx]) * v00[0] - (P[idx] - C[idx]) * v00[1];
  } else if (laTot == 0 && lbTot == 1) {
    int idx = (lb[0] == 1) ? 0 : ((lb[1] == 1) ? 1 : 2);
    val = (P[idx] - B[idx]) * v00[0] - (P[idx] - C[idx]) * v00[1];
  } else if (laTot == 1 && lbTot == 1) {
    int i = (la[0] == 1) ? 0 : ((la[1] == 1) ? 1 : 2);
    int j = (lb[0] == 1) ? 0 : ((lb[1] == 1) ? 1 : 2);
    const double delta_ij = (i == j) ? 1.0 : 0.0;
    const double v_i0_0 = (P[i] - A[i]) * v00[0] - (P[i] - C[i]) * v00[1];
    const double v_i0_1 = (P[i] - A[i]) * v00[1] - (P[i] - C[i]) * v00[2];
    val = (P[j] - B[j]) * v_i0_0 - (P[j] - C[j]) * v_i0_1 +
          (delta_ij / (2.0 * gamma)) * (v00[0] - v00[1]);
  }

  const double normA = primitiveNormFactor(a, la[0], la[1], la[2]);
  const double normB = primitiveNormFactor(b, lb[0], lb[1], lb[2]);
  return normA * normB * val;
}

double primitiveEri(
    double a, const std::array<int, 3> &la, const std::array<double, 3> &A,
    double b, const std::array<int, 3> &lb, const std::array<double, 3> &B,
    double c, const std::array<int, 3> &lc, const std::array<double, 3> &C,
    double d, const std::array<int, 3> &ld, const std::array<double, 3> &D
) {
  const double g1 = a + b;
  const double g2 = c + d;
  const double rho = g1 * g2 / (g1 + g2);

  std::array<double, 3> P{};
  std::array<double, 3> Q{};
  std::array<double, 3> W{};
  double rAB2 = 0.0;
  double rCD2 = 0.0;
  double rPQ2 = 0.0;

  for (int i = 0; i < 3; ++i) {
    P[i] = (a * A[i] + b * B[i]) / g1;
    Q[i] = (c * C[i] + d * D[i]) / g2;
    W[i] = (g1 * P[i] + g2 * Q[i]) / (g1 + g2);
    const double diffAB = A[i] - B[i];
    rAB2 += diffAB * diffAB;
    const double diffCD = C[i] - D[i];
    rCD2 += diffCD * diffCD;
    const double diffPQ = P[i] - Q[i];
    rPQ2 += diffPQ * diffPQ;
  }

  const double K_AB = std::exp(-a * b / g1 * rAB2);
  const double K_CD = std::exp(-c * d / g2 * rCD2);
  const double T = rho * rPQ2;

  const int mTot = (la[0] + la[1] + la[2]) + (lb[0] + lb[1] + lb[2]) +
                   (lc[0] + lc[1] + lc[2]) + (ld[0] + ld[1] + ld[2]);
  const auto F = boysFunctionArray(mTot, T);
  const double pre = 2.0 * std::pow(kPi, 2.5) / (g1 * g2 * std::sqrt(g1 + g2)) * K_AB * K_CD;

  // General Obara-Saika vertical recurrence (VRR) for l <= 1
  auto vrr = [&](auto &self, std::array<int, 3> cur_la, std::array<int, 3> cur_lb,
                 std::array<int, 3> cur_lc, std::array<int, 3> cur_ld, int m) -> double {
    if (cur_la == std::array<int, 3>{0, 0, 0} && cur_lb == std::array<int, 3>{0, 0, 0} &&
        cur_lc == std::array<int, 3>{0, 0, 0} && cur_ld == std::array<int, 3>{0, 0, 0}) {
      return pre * F[static_cast<std::size_t>(m)];
    }

    // VRR on center A if any l_A > 0
    for (int i = 0; i < 3; ++i) {
      if (cur_la[i] > 0) {
        auto la_minus = cur_la; la_minus[i] -= 1;
        const double t1 = (P[i] - A[i]) * self(self, la_minus, cur_lb, cur_lc, cur_ld, m);
        const double t2 = (W[i] - P[i]) * self(self, la_minus, cur_lb, cur_lc, cur_ld, m + 1);
        double t3 = 0.0;
        if (la_minus[i] > 0) {
          auto la_minus2 = la_minus; la_minus2[i] -= 1;
          t3 = (la_minus[i] / (2.0 * g1)) *
               (self(self, la_minus2, cur_lb, cur_lc, cur_ld, m) - (rho / g1) * self(self, la_minus2, cur_lb, cur_lc, cur_ld, m + 1));
        }
        double t4 = 0.0;
        if (cur_lb[i] > 0) {
          auto lb_minus = cur_lb; lb_minus[i] -= 1;
          t4 = (cur_lb[i] / (2.0 * g1)) *
               (self(self, la_minus, lb_minus, cur_lc, cur_ld, m) - (rho / g1) * self(self, la_minus, lb_minus, cur_lc, cur_ld, m + 1));
        }
        double t5 = 0.0;
        if (cur_lc[i] > 0) {
          auto lc_minus = cur_lc; lc_minus[i] -= 1;
          t5 = (cur_lc[i] / (2.0 * (g1 + g2))) * self(self, la_minus, cur_lb, lc_minus, cur_ld, m + 1);
        }
        double t6 = 0.0;
        if (cur_ld[i] > 0) {
          auto ld_minus = cur_ld; ld_minus[i] -= 1;
          t6 = (cur_ld[i] / (2.0 * (g1 + g2))) * self(self, la_minus, cur_lb, cur_lc, ld_minus, m + 1);
        }
        return t1 + t2 + t3 + t4 + t5 + t6;
      }
    }

    // VRR on center B if any l_B > 0
    for (int i = 0; i < 3; ++i) {
      if (cur_lb[i] > 0) {
        auto lb_minus = cur_lb; lb_minus[i] -= 1;
        const double t1 = (P[i] - B[i]) * self(self, cur_la, lb_minus, cur_lc, cur_ld, m);
        const double t2 = (W[i] - P[i]) * self(self, cur_la, lb_minus, cur_lc, cur_ld, m + 1);
        double t3 = 0.0;
        if (lb_minus[i] > 0) {
          auto lb_minus2 = lb_minus; lb_minus2[i] -= 1;
          t3 = (lb_minus[i] / (2.0 * g1)) *
               (self(self, cur_la, lb_minus2, cur_lc, cur_ld, m) - (rho / g1) * self(self, cur_la, lb_minus2, cur_lc, cur_ld, m + 1));
        }
        double t4 = 0.0;
        if (cur_la[i] > 0) {
          auto la_minus = cur_la; la_minus[i] -= 1;
          t4 = (cur_la[i] / (2.0 * g1)) *
               (self(self, la_minus, lb_minus, cur_lc, cur_ld, m) - (rho / g1) * self(self, la_minus, lb_minus, cur_lc, cur_ld, m + 1));
        }
        double t5 = 0.0;
        if (cur_lc[i] > 0) {
          auto lc_minus = cur_lc; lc_minus[i] -= 1;
          t5 = (cur_lc[i] / (2.0 * (g1 + g2))) * self(self, cur_la, lb_minus, lc_minus, cur_ld, m + 1);
        }
        double t6 = 0.0;
        if (cur_ld[i] > 0) {
          auto ld_minus = cur_ld; ld_minus[i] -= 1;
          t6 = (cur_ld[i] / (2.0 * (g1 + g2))) * self(self, cur_la, lb_minus, cur_lc, ld_minus, m + 1);
        }
        return t1 + t2 + t3 + t4 + t5 + t6;
      }
    }

    // VRR on center C if any l_C > 0
    for (int i = 0; i < 3; ++i) {
      if (cur_lc[i] > 0) {
        auto lc_minus = cur_lc; lc_minus[i] -= 1;
        const double t1 = (Q[i] - C[i]) * self(self, cur_la, cur_lb, lc_minus, cur_ld, m);
        const double t2 = (W[i] - Q[i]) * self(self, cur_la, cur_lb, lc_minus, cur_ld, m + 1);
        double t3 = 0.0;
        if (lc_minus[i] > 0) {
          auto lc_minus2 = lc_minus; lc_minus2[i] -= 1;
          t3 = (lc_minus[i] / (2.0 * g2)) *
               (self(self, cur_la, cur_lb, lc_minus2, cur_ld, m) - (rho / g2) * self(self, cur_la, cur_lb, lc_minus2, cur_ld, m + 1));
        }
        double t4 = 0.0;
        if (cur_ld[i] > 0) {
          auto ld_minus = cur_ld; ld_minus[i] -= 1;
          t4 = (cur_ld[i] / (2.0 * g2)) *
               (self(self, cur_la, cur_lb, lc_minus, ld_minus, m) - (rho / g2) * self(self, cur_la, cur_lb, lc_minus, ld_minus, m + 1));
        }
        double t5 = 0.0;
        if (cur_la[i] > 0) {
          auto la_minus = cur_la; la_minus[i] -= 1;
          t5 = (cur_la[i] / (2.0 * (g1 + g2))) * self(self, la_minus, cur_lb, lc_minus, cur_ld, m + 1);
        }
        double t6 = 0.0;
        if (cur_lb[i] > 0) {
          auto lb_minus = cur_lb; lb_minus[i] -= 1;
          t6 = (cur_lb[i] / (2.0 * (g1 + g2))) * self(self, cur_la, lb_minus, lc_minus, cur_ld, m + 1);
        }
        return t1 + t2 + t3 + t4 + t5 + t6;
      }
    }

    // VRR on center D if any l_D > 0
    for (int i = 0; i < 3; ++i) {
      if (cur_ld[i] > 0) {
        auto ld_minus = cur_ld; ld_minus[i] -= 1;
        const double t1 = (Q[i] - D[i]) * self(self, cur_la, cur_lb, cur_lc, ld_minus, m);
        const double t2 = (W[i] - Q[i]) * self(self, cur_la, cur_lb, cur_lc, ld_minus, m + 1);
        double t3 = 0.0;
        if (ld_minus[i] > 0) {
          auto ld_minus2 = ld_minus; ld_minus2[i] -= 1;
          t3 = (ld_minus[i] / (2.0 * g2)) *
               (self(self, cur_la, cur_lb, cur_lc, ld_minus2, m) - (rho / g2) * self(self, cur_la, cur_lb, ld_minus2, cur_ld, m + 1));
        }
        double t4 = 0.0;
        if (cur_lc[i] > 0) {
          auto lc_minus = cur_lc; lc_minus[i] -= 1;
          t4 = (cur_lc[i] / (2.0 * g2)) *
               (self(self, cur_la, cur_lb, lc_minus, ld_minus, m) - (rho / g2) * self(self, cur_la, cur_lb, lc_minus, ld_minus, m + 1));
        }
        double t5 = 0.0;
        if (cur_la[i] > 0) {
          auto la_minus = cur_la; la_minus[i] -= 1;
          t5 = (cur_la[i] / (2.0 * (g1 + g2))) * self(self, la_minus, cur_lb, cur_lc, ld_minus, m + 1);
        }
        double t6 = 0.0;
        if (cur_lb[i] > 0) {
          auto lb_minus = cur_lb; lb_minus[i] -= 1;
          t6 = (cur_lb[i] / (2.0 * (g1 + g2))) * self(self, cur_la, lb_minus, cur_lc, ld_minus, m + 1);
        }
        return t1 + t2 + t3 + t4 + t5 + t6;
      }
    }

    return 0.0;
  };

  const double raw = vrr(vrr, la, lb, lc, ld, 0);
  const double norm = primitiveNormFactor(a, la[0], la[1], la[2]) *
                      primitiveNormFactor(b, lb[0], lb[1], lb[2]) *
                      primitiveNormFactor(c, lc[0], lc[1], lc[2]) *
                      primitiveNormFactor(d, ld[0], ld[1], ld[2]);
  return norm * raw;
}

// Dynamic Jacobi diagonalization for small real symmetric matrices
void symmetricJacobiDiagonalize(
    const std::vector<double> &mat,
    int n,
    std::vector<double> &eigenvalues,
    std::vector<double> &eigenvectors
) {
  eigenvalues.resize(static_cast<std::size_t>(n));
  eigenvectors.assign(static_cast<std::size_t>(n * n), 0.0);
  for (int i = 0; i < n; ++i) {
    eigenvectors[static_cast<std::size_t>(i * n + i)] = 1.0;
  }
  std::vector<double> a = mat;

  constexpr double offDiagTol = 1e-14;
  const int maxRotations = 128 * n * n;
  bool converged = false;

  for (int rot = 0; rot < maxRotations; ++rot) {
    int p = 0;
    int q = 1;
    double maxVal = 0.0;
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        double val = std::abs(a[static_cast<std::size_t>(i * n + j)]);
        if (val > maxVal) {
          maxVal = val;
          p = i;
          q = j;
        }
      }
    }

    if (maxVal < offDiagTol) {
      converged = true;
      break;
    }

    const double app = a[static_cast<std::size_t>(p * n + p)];
    const double aqq = a[static_cast<std::size_t>(q * n + q)];
    const double apq = a[static_cast<std::size_t>(p * n + q)];

    double theta = 0.5 * (aqq - app) / apq;
    double t = 0.0;
    if (theta >= 0.0) {
      t = 1.0 / (theta + std::sqrt(1.0 + theta * theta));
    } else {
      t = -1.0 / (-theta + std::sqrt(1.0 + theta * theta));
    }
    const double c = 1.0 / std::sqrt(1.0 + t * t);
    const double s = t * c;
    const double tau = s / (1.0 + c);

    a[static_cast<std::size_t>(p * n + p)] -= t * apq;
    a[static_cast<std::size_t>(q * n + q)] += t * apq;
    a[static_cast<std::size_t>(p * n + q)] = 0.0;
    a[static_cast<std::size_t>(q * n + p)] = 0.0;

    for (int k = 0; k < n; ++k) {
      if (k != p && k != q) {
        const double akp = a[static_cast<std::size_t>(k * n + p)];
        const double akq = a[static_cast<std::size_t>(k * n + q)];
        a[static_cast<std::size_t>(k * n + p)] = a[static_cast<std::size_t>(p * n + k)] = akp - s * (akq + tau * akp);
        a[static_cast<std::size_t>(k * n + q)] = a[static_cast<std::size_t>(q * n + k)] = akq + s * (akp - tau * akq);
      }
    }

    for (int k = 0; k < n; ++k) {
      const double vkp = eigenvectors[static_cast<std::size_t>(k * n + p)];
      const double vkq = eigenvectors[static_cast<std::size_t>(k * n + q)];
      eigenvectors[static_cast<std::size_t>(k * n + p)] = vkp - s * (vkq + tau * vkp);
      eigenvectors[static_cast<std::size_t>(k * n + q)] = vkq + s * (vkp - tau * vkq);
    }
  }

  if (!converged) {
    throw std::runtime_error("Dynamic Jacobi eigensolver did not converge within rotation limit");
  }

  for (int i = 0; i < n; ++i) {
    eigenvalues[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i * n + i)];
  }

  // Sort ascending by eigenvalue
  for (int i = 0; i < n - 1; ++i) {
    int minIdx = i;
    for (int j = i + 1; j < n; ++j) {
      if (eigenvalues[static_cast<std::size_t>(j)] < eigenvalues[static_cast<std::size_t>(minIdx)]) {
        minIdx = j;
      }
    }
    if (minIdx != i) {
      std::swap(eigenvalues[static_cast<std::size_t>(i)], eigenvalues[static_cast<std::size_t>(minIdx)]);
      for (int k = 0; k < n; ++k) {
        std::swap(eigenvectors[static_cast<std::size_t>(k * n + i)], eigenvectors[static_cast<std::size_t>(k * n + minIdx)]);
      }
    }
  }
}

} // namespace

void MolecularGeometry::validate() const {
  if (atoms.empty()) {
    throw std::invalid_argument("Molecular geometry contains no atoms");
  }

  int totalZ = 0;
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const auto &atom = atoms[i];
    if (atom.element != "H" && atom.element != "He" && atom.element != "Li" && atom.element != "Be") {
      throw std::invalid_argument("Unsupported element: " + atom.element + " (STO-3G integral engine supports H, He, Li, Be)");
    }
    for (int d = 0; d < 3; ++d) {
      if (!std::isfinite(atom.coordinates[static_cast<std::size_t>(d)])) {
        throw std::invalid_argument("Atomic coordinates must be finite numbers");
      }
    }
    if (atom.element == "H") totalZ += 1;
    else if (atom.element == "He") totalZ += 2;
    else if (atom.element == "Li") totalZ += 3;
    else if (atom.element == "Be") totalZ += 4;

    for (std::size_t j = i + 1; j < atoms.size(); ++j) {
      double r2 = 0.0;
      for (int d = 0; d < 3; ++d) {
        double diff = atom.coordinates[static_cast<std::size_t>(d)] - atoms[j].coordinates[static_cast<std::size_t>(d)];
        r2 += diff * diff;
      }
      if (r2 < 1e-8) {
        throw std::invalid_argument("Atoms " + std::to_string(i) + " and " + std::to_string(j) + " are coincident");
      }
    }
  }

  const int totalElectrons = totalZ - charge;
  if (totalElectrons <= 0) {
    throw std::invalid_argument("Molecular system must have at least one electron (got " + std::to_string(totalElectrons) + ")");
  }
  if (multiplicity != 1) {
    throw std::invalid_argument("Currently only closed-shell singlet states (multiplicity=1) are supported");
  }
}

std::vector<MolecularAtom> MolecularGeometry::inBohr() const {
  validate();
  std::vector<MolecularAtom> res = atoms;
  const double scale = (unit == GeometryUnit::Angstrom) ? kAngstromToBohr : 1.0;
  for (auto &atom : res) {
    if (atom.element == "H") atom.atomicNumber = 1;
    else if (atom.element == "He") atom.atomicNumber = 2;
    else if (atom.element == "Li") atom.atomicNumber = 3;
    else if (atom.element == "Be") atom.atomicNumber = 4;
    for (int d = 0; d < 3; ++d) {
      atom.coordinates[static_cast<std::size_t>(d)] *= scale;
    }
  }
  return res;
}

double MolecularGeometry::nuclearRepulsionEnergy() const {
  const auto bohrAtoms = inBohr();
  double eNuc = 0.0;
  for (std::size_t i = 0; i < bohrAtoms.size(); ++i) {
    const double zi = static_cast<double>(bohrAtoms[i].atomicNumber);
    for (std::size_t j = i + 1; j < bohrAtoms.size(); ++j) {
      const double zj = static_cast<double>(bohrAtoms[j].atomicNumber);
      double r2 = 0.0;
      for (int d = 0; d < 3; ++d) {
        double diff = bohrAtoms[i].coordinates[static_cast<std::size_t>(d)] - bohrAtoms[j].coordinates[static_cast<std::size_t>(d)];
        r2 += diff * diff;
      }
      eNuc += (zi * zj) / std::sqrt(r2);
    }
  }
  return eNuc;
}

Sto3gBasis buildSto3gBasis(const MolecularGeometry &geometry) {
  const auto bohrAtoms = geometry.inBohr();
  Sto3gBasis basis;
  for (std::size_t i = 0; i < bohrAtoms.size(); ++i) {
    const auto &atom = bohrAtoms[i];
    if (atom.element == "H" || atom.element == "He") {
      ContractedGaussian cg;
      cg.atomIndex = static_cast<int>(i);
      cg.angularMomentum = 0;
      cg.lx = cg.ly = cg.lz = 0;
      cg.centerBohr = atom.coordinates;
      cg.primitives = (atom.element == "He") ? kHe_1s : kH_1s;
      basis.basisFunctions.push_back(std::move(cg));
    } else if (atom.element == "Li") {
      // 1s
      ContractedGaussian cg1s;
      cg1s.atomIndex = static_cast<int>(i);
      cg1s.angularMomentum = 0;
      cg1s.lx = cg1s.ly = cg1s.lz = 0;
      cg1s.centerBohr = atom.coordinates;
      cg1s.primitives = kLi_1s;
      basis.basisFunctions.push_back(std::move(cg1s));

      // 2s
      ContractedGaussian cg2s;
      cg2s.atomIndex = static_cast<int>(i);
      cg2s.angularMomentum = 0;
      cg2s.lx = cg2s.ly = cg2s.lz = 0;
      cg2s.centerBohr = atom.coordinates;
      cg2s.primitives = kLi_2s;
      basis.basisFunctions.push_back(std::move(cg2s));

      // 2px
      ContractedGaussian cg2px;
      cg2px.atomIndex = static_cast<int>(i);
      cg2px.angularMomentum = 1;
      cg2px.lx = 1; cg2px.ly = 0; cg2px.lz = 0;
      cg2px.centerBohr = atom.coordinates;
      cg2px.primitives = kLi_2p;
      basis.basisFunctions.push_back(std::move(cg2px));

      // 2py
      ContractedGaussian cg2py;
      cg2py.atomIndex = static_cast<int>(i);
      cg2py.angularMomentum = 1;
      cg2py.lx = 0; cg2py.ly = 1; cg2py.lz = 0;
      cg2py.centerBohr = atom.coordinates;
      cg2py.primitives = kLi_2p;
      basis.basisFunctions.push_back(std::move(cg2py));

      // 2pz
      ContractedGaussian cg2pz;
      cg2pz.atomIndex = static_cast<int>(i);
      cg2pz.angularMomentum = 1;
      cg2pz.lx = 0; cg2pz.ly = 0; cg2pz.lz = 1;
      cg2pz.centerBohr = atom.coordinates;
      cg2pz.primitives = kLi_2p;
      basis.basisFunctions.push_back(std::move(cg2pz));
    }
  }

  // Renormalize contracted basis functions by self-overlap: <chi | chi> = 1.0
  for (auto &cg : basis.basisFunctions) {
    std::array<int, 3> l{cg.lx, cg.ly, cg.lz};
    double selfS = 0.0;
    for (const auto &p1 : cg.primitives) {
      for (const auto &p2 : cg.primitives) {
        selfS += p1.contractionCoefficient * p2.contractionCoefficient *
                 primitiveOverlap(p1.exponent, l, cg.centerBohr, p2.exponent, l, cg.centerBohr);
      }
    }
    if (selfS > 1e-12) {
      const double renorm = 1.0 / std::sqrt(selfS);
      for (auto &p : cg.primitives) {
        p.contractionCoefficient *= renorm;
      }
    }
  }

  return basis;
}

MolecularIntegrals computeMolecularIntegrals(const MolecularGeometry &geometry) {
  geometry.validate();
  const auto bohrAtoms = geometry.inBohr();
  const auto basis = buildSto3gBasis(geometry);

  const int nao = static_cast<int>(basis.basisFunctions.size());
  MolecularIntegrals res;
  res.nao = nao;
  res.nuclearRepulsion = geometry.nuclearRepulsionEnergy();

  const std::size_t n2 = static_cast<std::size_t>(nao * nao);
  const std::size_t n4 = n2 * n2;
  res.overlap.assign(n2, 0.0);
  res.kinetic.assign(n2, 0.0);
  res.nuclear.assign(n2, 0.0);
  res.hCore.assign(n2, 0.0);
  res.eri.assign(n4, 0.0);

  // 1. One-electron integrals: S, T, V
  for (int i = 0; i < nao; ++i) {
    const auto &bfI = basis.basisFunctions[static_cast<std::size_t>(i)];
    std::array<int, 3> la{bfI.lx, bfI.ly, bfI.lz};

    for (int j = i; j < nao; ++j) {
      const auto &bfJ = basis.basisFunctions[static_cast<std::size_t>(j)];
      std::array<int, 3> lb{bfJ.lx, bfJ.ly, bfJ.lz};

      double sVal = 0.0;
      double tVal = 0.0;
      double vVal = 0.0;

      for (const auto &pI : bfI.primitives) {
        for (const auto &pJ : bfJ.primitives) {
          const double cProd = pI.contractionCoefficient * pJ.contractionCoefficient;
          sVal += cProd * primitiveOverlap(pI.exponent, la, bfI.centerBohr, pJ.exponent, lb, bfJ.centerBohr);
          tVal += cProd * primitiveKinetic(pI.exponent, la, bfI.centerBohr, pJ.exponent, lb, bfJ.centerBohr);

          for (const auto &nuc : bohrAtoms) {
            vVal += -static_cast<double>(nuc.atomicNumber) * cProd *
                    primitiveNuclearSingleCenter(pI.exponent, la, bfI.centerBohr,
                                                pJ.exponent, lb, bfJ.centerBohr,
                                                nuc.coordinates);
          }
        }
      }

      res.overlap[static_cast<std::size_t>(i * nao + j)] = res.overlap[static_cast<std::size_t>(j * nao + i)] = sVal;
      res.kinetic[static_cast<std::size_t>(i * nao + j)] = res.kinetic[static_cast<std::size_t>(j * nao + i)] = tVal;
      res.nuclear[static_cast<std::size_t>(i * nao + j)] = res.nuclear[static_cast<std::size_t>(j * nao + i)] = vVal;
      res.hCore[static_cast<std::size_t>(i * nao + j)] = res.hCore[static_cast<std::size_t>(j * nao + i)] = (tVal + vVal);
    }
  }

  // 2. Two-electron integrals: (ij|kl) in chemist notation with 8-fold canonical symmetry
  for (int i = 0; i < nao; ++i) {
    const auto &bfI = basis.basisFunctions[static_cast<std::size_t>(i)];
    std::array<int, 3> la{bfI.lx, bfI.ly, bfI.lz};

    for (int j = 0; j <= i; ++j) {
      const auto &bfJ = basis.basisFunctions[static_cast<std::size_t>(j)];
      std::array<int, 3> lb{bfJ.lx, bfJ.ly, bfJ.lz};
      const int ij = i * (i + 1) / 2 + j;

      for (int k = 0; k < nao; ++k) {
        const auto &bfK = basis.basisFunctions[static_cast<std::size_t>(k)];
        std::array<int, 3> lc{bfK.lx, bfK.ly, bfK.lz};

        for (int l = 0; l <= k; ++l) {
          const auto &bfL = basis.basisFunctions[static_cast<std::size_t>(l)];
          std::array<int, 3> ld{bfL.lx, bfL.ly, bfL.lz};
          const int kl = k * (k + 1) / 2 + l;

          if (ij < kl) continue; // 8-fold canonical loop

          double eriVal = 0.0;
          for (const auto &pI : bfI.primitives) {
            for (const auto &pJ : bfJ.primitives) {
              const double cIJ = pI.contractionCoefficient * pJ.contractionCoefficient;
              for (const auto &pK : bfK.primitives) {
                for (const auto &pL : bfL.primitives) {
                  const double cTotal = cIJ * pK.contractionCoefficient * pL.contractionCoefficient;
                  eriVal += cTotal * primitiveEri(pI.exponent, la, bfI.centerBohr,
                                                  pJ.exponent, lb, bfJ.centerBohr,
                                                  pK.exponent, lc, bfK.centerBohr,
                                                  pL.exponent, ld, bfL.centerBohr);
                }
              }
            }
          }

          auto setEri = [&](int u, int v, int w, int x) {
            res.eri[static_cast<std::size_t>(((u * nao + v) * nao + w) * nao + x)] = eriVal;
          };

          // Assign 8-fold permutation symmetry
          setEri(i, j, k, l);
          setEri(j, i, k, l);
          setEri(i, j, l, k);
          setEri(j, i, l, k);
          setEri(k, l, i, j);
          setEri(l, k, i, j);
          setEri(k, l, j, i);
          setEri(l, k, j, i);
        }
      }
    }
  }

  return res;
}

RhfResult computeRhfScf(const MolecularGeometry &geometry, const MolecularIntegrals &integrals) {
  const int nao = integrals.nao;
  if (nao <= 0) {
    throw std::invalid_argument("RHF SCF requires positive number of basis functions");
  }

  int totalZ = 0;
  for (const auto &atom : geometry.atoms) {
    if (atom.element == "H") totalZ += 1;
    else if (atom.element == "He") totalZ += 2;
    else if (atom.element == "Li") totalZ += 3;
    else if (atom.element == "Be") totalZ += 4;
  }
  const int nElectrons = totalZ - geometry.charge;
  if (nElectrons <= 0 || nElectrons % 2 != 0) {
    throw std::invalid_argument("Closed-shell RHF requires a positive even number of electrons");
  }
  const int nOcc = nElectrons / 2;

  // Orthogonalizer X = S^(-1/2) via Jacobi eigensolver
  std::vector<double> sVals;
  std::vector<double> sVecs;
  symmetricJacobiDiagonalize(integrals.overlap, nao, sVals, sVecs);

  std::vector<double> X(static_cast<std::size_t>(nao * nao), 0.0);
  for (int i = 0; i < nao; ++i) {
    for (int j = 0; j < nao; ++j) {
      double sum = 0.0;
      for (int k = 0; k < nao; ++k) {
        sum += sVecs[static_cast<std::size_t>(i * nao + k)] *
               (1.0 / std::sqrt(sVals[static_cast<std::size_t>(k)])) *
               sVecs[static_cast<std::size_t>(j * nao + k)];
      }
      X[static_cast<std::size_t>(i * nao + j)] = sum;
    }
  }

  // Initial core-Hamiltonian guess
  std::vector<double> P(static_cast<std::size_t>(nao * nao), 0.0);
  std::vector<double> C(static_cast<std::size_t>(nao * nao), 0.0);
  std::vector<double> F(static_cast<std::size_t>(nao * nao), 0.0);
  std::vector<double> moE(static_cast<std::size_t>(nao), 0.0);

  bool converged = false;
  int iter = 0;
  constexpr int maxIter = 100;
  constexpr double convTol = 1e-12;

  for (iter = 0; iter < maxIter; ++iter) {
    // Form Fock matrix: F_uv = Hcore_uv + sum_ls P_ls [ (uv|ls) - 0.5 (us|lv) ]
    for (int u = 0; u < nao; ++u) {
      for (int v = 0; v < nao; ++v) {
        double fVal = integrals.getHCore(u, v);
        for (int l = 0; l < nao; ++l) {
          for (int s = 0; s < nao; ++s) {
            const double pls = P[static_cast<std::size_t>(l * nao + s)];
            fVal += pls * (integrals.getEri(u, v, l, s) - 0.5 * integrals.getEri(u, s, l, v));
          }
        }
        F[static_cast<std::size_t>(u * nao + v)] = fVal;
      }
    }

    // F' = X^T F X
    std::vector<double> FX(static_cast<std::size_t>(nao * nao), 0.0);
    for (int i = 0; i < nao; ++i) {
      for (int j = 0; j < nao; ++j) {
        double sum = 0.0;
        for (int k = 0; k < nao; ++k) {
          sum += F[static_cast<std::size_t>(i * nao + k)] * X[static_cast<std::size_t>(k * nao + j)];
        }
        FX[static_cast<std::size_t>(i * nao + j)] = sum;
      }
    }

    std::vector<double> Fprime(static_cast<std::size_t>(nao * nao), 0.0);
    for (int i = 0; i < nao; ++i) {
      for (int j = 0; j < nao; ++j) {
        double sum = 0.0;
        for (int k = 0; k < nao; ++k) {
          sum += X[static_cast<std::size_t>(k * nao + i)] * FX[static_cast<std::size_t>(k * nao + j)];
        }
        Fprime[static_cast<std::size_t>(i * nao + j)] = sum;
      }
    }

    // Diagonalize F' -> C' and mo energies
    std::vector<double> Cprime;
    symmetricJacobiDiagonalize(Fprime, nao, moE, Cprime);

    // C = X C'
    std::vector<double> Cnew(static_cast<std::size_t>(nao * nao), 0.0);
    for (int i = 0; i < nao; ++i) {
      for (int j = 0; j < nao; ++j) {
        double sum = 0.0;
        for (int k = 0; k < nao; ++k) {
          sum += X[static_cast<std::size_t>(i * nao + k)] * Cprime[static_cast<std::size_t>(k * nao + j)];
        }
        Cnew[static_cast<std::size_t>(i * nao + j)] = sum;
      }
    }

    // Update density matrix P_uv = 2 sum_{i in occ} C_ui C_vi
    std::vector<double> Pnew(static_cast<std::size_t>(nao * nao), 0.0);
    double maxDelta = 0.0;
    for (int u = 0; u < nao; ++u) {
      for (int v = 0; v < nao; ++v) {
        double pVal = 0.0;
        for (int i = 0; i < nOcc; ++i) {
          pVal += 2.0 * Cnew[static_cast<std::size_t>(u * nao + i)] * Cnew[static_cast<std::size_t>(v * nao + i)];
        }
        Pnew[static_cast<std::size_t>(u * nao + v)] = pVal;
        maxDelta = std::max(maxDelta, std::abs(pVal - P[static_cast<std::size_t>(u * nao + v)]));
      }
    }

    P = Pnew;
    C = Cnew;
    if (maxDelta < convTol) {
      converged = true;
      break;
    }
  }

  // Electronic Energy = 0.5 * sum_uv P_uv (Hcore_uv + F_uv)
  double eElec = 0.0;
  for (int u = 0; u < nao; ++u) {
    for (int v = 0; v < nao; ++v) {
      eElec += 0.5 * P[static_cast<std::size_t>(u * nao + v)] *
               (integrals.getHCore(u, v) + F[static_cast<std::size_t>(u * nao + v)]);
    }
  }

  RhfResult res;
  res.converged = converged;
  res.iterations = iter + 1;
  res.electronicEnergy = eElec;
  res.totalEnergy = eElec + integrals.nuclearRepulsion;
  res.nElectrons = nElectrons;
  res.nOccupied = nOcc;
  res.nOrbitals = nao;
  res.moCoefficients = C;
  res.moEnergies = moE;

  // 1-electron transformation: h1Mo_pq = sum_uv C_up Hcore_uv C_vq
  const std::size_t nmo2 = static_cast<std::size_t>(nao * nao);
  const std::size_t nmo4 = nmo2 * nmo2;
  res.h1Mo.assign(nmo2, 0.0);
  for (int p = 0; p < nao; ++p) {
    for (int q = 0; q < nao; ++q) {
      double sum = 0.0;
      for (int u = 0; u < nao; ++u) {
        for (int v = 0; v < nao; ++v) {
          sum += C[static_cast<std::size_t>(u * nao + p)] * integrals.getHCore(u, v) * C[static_cast<std::size_t>(v * nao + q)];
        }
      }
      res.h1Mo[static_cast<std::size_t>(p * nao + q)] = sum;
    }
  }

  // 4-index transformation to MO ERI: (pq|rs) using quarter-transformations O(N^5)
  std::vector<double> temp1(nmo4, 0.0);
  for (int p = 0; p < nao; ++p) {
    for (int v = 0; v < nao; ++v) {
      for (int w = 0; w < nao; ++w) {
        for (int x = 0; x < nao; ++x) {
          double sum = 0.0;
          for (int u = 0; u < nao; ++u) {
            sum += C[static_cast<std::size_t>(u * nao + p)] * integrals.getEri(u, v, w, x);
          }
          temp1[static_cast<std::size_t>(((p * nao + v) * nao + w) * nao + x)] = sum;
        }
      }
    }
  }

  std::vector<double> temp2(nmo4, 0.0);
  for (int p = 0; p < nao; ++p) {
    for (int q = 0; q < nao; ++q) {
      for (int w = 0; w < nao; ++w) {
        for (int x = 0; x < nao; ++x) {
          double sum = 0.0;
          for (int v = 0; v < nao; ++v) {
            sum += C[static_cast<std::size_t>(v * nao + q)] * temp1[static_cast<std::size_t>(((p * nao + v) * nao + w) * nao + x)];
          }
          temp2[static_cast<std::size_t>(((p * nao + q) * nao + w) * nao + x)] = sum;
        }
      }
    }
  }

  std::vector<double> temp3(nmo4, 0.0);
  for (int p = 0; p < nao; ++p) {
    for (int q = 0; q < nao; ++q) {
      for (int r = 0; r < nao; ++r) {
        for (int x = 0; x < nao; ++x) {
          double sum = 0.0;
          for (int w = 0; w < nao; ++w) {
            sum += C[static_cast<std::size_t>(w * nao + r)] * temp2[static_cast<std::size_t>(((p * nao + q) * nao + w) * nao + x)];
          }
          temp3[static_cast<std::size_t>(((p * nao + q) * nao + r) * nao + x)] = sum;
        }
      }
    }
  }

  res.eriMo.assign(nmo4, 0.0);
  for (int p = 0; p < nao; ++p) {
    for (int q = 0; q < nao; ++q) {
      for (int r = 0; r < nao; ++r) {
        for (int s = 0; s < nao; ++s) {
          double sum = 0.0;
          for (int x = 0; x < nao; ++x) {
            sum += C[static_cast<std::size_t>(x * nao + s)] * temp3[static_cast<std::size_t>(((p * nao + q) * nao + r) * nao + x)];
          }
          res.eriMo[static_cast<std::size_t>(((p * nao + q) * nao + r) * nao + s)] = sum;
        }
      }
    }
  }

  return res;
}

namespace {

// Ladder operator matrix builder on 4 modes (16 basis states)
void applyLadderTerm(
    std::array<std::array<std::complex<double>, 16>, 16> &mat,
    const std::vector<std::pair<char, int>> &ops,
    double coeff
) {
  for (unsigned int ket = 0; ket < 16; ++ket) {
    unsigned int state = ket;
    double phase = 1.0;
    bool valid = true;

    for (auto it = ops.rbegin(); it != ops.rend(); ++it) {
      const char act = it->first;
      const int idx = it->second;
      const unsigned int mask = (1U << idx);
      const bool occ = ((state & mask) != 0U);

      if (act == '+' && occ) { valid = false; break; }
      if (act == '-' && !occ) { valid = false; break; }

      const unsigned int count = (state & (mask - 1U));
      int ones = 0;
      for (unsigned int c = count; c != 0U; c >>= 1U) {
        if ((c & 1U) != 0U) ++ones;
      }
      if ((ones % 2) != 0) {
        phase = -phase;
      }
      state ^= mask;
    }

    if (valid) {
      mat[state][ket] += coeff * phase;
    }
  }
}

} // namespace

GeneratedHamiltonianResult buildMolecularHamiltonian(
    const MolecularGeometry &geometry,
    MappingType mapping,
    const ActiveSpaceConfig &activeSpace
) {
  const auto integrals = computeMolecularIntegrals(geometry);
  const auto rhf = computeRhfScf(geometry, integrals);
  if (!rhf.converged) {
    throw std::runtime_error("RHF SCF calculation failed to converge");
  }

  const int nao = integrals.nao;
  const int nCore = activeSpace.nCore;
  const int nActive = activeSpace.nActive;

  if (nCore + nActive > nao) {
    throw std::invalid_argument("ActiveSpaceConfig: nCore + nActive exceeds number of molecular orbitals");
  }

  // 1. Core energy shift: E_core = E_nuc + 2 sum_c h_cc + sum_cd [ 2(cc|dd) - (cd|dc) ]
  double eCore = integrals.nuclearRepulsion;
  for (int c = 0; c < nCore; ++c) {
    eCore += 2.0 * rhf.h1Mo[static_cast<std::size_t>(c * nao + c)];
  }
  for (int c = 0; c < nCore; ++c) {
    for (int d = 0; d < nCore; ++d) {
      const double eri_ccdd = rhf.eriMo[static_cast<std::size_t>(((c * nao + c) * nao + d) * nao + d)];
      const double eri_cddc = rhf.eriMo[static_cast<std::size_t>(((c * nao + d) * nao + d) * nao + c)];
      eCore += (2.0 * eri_ccdd - eri_cddc);
    }
  }

  // 2. Effective 1-electron integrals for active space: h_eff_pq = h_mo_pq + sum_c [ 2(pq|cc) - (pc|cq) ]
  std::vector<double> hEff(static_cast<std::size_t>(nActive * nActive), 0.0);
  for (int p = 0; p < nActive; ++p) {
    const int moP = nCore + p;
    for (int q = 0; q < nActive; ++q) {
      const int moQ = nCore + q;
      double val = rhf.h1Mo[static_cast<std::size_t>(moP * nao + moQ)];
      for (int c = 0; c < nCore; ++c) {
        const double eri_pqcc = rhf.eriMo[static_cast<std::size_t>(((moP * nao + moQ) * nao + c) * nao + c)];
        const double eri_pccq = rhf.eriMo[static_cast<std::size_t>(((moP * nao + c) * nao + c) * nao + moQ)];
        val += (2.0 * eri_pqcc - eri_pccq);
      }
      hEff[static_cast<std::size_t>(p * nActive + q)] = val;
    }
  }

  // 3. Active 2-electron integrals: (pq|rs)_act
  const std::size_t nAct4 = static_cast<std::size_t>(nActive * nActive * nActive * nActive);
  std::vector<double> eriAct(nAct4, 0.0);
  for (int p = 0; p < nActive; ++p) {
    const int moP = nCore + p;
    for (int q = 0; q < nActive; ++q) {
      const int moQ = nCore + q;
      for (int r = 0; r < nActive; ++r) {
        const int moR = nCore + r;
        for (int s = 0; s < nActive; ++s) {
          const int moS = nCore + s;
          eriAct[static_cast<std::size_t>(((p * nActive + q) * nActive + r) * nActive + s)] =
              rhf.eriMo[static_cast<std::size_t>(((moP * nao + moQ) * nao + moR) * nao + moS)];
        }
      }
    }
  }

  // 4. Build 16x16 determinant matrix for active space (2 spatial orbitals x 2 spins = 4 spin-orbitals)
  std::array<std::array<std::complex<double>, 16>, 16> hDet{};

  for (int spin = 0; spin < 2; ++spin) {
    for (int p = 0; p < nActive; ++p) {
      for (int q = 0; q < nActive; ++q) {
        const double hval = hEff[static_cast<std::size_t>(p * nActive + q)];
        if (std::abs(hval) > 1e-12) {
          applyLadderTerm(hDet, {{'+', p + spin * 2}, {'-', q + spin * 2}}, hval);
        }
      }
    }
  }

  for (int s1 = 0; s1 < 2; ++s1) {
    for (int s2 = 0; s2 < 2; ++s2) {
      for (int p = 0; p < nActive; ++p) {
        for (int q = 0; q < nActive; ++q) {
          for (int r = 0; r < nActive; ++r) {
            for (int s = 0; s < nActive; ++s) {
              const double erival = 0.5 * eriAct[static_cast<std::size_t>(((p * nActive + q) * nActive + r) * nActive + s)];
              if (std::abs(erival) > 1e-12) {
                applyLadderTerm(hDet, {{'+', p + s1 * 2}, {'+', r + s2 * 2}, {'-', s + s2 * 2}, {'-', q + s1 * 2}}, erival);
              }
            }
          }
        }
      }
    }
  }

  // Calculate physical singlet FCI ground state by projecting onto 3 singlet CSFs
  const double invSqrt2 = 1.0 / std::sqrt(2.0);
  std::array<std::array<double, 16>, 3> csf{};
  csf[0][5] = 1.0;
  csf[1][6] = invSqrt2;
  csf[1][9] = invSqrt2;
  csf[2][10] = 1.0;

  std::array<std::array<double, 3>, 3> hCsf{};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      std::complex<double> sum(0.0, 0.0);
      for (int b1 = 0; b1 < 16; ++b1) {
        for (int b2 = 0; b2 < 16; ++b2) {
          sum += csf[static_cast<std::size_t>(i)][static_cast<std::size_t>(b1)] *
                 hDet[static_cast<std::size_t>(b1)][static_cast<std::size_t>(b2)] *
                 csf[static_cast<std::size_t>(j)][static_cast<std::size_t>(b2)];
        }
      }
      hCsf[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = sum.real();
    }
  }

  const double exactFciElec = lowestSymmetricEigenvalue<3>(hCsf);
  const double exactFciTotal = eCore + exactFciElec;

  GeneratedHamiltonianResult res;
  res.nuclearRepulsion = integrals.nuclearRepulsion;
  res.coreEnergy = eCore;
  res.rhfTotalEnergy = rhf.totalEnergy;
  res.fciExactElectronicEnergy = exactFciElec;
  res.fciExactTotalEnergy = exactFciTotal;
  res.activeElectrons = 2;
  res.activeSpatialOrbitals = nActive;

  if (mapping == MappingType::Parity) {
    std::array<std::array<std::complex<double>, 16>, 16> hParity{};
    auto parityOf = [](unsigned int occ) {
      unsigned int p = 0;
      unsigned int run = 0;
      for (int m = 0; m < 4; ++m) {
        run ^= (occ >> m) & 1U;
        p |= (run << m);
      }
      return p;
    };

    for (unsigned int b = 0; b < 16; ++b) {
      for (unsigned int k = 0; k < 16; ++k) {
        hParity[parityOf(b)][parityOf(k)] = hDet[b][k];
      }
    }

    auto fullIdx = [](unsigned int r) {
      unsigned int p0 = r & 1U;
      unsigned int p2 = (r >> 1U) & 1U;
      return p0 | (1U << 1) | (p2 << 2);
    };

    std::array<std::array<std::complex<double>, 4>, 4> hRed{};
    for (unsigned int i = 0; i < 4; ++i) {
      for (unsigned int j = 0; j < 4; ++j) {
        hRed[i][j] = hParity[fullIdx(i)][fullIdx(j)];
      }
    }

    std::vector<PauliTerm> pauliTerms;
    const std::array<Pauli, 4> paulis{Pauli::I, Pauli::X, Pauli::Y, Pauli::Z};
    const auto pauliMatrix = [](Pauli p, unsigned int row, unsigned int col) -> std::complex<double> {
      if (p == Pauli::I) return (row == col) ? 1.0 : 0.0;
      if (p == Pauli::X) return (row != col) ? 1.0 : 0.0;
      if (p == Pauli::Y) return (row == 0 && col == 1) ? std::complex<double>(0.0, -1.0) : ((row == 1 && col == 0) ? std::complex<double>(0.0, 1.0) : 0.0);
      if (p == Pauli::Z) return (row == 0 && col == 0) ? 1.0 : ((row == 1 && col == 1) ? -1.0 : 0.0);
      return 0.0;
    };

    for (Pauli p1 : paulis) {
      for (Pauli p0 : paulis) {
        std::complex<double> tr(0.0, 0.0);
        for (unsigned int row = 0; row < 4; ++row) {
          const unsigned int r0 = row & 1U;
          const unsigned int r1 = (row >> 1U) & 1U;
          for (unsigned int col = 0; col < 4; ++col) {
            const unsigned int c0 = col & 1U;
            const unsigned int c1 = (col >> 1U) & 1U;
            tr += hRed[row][col] * pauliMatrix(p1, col != 0 ? c1 : c1, r1) * pauliMatrix(p0, c0, r0);
          }
        }
        const double coeff = tr.real() / 4.0;
        if (std::abs(coeff) > 1e-10) {
          pauliTerms.push_back({coeff, {p0, p1}});
        }
      }
    }

    res.qubitCount = 2;
    res.hamiltonian = QubitHamiltonian(2, std::move(pauliTerms));
    res.provenance = (nCore > 0)
        ? "Generated C++ STO-3G Integral Engine // Active-Space Parity Mapping (1,1) // 2 Qubits"
        : "Generated C++ STO-3G Integral Engine // Parity Mapping (1,1) // 2 Qubits";
  } else {
    std::vector<PauliTerm> pauliTerms;
    const std::array<Pauli, 4> paulis{Pauli::I, Pauli::X, Pauli::Y, Pauli::Z};
    const auto pauliMatrix = [](Pauli p, unsigned int row, unsigned int col) -> std::complex<double> {
      if (p == Pauli::I) return (row == col) ? 1.0 : 0.0;
      if (p == Pauli::X) return (row != col) ? 1.0 : 0.0;
      if (p == Pauli::Y) return (row == 0 && col == 1) ? std::complex<double>(0.0, -1.0) : ((row == 1 && col == 0) ? std::complex<double>(0.0, 1.0) : 0.0);
      if (p == Pauli::Z) return (row == 0 && col == 0) ? 1.0 : ((row == 1 && col == 1) ? -1.0 : 0.0);
      return 0.0;
    };

    for (Pauli p3 : paulis) {
      for (Pauli p2 : paulis) {
        for (Pauli p1 : paulis) {
          for (Pauli p0 : paulis) {
            std::complex<double> tr(0.0, 0.0);
            for (unsigned int row = 0; row < 16; ++row) {
              const unsigned int r0 = row & 1U;
              const unsigned int r1 = (row >> 1U) & 1U;
              const unsigned int r2 = (row >> 2U) & 1U;
              const unsigned int r3 = (row >> 3U) & 1U;
              for (unsigned int col = 0; col < 16; ++col) {
                const unsigned int c0 = col & 1U;
                const unsigned int c1 = (col >> 1U) & 1U;
                const unsigned int c2 = (col >> 2U) & 1U;
                const unsigned int c3 = (col >> 3U) & 1U;
                tr += hDet[row][col] *
                      pauliMatrix(p3, c3, r3) *
                      pauliMatrix(p2, c2, r2) *
                      pauliMatrix(p1, c1, r1) *
                      pauliMatrix(p0, c0, r0);
              }
            }
            const double coeff = tr.real() / 16.0;
            if (std::abs(coeff) > 1e-10) {
              pauliTerms.push_back({coeff, {p0, p1, p2, p3}});
            }
          }
        }
      }
    }

    res.qubitCount = 4;
    res.hamiltonian = QubitHamiltonian(4, std::move(pauliTerms));
    res.provenance = (nCore > 0)
        ? "Generated C++ STO-3G Integral Engine // Active-Space Jordan-Wigner Mapping // 4 Qubits"
        : "Generated C++ STO-3G Integral Engine // Jordan-Wigner Mapping // 4 Qubits";
  }

  return res;
}
