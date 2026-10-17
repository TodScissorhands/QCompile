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

struct ElementBasisInfo {
  int atomicNumber;
  std::vector<std::vector<GaussianPrimitive>> shells;
};

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

double boysFunction0(double t) {
  if (t < 1e-7) {
    return 1.0 - t / 3.0 + (t * t) / 10.0;
  }
  return 0.5 * std::sqrt(kPi / t) * std::erf(std::sqrt(t));
}

double primitiveOverlap(double a, const std::array<double, 3> &A,
                        double b, const std::array<double, 3> &B) {
  const double gamma = a + b;
  const double dx = A[0] - B[0];
  const double dy = A[1] - B[1];
  const double dz = A[2] - B[2];
  const double r2 = dx * dx + dy * dy + dz * dz;
  return std::pow(2.0 * std::sqrt(a * b) / gamma, 1.5) * std::exp(-a * b / gamma * r2);
}

double primitiveKinetic(double a, const std::array<double, 3> &A,
                        double b, const std::array<double, 3> &B) {
  const double gamma = a + b;
  const double dx = A[0] - B[0];
  const double dy = A[1] - B[1];
  const double dz = A[2] - B[2];
  const double r2 = dx * dx + dy * dy + dz * dz;
  const double s = primitiveOverlap(a, A, b, B);
  return (a * b / gamma) * (3.0 - 2.0 * (a * b / gamma) * r2) * s;
}

double primitiveNuclear(double a, const std::array<double, 3> &A,
                        double b, const std::array<double, 3> &B,
                        const std::array<double, 3> &C, double Zc) {
  const double gamma = a + b;
  const double Px = (a * A[0] + b * B[0]) / gamma;
  const double Py = (a * A[1] + b * B[1]) / gamma;
  const double Pz = (a * A[2] + b * B[2]) / gamma;
  const double dx = Px - C[0];
  const double dy = Py - C[1];
  const double dz = Pz - C[2];
  const double rpc2 = dx * dx + dy * dy + dz * dz;
  const double s = primitiveOverlap(a, A, b, B);
  return -Zc * (2.0 * std::sqrt(gamma / kPi)) * s * boysFunction0(gamma * rpc2);
}

double primitiveEri(double a, const std::array<double, 3> &A,
                    double b, const std::array<double, 3> &B,
                    double c, const std::array<double, 3> &C,
                    double d, const std::array<double, 3> &D) {
  const double g1 = a + b;
  const double g2 = c + d;
  const double Px = (a * A[0] + b * B[0]) / g1;
  const double Py = (a * A[1] + b * B[1]) / g1;
  const double Pz = (a * A[2] + b * B[2]) / g1;
  const double Qx = (c * C[0] + d * D[0]) / g2;
  const double Qy = (c * C[1] + d * D[1]) / g2;
  const double Qz = (c * C[2] + d * D[2]) / g2;
  const double dx = Px - Qx;
  const double dy = Py - Qy;
  const double dz = Pz - Qz;
  const double rpq2 = dx * dx + dy * dy + dz * dz;
  const double rho = g1 * g2 / (g1 + g2);
  const double s1 = primitiveOverlap(a, A, b, B);
  const double s2 = primitiveOverlap(c, C, d, D);
  return s1 * s2 * (2.0 * std::sqrt(rho / kPi)) * boysFunction0(rho * rpq2);
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
    if (maxVal <= offDiagTol) {
      converged = true;
      break;
    }

    const double app = a[static_cast<std::size_t>(p * n + p)];
    const double aqq = a[static_cast<std::size_t>(q * n + q)];
    const double apq = a[static_cast<std::size_t>(p * n + q)];
    const double theta = 0.5 * std::atan2(2.0 * apq, aqq - app);
    const double c = std::cos(theta);
    const double s = std::sin(theta);

    for (int k = 0; k < n; ++k) {
      if (k != p && k != q) {
        const double akp = a[static_cast<std::size_t>(k * n + p)];
        const double akq = a[static_cast<std::size_t>(k * n + q)];
        a[static_cast<std::size_t>(k * n + p)] = a[static_cast<std::size_t>(p * n + k)] = c * akp - s * akq;
        a[static_cast<std::size_t>(k * n + q)] = a[static_cast<std::size_t>(q * n + k)] = s * akp + c * akq;
      }
      const double vkp = eigenvectors[static_cast<std::size_t>(k * n + p)];
      const double vkq = eigenvectors[static_cast<std::size_t>(k * n + q)];
      eigenvectors[static_cast<std::size_t>(k * n + p)] = c * vkp - s * vkq;
      eigenvectors[static_cast<std::size_t>(k * n + q)] = s * vkp + c * vkq;
    }
    a[static_cast<std::size_t>(p * n + p)] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
    a[static_cast<std::size_t>(q * n + q)] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
    a[static_cast<std::size_t>(p * n + q)] = a[static_cast<std::size_t>(q * n + p)] = 0.0;
  }

  if (!converged) {
    throw std::runtime_error("Dynamic Jacobi eigensolver did not converge");
  }

  for (int i = 0; i < n; ++i) {
    eigenvalues[static_cast<std::size_t>(i)] = a[static_cast<std::size_t>(i * n + i)];
  }

  // Sort eigenvalues ascending and permute eigenvector columns
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
        std::swap(eigenvectors[static_cast<std::size_t>(k * n + i)],
                  eigenvectors[static_cast<std::size_t>(k * n + minIdx)]);
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
    if (atom.element != "H" && atom.element != "He") {
      throw std::invalid_argument("Unsupported element: " + atom.element + " (STO-3G integral engine supports H, He)");
    }
    for (int d = 0; d < 3; ++d) {
      if (!std::isfinite(atom.coordinates[static_cast<std::size_t>(d)])) {
        throw std::invalid_argument("Atomic coordinates must be finite numbers");
      }
    }
    totalZ += (atom.element == "He" ? 2 : 1);

    // Coincident nuclei check
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
    atom.atomicNumber = (atom.element == "He" ? 2 : 1);
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
    ContractedGaussian cg;
    cg.atomIndex = static_cast<int>(i);
    cg.angularMomentum = 0;
    cg.lx = cg.ly = cg.lz = 0;
    cg.centerBohr = atom.coordinates;
    cg.primitives = (atom.element == "He") ? kHe_1s : kH_1s;
    basis.basisFunctions.push_back(std::move(cg));
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
    for (int j = i; j < nao; ++j) {
      const auto &bfJ = basis.basisFunctions[static_cast<std::size_t>(j)];

      double sVal = 0.0;
      double tVal = 0.0;
      double vVal = 0.0;

      for (const auto &pI : bfI.primitives) {
        for (const auto &pJ : bfJ.primitives) {
          const double cProd = pI.contractionCoefficient * pJ.contractionCoefficient;
          sVal += cProd * primitiveOverlap(pI.exponent, bfI.centerBohr, pJ.exponent, bfJ.centerBohr);
          tVal += cProd * primitiveKinetic(pI.exponent, bfI.centerBohr, pJ.exponent, bfJ.centerBohr);

          for (const auto &nuc : bohrAtoms) {
            vVal += cProd * primitiveNuclear(pI.exponent, bfI.centerBohr, pJ.exponent, bfJ.centerBohr,
                                             nuc.coordinates, static_cast<double>(nuc.atomicNumber));
          }
        }
      }

      res.overlap[static_cast<std::size_t>(i * nao + j)] = res.overlap[static_cast<std::size_t>(j * nao + i)] = sVal;
      res.kinetic[static_cast<std::size_t>(i * nao + j)] = res.kinetic[static_cast<std::size_t>(j * nao + i)] = tVal;
      res.nuclear[static_cast<std::size_t>(i * nao + j)] = res.nuclear[static_cast<std::size_t>(j * nao + i)] = vVal;
      res.hCore[static_cast<std::size_t>(i * nao + j)] = res.hCore[static_cast<std::size_t>(j * nao + i)] = (tVal + vVal);
    }
  }

  // 2. Two-electron integrals: (ij|kl) in chemist notation with 8-fold symmetry
  for (int i = 0; i < nao; ++i) {
    const auto &bfI = basis.basisFunctions[static_cast<std::size_t>(i)];
    for (int j = 0; j <= i; ++j) {
      const auto &bfJ = basis.basisFunctions[static_cast<std::size_t>(j)];
      const int ij = i * (i + 1) / 2 + j;

      for (int k = 0; k < nao; ++k) {
        const auto &bfK = basis.basisFunctions[static_cast<std::size_t>(k)];
        for (int l = 0; l <= k; ++l) {
          const auto &bfL = basis.basisFunctions[static_cast<std::size_t>(l)];
          const int kl = k * (k + 1) / 2 + l;

          if (ij < kl) continue; // 8-fold symmetry canonical loop

          double eriVal = 0.0;
          for (const auto &pI : bfI.primitives) {
            for (const auto &pJ : bfJ.primitives) {
              const double cIJ = pI.contractionCoefficient * pJ.contractionCoefficient;
              for (const auto &pK : bfK.primitives) {
                for (const auto &pL : bfL.primitives) {
                  const double cTotal = cIJ * pK.contractionCoefficient * pL.contractionCoefficient;
                  eriVal += cTotal * primitiveEri(pI.exponent, bfI.centerBohr, pJ.exponent, bfJ.centerBohr,
                                                  pK.exponent, bfK.centerBohr, pL.exponent, bfL.centerBohr);
                }
              }
            }
          }

          auto setEri = [&](int u, int v, int w, int x) {
            res.eri[static_cast<std::size_t>(((u * nao + v) * nao + w) * nao + x)] = eriVal;
          };

          // Assign all 8 permutations
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
  (void)geometry;
  const int nao = integrals.nao;
  if (nao != 2) {
    throw std::invalid_argument("Current RHF SCF implementation is specialized for 2 spatial orbitals (nao=2)");
  }
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

  // SCF loop: closed-shell 2 electrons (1 occupied spatial orbital)
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

    // Update density matrix P_uv = 2 * C_u0 * C_v0 (1 closed-shell orbital)
    std::vector<double> Pnew(static_cast<std::size_t>(nao * nao), 0.0);
    double maxDelta = 0.0;
    for (int u = 0; u < nao; ++u) {
      for (int v = 0; v < nao; ++v) {
        Pnew[static_cast<std::size_t>(u * nao + v)] = 2.0 * Cnew[static_cast<std::size_t>(u * nao + 0)] * Cnew[static_cast<std::size_t>(v * nao + 0)];
        maxDelta = std::max(maxDelta, std::abs(Pnew[static_cast<std::size_t>(u * nao + v)] - P[static_cast<std::size_t>(u * nao + v)]));
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
  res.nElectrons = 2;
  res.nOccupied = 1;
  res.nOrbitals = nao;
  res.moCoefficients = C;
  res.moEnergies = moE;

  // Transform AO integrals to MO integrals:
  // h1Mo_pq = sum_uv C_up Hcore_uv C_vq
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

  // 4-index transformation to MO ERI: (pq|rs)
  res.eriMo.assign(nmo4, 0.0);
  for (int p = 0; p < nao; ++p) {
    for (int q = 0; q < nao; ++q) {
      for (int r = 0; r < nao; ++r) {
        for (int s = 0; s < nao; ++s) {
          double sum = 0.0;
          for (int u = 0; u < nao; ++u) {
            for (int v = 0; v < nao; ++v) {
              for (int w = 0; w < nao; ++w) {
                for (int x = 0; x < nao; ++x) {
                  sum += C[static_cast<std::size_t>(u * nao + p)] *
                         C[static_cast<std::size_t>(v * nao + q)] *
                         C[static_cast<std::size_t>(w * nao + r)] *
                         C[static_cast<std::size_t>(x * nao + s)] *
                         integrals.getEri(u, v, w, x);
                }
              }
            }
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

      // Jordan-Wigner phase (-1)^count of occupied modes before idx
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
    MappingType mapping
) {
  const auto integrals = computeMolecularIntegrals(geometry);
  const auto rhf = computeRhfScf(geometry, integrals);
  if (!rhf.converged) {
    throw std::runtime_error("RHF SCF calculation failed to converge");
  }

  // Build 16x16 determinant matrix for 4 modes (2 spatial orbitals x 2 spins)
  std::array<std::array<std::complex<double>, 16>, 16> hDet{};

  // 1-electron terms: h_pq a_{p,s}^dag a_{q,s}
  for (int spin = 0; spin < 2; ++spin) {
    for (int p = 0; p < 2; ++p) {
      for (int q = 0; q < 2; ++q) {
        const double hval = rhf.h1Mo[static_cast<std::size_t>(p * 2 + q)];
        if (std::abs(hval) > 1e-12) {
          applyLadderTerm(hDet, {{'+', p + spin * 2}, {'-', q + spin * 2}}, hval);
        }
      }
    }
  }

  // 2-electron terms: 0.5 * (pq|rs) a_{p,s1}^dag a_{r,s2}^dag a_{s,s2} a_{q,s1}
  for (int s1 = 0; s1 < 2; ++s1) {
    for (int s2 = 0; s2 < 2; ++s2) {
      for (int p = 0; p < 2; ++p) {
        for (int q = 0; q < 2; ++q) {
          for (int r = 0; r < 2; ++r) {
            for (int s = 0; s < 2; ++s) {
              const double erival = 0.5 * rhf.eriMo[static_cast<std::size_t>(((p * 2 + q) * 2 + r) * 2 + s)];
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
  // CSF 0: |0101> (index 5)
  // CSF 1: (|0110> + |1001>)/sqrt(2) (indices 6 and 9)
  // CSF 2: |1010> (index 10)
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
  const double exactFciTotal = exactFciElec + integrals.nuclearRepulsion;

  GeneratedHamiltonianResult res;
  res.nuclearRepulsion = integrals.nuclearRepulsion;
  res.rhfTotalEnergy = rhf.totalEnergy;
  res.fciExactElectronicEnergy = exactFciElec;
  res.fciExactTotalEnergy = exactFciTotal;

  if (mapping == MappingType::Parity) {
    // Parity transformation: p_i = sum_{k=0}^i occ_k mod 2
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

    // Two-electron reduction: p1 = 1, p3 = 0.
    // Reduced qubits: q0 = p0, q1 = p2 (indices 2, 3, 6, 7)
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

    // Extract 2-qubit Pauli terms
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
    res.provenance = "Generated C++ STO-3G Integral Engine // Parity Mapping (1,1) // 2 Qubits";
  } else {
    // Jordan-Wigner 4-qubit mapping
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
    res.provenance = "Generated C++ STO-3G Integral Engine // Jordan-Wigner Mapping // 4 Qubits";
  }

  return res;
}
