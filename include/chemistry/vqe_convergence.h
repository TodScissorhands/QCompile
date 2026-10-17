#pragma once

#include <functional>

struct VqeConvergencePoint {
  int evaluation = 0;
  double parameter0 = 0.0;
  double parameter1 = 0.0;
  double electronicEnergy = 0.0;
  double totalEnergy = 0.0;
  bool isBest = false;
};

using VqeCallback = std::function<void(const VqeConvergencePoint &)>;
