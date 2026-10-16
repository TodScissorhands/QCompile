#pragma once

#include "chemistry/hamiltonian.h"

#include <string>

struct GeneratedH2 {
  double distanceAngstrom;
  double distanceBohr;
  double nuclearRepulsionEnergy;
  double referenceElectronicEnergy;
  QubitHamiltonian electronicHamiltonian;

  GeneratedH2(double distanceAngstromValue, double distanceBohrValue,
              double nuclearRepulsionEnergyValue,
              double referenceElectronicEnergyValue,
              QubitHamiltonian hamiltonian);
};

// Reads only the versioned two-qubit neutral singlet H2 format emitted by
// benchmarks/h2_geometry/generate_h2.py. Unsupported model metadata is rejected.
GeneratedH2 loadGeneratedH2(const std::string &path);
