#pragma once

#include "compiler/code_generator.h"
#include "chemistry/hamiltonian.h"

#include <vector>

struct HeHResult {
  double theta;
  double phi;
  double electronicEnergy;
  double exactElectronicEnergy;
  double nuclearRepulsionEnergy;
  double totalEnergy;
  int energyEvaluations;
};

const QubitHamiltonian &hehHamiltonian();
std::vector<Instruction> hehSingletTrialCircuit(double theta, double phi);
double evaluateHeHEnergy(double theta, double phi);
double exactHeHSingletEnergy();
HeHResult optimizeHeH();

constexpr double HEH_NUCLEAR_REPULSION_HARTREE = 1.3667896093920315;
constexpr double HEH_ELECTRONIC_REFERENCE_HARTREE = -4.218264154490385;
constexpr double HEH_TOTAL_REFERENCE_HARTREE = -2.8514745450983536;
