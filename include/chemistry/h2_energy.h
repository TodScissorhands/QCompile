#pragma once

#include "chemistry/hamiltonian.h"
#include "compiler/code_generator.h"
#include "chemistry/h2_geometry.h"

#include <array>


constexpr double H2_ELECTRONIC_REFERENCE_HARTREE = -1.857275030202;
constexpr double H2_NUCLEAR_REPULSION_HARTREE = 0.719968994449;

constexpr double H2_TOTAL_REFERENCE_HARTREE = -1.137306035753;

struct H2Result {
  double parameter;
  double electronicEnergy;
  double exactElectronicEnergy;
  double nuclearRepulsionEnergy;
  double totalEnergy;
  int energyEvaluations;
};

struct GeneratedH2EnergyResult {
  double parameter;
  double electronicEnergy;
  double exactElectronicEnergy;
  double referenceElectronicEnergy;
  double nuclearRepulsionEnergy;
  double totalEnergy;
  double exactTotalEnergy;
  double variationalError;
  int energyEvaluations;
  bool converged;
};

double evaluateGeneratedH2Energy(const GeneratedH2 &model, double parameter);
double exactGeneratedH2GroundEnergy(const GeneratedH2 &model);
std::vector<Instruction> h2TrialInstructions(double parameter);
GeneratedH2EnergyResult optimizeGeneratedH2(const GeneratedH2 &model);

QubitHamiltonian h2Hamiltonian();
double evaluateH2Energy(double parameter);
std::array<std::array<double, 4>, 4> h2HamiltonianMatrix();
double exactH2GroundEnergy();
H2Result optimizeH2();
