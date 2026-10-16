#include "chemistry/heh_energy.h"

#include <iomanip>
#include <iostream>

int main() {
  const HeHResult result = optimizeHeH();
  std::cout << std::fixed << std::setprecision(15)
            << "HeH+ STO-3G Jordan-Wigner singlet benchmark\n"
            << "Parameters: theta=" << result.theta << " phi=" << result.phi << "\n"
            << "Electronic energy: " << result.electronicEnergy << " Ha\n"
            << "Exact singlet electronic energy: "
            << result.exactElectronicEnergy << " Ha\n"
            << "Nuclear repulsion: " << result.nuclearRepulsionEnergy << " Ha\n"
            << "Total Born-Oppenheimer energy: " << result.totalEnergy << " Ha\n"
            << "Energy evaluations: " << result.energyEvaluations << "\n";
  return 0;
}
