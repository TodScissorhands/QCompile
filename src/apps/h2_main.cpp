#include "chemistry/h2_energy.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>

int main() {
  try {
    const H2Result result = optimizeH2();
    std::cout << std::fixed << std::setprecision(12);
    std::cout << "H2: neutral singlet, 0.735 Angstrom, STO-3G, 2 electrons\n";
    std::cout << "Encoding: 2-qubit parity mapping, particle-number reduced\n";
    std::cout << "Optimized parameter (radians): " << result.parameter << '\n';
    std::cout << "Electronic energy (Hartree): "
              << result.electronicEnergy << '\n';
    std::cout << "Exact mapped-Hamiltonian minimum (Hartree): "
              << result.exactElectronicEnergy << '\n';
    std::cout << "Nuclear repulsion (Hartree): "
              << result.nuclearRepulsionEnergy << '\n';
    std::cout << "Total Born-Oppenheimer energy (Hartree): "
              << result.totalEnergy << '\n';
    std::cout << "Energy evaluations: " << result.energyEvaluations << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
