#include "test_support.h"
#include <iostream>

namespace test {
void testPauliExpectations();
void testHamiltonianExpectations();
void testSymmetricEigensolverValidation();
void testBenchmarkCatalog();
void testHeHEnergy();
void testH2Energy();
void testGeneratedH2();

void testLexer();
void testParser();
void testSemantics();
void testInstructionOptimization();
void testCodeGeneration();
void testSimulator();
void testCapacity();
void testCompilerFacade();
void testSimulationSession();
void testWorkbenchIntegration();
void testExperimentBundle();
void testMolecularIntegralEngine();
void testCircuitDag();
} // namespace test
int main() {
  try {
    test::testPauliExpectations();
    test::testHamiltonianExpectations();
    test::testSymmetricEigensolverValidation();
    test::testBenchmarkCatalog();
    test::testHeHEnergy();
    test::testH2Energy();
    test::testGeneratedH2();
    test::testLexer();
    test::testParser();
    test::testSemantics();
    test::testInstructionOptimization();
    test::testCodeGeneration();
    test::testSimulator();
    test::testCapacity();
    test::testCompilerFacade();
    test::testSimulationSession();
    test::testWorkbenchIntegration();
    test::testExperimentBundle();
    test::testMolecularIntegralEngine();
    test::testCircuitDag();
  } catch (const std::exception &error) {
    std::cerr << "FAIL after " << test::assertions << " assertions: " << error.what() << '\n';
    return 1;
  }
  std::cout << "Passed " << test::assertions << " assertions.\n";
  return 0;
}
