#include "test_support.h"

#include <string>
#include <vector>

namespace test {

void testSemantics() {
  Program program;
  auto result = check("qubit q0\nqubit q0\n", program);
  require(!result.errors.empty(), "duplicate qubit declaration rejected");
  result = check("H missing\n", program);
  require(!result.errors.empty(), "undeclared qubit rejected");
  result = check("RZ(theta) q0\nqubit q0\n", program);
  require(!result.errors.empty(), "qubit and variable use before declaration rejected");
  result = check("qubit q0, q1\nH q0, q1\n", program);
  require(!result.errors.empty(), "invalid gate arity rejected");
  result = check("qubit q0, q1\nCNOT q0\n", program);
  require(!result.errors.empty(), "invalid CNOT arity rejected");
  result = check("qubit q0\nCNOT q0, q0\n", program);
  require(!result.errors.empty(), "identical CNOT operands rejected");

  // Phase 5 C1: Duplicate qubit operands within a single measurement statement rejected
  result = check("qubit q0\nmeasure q0, q0\n", program);
  require(!result.errors.empty(), "duplicate measurement operands rejected");
  require(result.errors[0].find("duplicate qubit 'q0' in measure statement at line 2") != std::string::npos,
          "expected line-located diagnostic for duplicate measurement qubit");

  // Distinct operands in measurement remain valid
  result = check("qubit q0, q1\nmeasure q0, q1\n", program);
  require(result.errors.empty(), "distinct measurement operands accepted");
  result = check("qubit q0\nRZ(1 / 0) q0\n", program);
  require(!result.errors.empty(), "division by zero rejected");
  require(result.errors[0].find("division by zero in gate parameter") != std::string::npos,
          "division by zero context in gate parameter");
  result = check("let theta = 1 / 0\n", program);
  require(!result.errors.empty(), "division by zero in declaration rejected");
  require(result.errors[0].find("division by zero in declaration") != std::string::npos,
          "division by zero context in declaration");

  // Phase 1 numeric safety tests
  // Oversized literal: out of range
  std::string hugeNum(400, '9');
  result = check("qubit q0\nlet theta = " + hugeNum + "\n", program);
  require(!result.errors.empty(), "oversized literal rejected");
  require(result.errors[0].find("out of range") != std::string::npos,
          "oversized literal error message");

  // Arithmetic overflow resulting in non-finite
  std::string largeNum = "1" + std::string(300, '0');
  result = check("qubit q0\nlet theta = " + largeNum + " * " + largeNum + "\n", program);
  require(!result.errors.empty(), "arithmetic overflow rejected");
  require(result.errors[0].find("non-finite arithmetic result") != std::string::npos,
          "arithmetic overflow error message");

  // Non-finite gate parameter
  result = check("qubit q0\nRZ(" + largeNum + " * " + largeNum + ") q0\n", program);
  require(!result.errors.empty(), "non-finite gate parameter rejected");

  // Valid finite expressions
  result = check("qubit q0\nlet x = 3.1415926535\nlet y = (x + 2.0) * 0.5 - 1.0 / 4.0\nRZ(y) q0\n", program);
  require(result.errors.empty(), "valid finite expressions accepted");
  require(std::abs(result.variableValues[1] - ((3.1415926535 + 2.0) * 0.5 - 0.25)) < 1e-12,
          "evaluated variable matches expected value");

  // Unary minus with non-finite
  result = check("qubit q0\nlet x = -(" + largeNum + " * " + largeNum + ")\n", program);
  require(!result.errors.empty(), "unary minus on non-finite rejected");

  result = check("let theta = 2\nlet theta = 3\n", program);
  require(!result.errors.empty(), "duplicate variable declaration rejected");
  result = check("let later = earlier\n", program);
  require(!result.errors.empty(), "variable use before declaration rejected");
}

} // namespace test
