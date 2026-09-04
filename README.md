# QCompile

QCompile is a compiler for a small domain-specific language (DSL) that describes quantum circuits. It is intended to compile those descriptions into instructions for a classical, complex-valued state-vector simulator—not real quantum hardware. Its overall design is similar in spirit to quantum software frameworks such as Qiskit and Cirq.

## Intended pipeline

QCompile has a five-stage compiler pipeline:

1. **Lexer — DONE**: Converts source characters into tokens such as keywords, gate names, identifiers, commas, and newlines.
2. **Parser — DONE**: Converts the token list into an abstract syntax tree (AST) containing qubit declarations and gate operations.
3. **Semantic Checker — DONE**: Validates declared qubits, gate arity, CNOT control/target pairs, and numeric gate parameters before code generation.
4. **Code Generator — DONE**: Lowers the validated AST into a flat instruction list with zero-based qubit indices.
5. **Simulator — DONE**: Executes instructions using a complex-valued state vector, gate matrices, CNOT swaps, and seeded measurement sampling.

## Language specification

The language has one statement per line:

```text
qubit q0, q1
H q0
X q0
Z q0
CNOT q0, q1
measure q0, q1
```

The keywords are `qubit` and `measure`. The supported gate names are `H`, `X`, `Z`, `CNOT`, and `RZ`. Identifiers such as `q0` and `q1` name qubits. Commas separate multiple names, and newlines separate statements.

Variables are declared with `let`, followed by a name, `=`, and an arithmetic expression. A variable must be declared before it is used in a later variable declaration or a gate parameter.

```text
let theta = 0.5
let phi = 2 * 3.14159 / 4
```

Gates may optionally have one arithmetic expression in parentheses before the qubit list. Expressions can use number literals, identifiers such as `theta`, the operators `+`, `-`, `*`, and `/`, and parentheses. Multiplication and division have higher precedence than addition and subtraction. Unary minus is also supported. The current simulator requires parameters to evaluate to numbers, so an identifier without an assigned value is a semantic error.

```text
qubit q0
RZ(0.5 + 0.25) q0
RZ(2 * theta - 0.1) q0
```

Example Bell-state program:

```text
qubit q0, q1
H q0
CNOT q0, q1
measure q0, q1
```

## Applications and motivation

Compiling a high-level circuit description into executable simulation instructions mirrors the workflow used by real quantum SDKs. This project demonstrates core compiler concepts—lexing, parsing, static validation, and code generation—in an emerging technical domain, while using a classical simulator that can be run and inspected on an ordinary computer.

## Comparison to Related Work

QCompile is a small-scale but functionally complete compiler for a custom-designed quantum-circuit DSL. It implements all five classic stages used in this project—lexing, parsing, semantic checking, code generation, and simulation. It is not a novel research contribution; it is a genuine working implementation of established compiler techniques applied to quantum circuits.

Its workflow is similar in spirit to IBM's Qiskit and Google's Cirq: a circuit is described in code, compiled into an executable representation, and run on a simulator or, in larger platforms, real hardware. QCompile's DSL is also a simplified, custom cousin of OpenQASM, an industry-standard quantum assembly language. Both share core ideas such as qubits, gates, and measurement, while QCompile intentionally uses a smaller self-designed syntax.

Several scope decisions keep the project focused and fully demonstrable. QCompile does not perform circuit optimization, such as reducing gate counts or adapting a circuit to hardware connectivity, because those features are outside this project's size. Its simulator uses a straightforward state vector of size 2^n: the textbook-correct method and an appropriate choice for the small qubit counts in the examples, rather than a large-scale simulation technique. The DSL also omits control flow, registers, and subroutines so its grammar can be fully understood and demonstrated within the project timeline instead of partially implementing a broader language.

Compared with a minimal version of this assignment, QCompile includes real recursive-descent expression parsing with operator precedence through `parseExpression`, `parseTerm`, and `parseFactor`, rather than only a flat command list. Its semantic checker evaluates arithmetic expressions and reports undefined variables, not just structural errors. The result is a complete connected five-stage pipeline rather than a partial implementation.

## Current status

All five stages are implemented and runnable. The program accepts a source filename, tokenizes it, builds an AST, validates semantics, prints generated instructions, and simulates the resulting circuit. It includes sample programs and line-numbered parser and semantic errors.

## Project files

`include/ast.h` defines the AST data structures: statement types, statements, and programs. The parser fills these structures, and the driver prints them.

`include/lexer.h` declares token types, the `Token` struct, and the lexer interface used by the driver.

`src/lexer.cpp` reads raw source text and produces the token list declared in `include/lexer.h`. Its output feeds the parser.

`include/parser.h` declares the parser and its helper functions. It connects the token definitions to the AST definitions.

`src/parser.cpp` consumes tokens from the lexer and builds the AST structures from `include/ast.h`. It checks statement structure but does not perform semantic validation.

`include/semantic_checker.h` and `src/semantic_checker.cpp` validate the parsed AST, collect all semantic errors, evaluate variable declarations and gate parameter expressions, and record declared-qubit order.

`include/code_generator.h` and `src/code_generator.cpp` convert a valid AST into simple simulator instructions that use integer qubit indices.

`include/simulator.h` and `src/simulator.cpp` maintain the complex-valued state vector, apply supported gates, and sample measurement outcomes.

`src/main.cpp` reads the source filename and connects all five stages: lexer, parser, semantic checker, code generator, and simulator.

`samples/` contains example input programs. Its `README.md` describes what each sample demonstrates.

## Build

From the QCompile root directory, run:

```bash
make
```

The Makefile compiles the `.cpp` files in `src/`, looks for headers in `include/`, and creates `qcompile` in the project root. To remove the generated executable, run:

```bash
make clean
```

If `make` is not available, use the equivalent compiler command:

```bash
g++ -std=c++17 -Wall -Wextra -pedantic -Iinclude src/main.cpp src/lexer.cpp src/parser.cpp src/semantic_checker.cpp src/code_generator.cpp src/simulator.cpp -o qcompile
```

## Run

Pass a DSL source file as the command-line argument:

```bash
./qcompile samples/bell_state.txt
```

Other samples can be run in the same way:

```bash
./qcompile samples/single_gate.txt
./qcompile samples/sequential_gates.txt
./qcompile samples/three_qubit.txt
./qcompile samples/invalid_gate.txt
./qcompile samples/parameterized_gate.txt
./qcompile samples/precedence_check.txt
./qcompile samples/expression_with_parens.txt
./qcompile samples/invalid_expression.txt
./qcompile samples/undefined_variable.txt
```
