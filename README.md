# QCompile

QCompile is a compiler for a small domain-specific language (DSL) that describes quantum circuits. It is intended to compile those descriptions into instructions for a classical, complex-valued state-vector simulator—not real quantum hardware. Its overall design is similar in spirit to quantum software frameworks such as Qiskit and Cirq.

## Intended pipeline

QCompile is planned as a five-stage compiler:

1. **Lexer — DONE**: Converts source characters into tokens such as keywords, gate names, identifiers, commas, and newlines.
2. **Parser — DONE**: Converts the token list into an abstract syntax tree (AST) containing qubit declarations and gate operations.
3. **Semantic Checker — PLANNED**: Will verify that qubits were declared, gates have the correct number of operands, and measurement rules are followed before code generation.
4. **Code Generator — PLANNED**: Will lower the AST into a flat instruction list that the simulator can execute.
5. **Simulator — PLANNED**: Will execute the instruction list using a complex-valued state vector and gate matrices, including simulated measurement.

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

Gates may optionally have one arithmetic expression in parentheses before the qubit list. Expressions can use number literals, identifiers such as `theta`, the operators `+`, `-`, `*`, and `/`, and parentheses. Multiplication and division have higher precedence than addition and subtraction. Unary minus is also supported.

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

## Current status

The lexer and parser are implemented and runnable. The current program accepts a source filename, tokenizes its contents, builds an AST, and prints the tokens and AST. It includes sample programs and error handling with line numbers. Semantic checking, code generation, and state-vector simulation are not implemented yet.

## Project files

`include/ast.h` defines the AST data structures: statement types, statements, and programs. The parser fills these structures, and the driver prints them.

`include/lexer.h` declares token types, the `Token` struct, and the lexer interface used by the driver.

`src/lexer.cpp` reads raw source text and produces the token list declared in `include/lexer.h`. Its output feeds the parser.

`include/parser.h` declares the parser and its helper functions. It connects the token definitions to the AST definitions.

`src/parser.cpp` consumes tokens from the lexer and builds the AST structures from `include/ast.h`. It checks statement structure but does not perform semantic validation.

`src/main.cpp` reads the source filename, runs the file through the lexer and parser, and prints the tokens and AST.

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
g++ -std=c++17 -Wall -Wextra -pedantic -Iinclude src/main.cpp src/lexer.cpp src/parser.cpp -o qcompile
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
```
