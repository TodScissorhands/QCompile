# QCompile Technical Documentation

## 1. Project Overview

QCompile is a small C++ compiler for a custom domain-specific language (DSL) that describes quantum circuits. It compiles a text program into instructions for a classical state-vector simulator. It does not communicate with real quantum hardware.

The implemented pipeline has six stages:

1. **Lexer** — Reads source characters and creates tokens such as `KeywordQubit`, `Gate`, `Identifier`, `Number`, `Comma`, `Newline`, and arithmetic operators. Tokens preserve their original text and source location.
2. **Parser** — Reads the token list and creates an abstract syntax tree (AST). The AST represents qubit declarations, variable declarations, gate operations, and arithmetic expression trees.
3. **Semantic Checker** — Validates information that cannot be checked by token shape alone: qubits must be declared before use, gates must receive the right number of qubits, CNOT must use different control and target qubits, and parameter expressions must evaluate to numbers using already-declared variables.
4. **Code Generator** — Converts valid AST statements into a flat vector of instructions. It changes qubit names such as `q0` into zero-based integer indices such as `0`.
5. **Instruction optimizer** — Removes adjacent identical H, X, Z, and ordered CNOT pairs. It leaves RZ, measurements, unknown gates, and malformed instructions unchanged.
6. **Simulator** — Maintains a vector of complex amplitudes, applies the generated gates, calculates measurement probabilities, and samples a measurement result with a fixed random seed.

The driver in `src/apps/main.cpp` connects these stages in order. If parsing fails, later stages do not run. If semantic checking reports one or more errors, code generation and simulation do not run.

## 2. Language Specification

### Statements

QCompile uses one statement per line. Blank lines are accepted. The currently implemented statement forms are:

```text
qubit q0, q1
let theta = 0.5
let phi = 2 * 3.14159 / 4
H q0
X q0
Z q0
CNOT q0, q1
RZ(0.5) q0
RZ(2 * theta - 0.1) q0
measure q0, q1
```

`qubit` declares one or more qubit names. `let` declares a numeric variable. `H`, `X`, `Z`, `CNOT`, and `RZ` are recognized gate names. `measure` is parsed as a gate-operation-style statement but generates a measurement instruction.

Parameters are structurally optional after any parsed gate operation, because the parser accepts an optional parenthesized expression after a gate name or `measure`. Semantic checking then enforces the current rule: only `RZ` may have a parameter, and `RZ` must have one. Variables must already have been declared by an earlier `let` statement before they can be used.

### Arithmetic expressions

Expressions support decimal or integer number literals, variables, unary minus, `+`, `-`, `*`, `/`, and parentheses. Multiplication and division bind more tightly than addition and subtraction.

```text
RZ(2 + 3 * 4) q0
RZ((2 + 3) * 4) q0
RZ(-0.5) q0
```

The lexer accepts a number that starts with a digit, followed by digits, with an optional decimal point and more digits. For example, `0.5`, `2`, and `3.14` are valid number tokens. An identifier starts with a letter or underscore and may continue with letters, digits, or underscores.

### Implemented BNF grammar

This grammar describes the structure implemented by `parser.cpp`. Newlines separate statements; the parser accepts leading, repeated, and trailing newlines, and it also accepts an EOF immediately after the final statement.

```bnf
<program>              ::= { <newline> } { <statement> <line-end> }
<line-end>             ::= <newline> { <newline> } | <EOF>

<statement>            ::= <qubit-declaration>
                         | <variable-declaration>
                         | <gate-operation>

<qubit-declaration>    ::= "qubit" <name-list>
<variable-declaration> ::= "let" <identifier> "=" <expression>

<gate-operation>       ::= <operation-name> [ "(" <expression> ")" ] <name-list>
<operation-name>       ::= "H" | "X" | "Z" | "CNOT" | "RZ" | "measure"

<name-list>            ::= <identifier> { "," <identifier> }

<expression>           ::= <term> { ( "+" | "-" ) <term> }
<term>                 ::= <factor> { ( "*" | "/" ) <factor> }
<factor>               ::= <number>
                         | <identifier>
                         | "-" <factor>
                         | "(" <expression> ")"
```

The parser accepts the optional parameter syntax after every `<operation-name>`. The semantic checker narrows that syntax to a meaningful subset by allowing it only for `RZ`.

### Worked Bell-state example and real pipeline trace

Input file `samples/bell_state.txt`:

```text
qubit q0, q1
H q0
CNOT q0, q1
measure q0, q1
```

The following is the exact output captured by running the current executable with `./qcompile samples/bell_state.txt`:

```text
Tokens:
  KeywordQubit("qubit", line 1)
  Identifier("q0", line 1)
  Comma(",", line 1)
  Identifier("q1", line 1)
  Newline("\n", line 1)
  Gate("H", line 2)
  Identifier("q0", line 2)
  Newline("\n", line 2)
  Gate("CNOT", line 3)
  Identifier("q0", line 3)
  Comma(",", line 3)
  Identifier("q1", line 3)
  Newline("\n", line 3)
  KeywordMeasure("measure", line 4)
  Identifier("q0", line 4)
  Comma(",", line 4)
  Identifier("q1", line 4)
  Newline("\n", line 4)
  EndOfFile("", line 5)

AST:
QubitDecl: q0, q1
GateOp: H [q0]
GateOp: CNOT [q0, q1]
GateOp: measure [q0, q1]
Parsed 4 statements from 4 lines.

Instructions:
APPLY H [0]
APPLY CNOT [0, 1]
MEASURE [0, 1]

Probabilities:
  00: 0.500000
  11: 0.500000
Result: 11
```

The simulator uses `std::mt19937` seeded with `12345`, so this build produces the shown sampled result consistently. The significant physical result is the displayed equal probability of `00` and `11`.

## 3. Header File Reference

### `include/compiler/ast.h`

This header defines the structures shared between parsing, semantic checking, code generation, and printing.

#### `enum class StatementType`

`StatementType` labels the meaning of a `Statement` without using an inheritance hierarchy.

- `QubitDeclaration` means the statement stores declared qubit names in `names`.
- `VariableDeclaration` means the statement stores a variable name and expression in `variableName` and `variableExpression`.
- `GateOperation` means the statement stores a gate or measurement name and its qubit operands.

#### `enum class ExpressionType`

`ExpressionType` labels the kind of expression node.

- `NumberLiteral` stores source number text.
- `Identifier` stores a variable name.
- `BinaryOperation` stores an operator and two child expressions.
- `UnaryOperation` currently represents unary minus and stores one child expression.

#### `struct Expression`

- `type` identifies which expression form this node represents.
- `text` stores the literal number text or identifier text for leaf nodes.
- `operation` stores an operator character such as `+`, `*`, or `-` for operation nodes.
- `children` stores nested expressions. A binary node has two children in left-then-right order; a unary-minus node has one child. Using a vector permits a recursive tree without smart pointers.

#### `struct Statement`

One struct is reused for all statement kinds; fields relevant to the selected `type` are filled by the parser.

- `type` identifies the statement category.
- `line` stores the source line used by semantic error messages.
- `names` holds qubit names for `QubitDeclaration`.
- `variableName` holds the left side of a `let` declaration.
- `variableExpression` holds the parsed right side of a `let` declaration.
- `gate` holds a gate name such as `H`, `CNOT`, or the text `measure`.
- `qubits` holds gate or measurement operand names.
- `hasParameter` states whether the operation had parentheses containing an expression.
- `parameter` holds that parsed parameter expression when `hasParameter` is true.

#### `struct Program`

- `statements` holds every parsed statement in source order. Later stages walk this vector from beginning to end, which is important because declarations must come before uses.

### `include/compiler/lexer.h`

#### `enum class TokenType`

`TokenType` lists every token category produced by the lexer: `KeywordQubit`, `KeywordMeasure`, `KeywordLet`, `Gate`, `Identifier`, `Number`, `Comma`, `Newline`, `EndOfFile`, `LeftParen`, `RightParen`, `Plus`, `Minus`, `Star`, `Slash`, and `Equal`.

#### `struct Token`

- `type` is the category from `TokenType`.
- `text` is the original spelling, for example `"CNOT"`, `"0.5"`, or `"+"`.
- `line` stores the one-based input line number.
- `column` stores the one-based input column. The current displayed error messages use the line number, but the lexer still records the column.

#### `class Lexer`

- `Lexer(const std::string &sourceText)` stores the complete source text in the private `source` field.
- `std::vector<Token> tokenize() const` scans `source` and returns the complete token list, ending with `EndOfFile`.
- `source` is the private source-text copy scanned by `tokenize`.

#### `tokenTypeName`

- `const char *tokenTypeName(TokenType type)` converts a token enum value into readable text for the token trace printed by `main.cpp`.

### `include/compiler/parser.h`

#### `class Parser`

- `Parser(const std::vector<Token> &tokenList)` stores a reference to the lexer output and begins at token index zero.
- `Program parse()` parses the complete token list into a `Program`.
- `tokens` is the private reference to the token vector; it is not copied.
- `current` is the private index of the next token to inspect.
- `peek()` returns the current token without advancing.
- `advance()` returns the current token and increments `current`.
- `check(TokenType type)` reports whether the current token has a requested type.
- `match(TokenType type)` consumes the current token only when it has that type.
- `skipNewlines()` consumes one or more newline tokens so blank lines are ignored.
- `expect(TokenType type)` requires a token type and throws a line-numbered parser error when it is absent.
- `parseNameList()` parses one identifier followed by zero or more comma-and-identifier pairs.
- `parseExpression()` parses addition and subtraction expressions.
- `parseTerm()` parses multiplication and division expressions.
- `parseFactor()` parses a number, identifier, unary minus, or parenthesized expression.
- `parseStatement()` chooses and parses one qubit declaration, variable declaration, or gate operation.

### `include/compiler/semantic_checker.h`

#### `struct SemanticResult`

`SemanticResult` returns both errors and the information needed by code generation.

- `errors` collects all semantic error strings instead of stopping after the first one.
- `qubitNames` stores declared qubits in declaration order, which becomes code-generation index order.
- `variableNames` stores valid variable names in declaration order.
- `variableValues` stores the corresponding evaluated `double` values at the same indices as `variableNames`.
- `parameterValues` stores evaluated parameter values aligned with `Program::statements`; non-parameter statements keep `0.0`.

#### `checkSemantics`

- `SemanticResult checkSemantics(const Program &program)` walks the AST in order, collects all semantic errors, evaluates valid variables and parameters, and returns the result.

### `include/compiler/code_generator.h`

#### `enum class InstructionType`

- `ApplyGate` marks an executable gate instruction.
- `Measure` marks a measurement instruction.

#### `struct Instruction`

- `type` identifies whether the instruction applies a gate or measures.
- `gateName` stores a gate name for `ApplyGate`; measurement instructions use the original statement text but printing uses `MEASURE`.
- `qubitIndices` stores zero-based qubit indices rather than source names.
- `parameterValue` stores the evaluated numeric parameter; it is `0.0` when no parameter is used.

#### Code-generation functions

- `generateInstructions(const Program &program, const SemanticResult &result)` converts non-declaration AST statements into simulator instructions. It assumes semantic checking has already succeeded.
- `printInstructions(const std::vector<Instruction> &instructions)` prints the flat instruction list for the user.

### `include/simulator/simulator.h`

#### `class Simulator`

- `Simulator(int qubitCount)` allocates a state vector of size `2^qubitCount`, initializes it to the all-zero state, and seeds the random generator with `12345`.
- `run(const std::vector<Instruction> &instructions)` executes each generated instruction in order.
- `qubitCount` stores how many declared qubits exist.
- `state` stores complex amplitudes for every computational basis state.
- `randomGenerator` is the seeded `std::mt19937` used by measurement sampling.
- `applySingleQubitGate(...)` applies a supplied 2x2 complex matrix to a chosen qubit through amplitude pairs.
- `applyCNOT(int control, int target)` performs the controlled-not swap operation.
- `measure(const std::vector<int> &measuredQubits)` calculates probabilities, samples a full basis-state outcome, and prints bits for the requested measured qubits.
- `bitString(int stateIndex) const` converts a basis-state vector index into a displayed bit string in qubit-index order.

## 4. Source File Reference

### `src/compiler/lexer.cpp`

`classifyWord` receives a complete word and its location. It returns keyword tokens for `qubit`, `measure`, and `let`, a `Gate` token for `H`, `X`, `Z`, `CNOT`, and `RZ`, and an `Identifier` for every other valid word. This separation lets the parser distinguish statement starts without needing to compare raw strings everywhere.

`Lexer::Lexer` copies the supplied source text into the object's `source` field.

`Lexer::tokenize` is the main character scanner. It maintains `position`, `line`, and `column`. It skips spaces, tabs, and carriage returns; emits a `Newline` token for `\n`; emits individual tokens for punctuation and operators; scans a run of digits with an optional decimal point as a `Number`; and scans valid word characters as a word before passing it to `classifyWord`. An unrecognized character throws a line-numbered lexer error. At the end, it appends an `EndOfFile` token so the parser has an explicit stopping marker.

`tokenTypeName` is a `switch` that maps each `TokenType` enum member to the readable names shown in the driver output.

### `src/compiler/parser.cpp`

`Parser::Parser` stores a reference to the complete token list and initializes `current` to zero.

`peek` returns `tokens[current]` without changing the parser position. `advance` returns that token and then increments `current`. `check` compares `peek().type` to a requested type. `match` combines checking and consuming: it returns false without moving when the token does not match, and consumes one token when it does.

`skipNewlines` repeatedly calls `match(Newline)`. This is why blank lines do not become statements.

`expect` requires a particular token type. On failure it throws `Error: unexpected token '...' at line ...`; on success it calls `advance`.

`parseNameList` first requires an identifier, then saves `tokens[current - 1].text`. The `current - 1` expression is intentional: `expect` already advanced past the identifier, so the immediately previous token is the identifier that was just accepted. The function then repeats the same process for each comma followed by an identifier.

`parseExpression` implements the lowest arithmetic precedence. It starts with a complete term, then repeatedly consumes `+` or `-` and another term. Each loop creates a binary `Expression` whose two children are the previous left expression and the newly parsed right expression.

`parseTerm` has the same shape but handles `*` and `/`. Because `parseExpression` asks `parseTerm` for both sides of addition and subtraction, multiplication and division are fully grouped first.

`parseFactor` handles expression leaves and grouping. A number becomes a `NumberLiteral`; an identifier becomes an `Identifier`; a leading minus becomes a one-child `UnaryOperation`; and parentheses recursively call `parseExpression` before requiring `)`. If none of those forms match, it throws a parser error.

`parseStatement` selects a statement by its first token. For `qubit`, it creates a `QubitDeclaration`, saves the keyword's line number, and calls `parseNameList`. For `let`, it creates a `VariableDeclaration`, requires an identifier and `=`, and parses the right side with `parseExpression`. For a recognized gate or `measure`, it records the operation line and name, optionally parses a parenthesized parameter expression, then parses the operand list. The line is stored in every statement because semantic errors need the original source location.

`parse` skips initial blank lines, repeatedly parses one statement, requires a newline unless EOF follows immediately, and then skips extra blank lines. Its EOF exception permits an input file without a final newline.

### `src/compiler/semantic_checker.cpp`

The anonymous-namespace helper `findNameIndex` searches a vector of names with a plain loop and returns its index or `-1`. Both qubit and variable lookup use this simple approach, which keeps the symbol tables beginner-readable.

The anonymous-namespace helper `evaluateExpression` recursively evaluates an `Expression` into a `double`. Number nodes use `std::stod`. Identifier nodes look up the name in the already-built variable vectors; if absent, the function adds `Error: undefined variable 'name' used at line N`. Unary nodes evaluate their one child and negate it. Binary nodes evaluate both children first, allowing errors from both sides to be collected, then apply the stored operator. Division by zero is reported as a semantic error. This stage intentionally evaluates expressions to ordinary doubles rather than keeping symbolic values because the simulator needs a concrete numeric angle for `RZ`.

`checkSemantics` creates a `SemanticResult` and sizes `parameterValues` to match the statement vector. It walks statements in source order. For qubit declarations, it rejects duplicate names and otherwise appends names to `qubitNames`. For variable declarations, it rejects duplicate variable names, evaluates the declaration using only variables already stored, and adds a valid new name/value pair afterward. This ordering enforces the rule that a variable cannot be used before its `let` declaration.

For gate operations, `checkSemantics` reports every undeclared operand, checks required arity for `H`, `X`, `Z`, `RZ`, and `CNOT`, rejects a CNOT that uses the same qubit twice, requires a parameter for `RZ`, and rejects parameters for all other operation names. If a parameter exists, it evaluates the expression and saves the result at the same index as the statement. The function returns every collected error rather than throwing after the first one.

### `src/compiler/code_generator.cpp`

The anonymous-namespace helper `findQubitIndex` searches the validated declaration-order qubit list and returns a zero-based index. Semantic checking has already guaranteed that each use exists, so generation can safely use this result.

`generateInstructions` walks the AST in order. It skips qubit and variable declarations because they are compile-time information, not simulator actions. For each gate operation, it creates an `Instruction`, copies the gate name and evaluated parameter value, chooses `Measure` when the source name is `measure`, converts every source qubit name to an integer index, and appends the instruction. Converting names here separates the user-facing source language from the simulator's fast numeric representation.

`printInstructions` prints each instruction in a demonstration-friendly form. Gate instructions use `APPLY`; `RZ` includes its evaluated numeric value; measurement instructions use `MEASURE`; and all operands print as bracketed indices.
`optimizeInstructions` receives the validated circuit qubit count and performs a single left-to-right vector pass. It cancels only adjacent `ApplyGate` instructions with identical gate names and ordered qubit indices for H, X, Z, or CNOT. The supported single-qubit gates require exactly one index in `[0, qubitCount)`; CNOT requires two distinct indices in that range. These gates must have the unused `parameterValue` equal to zero. The output stack allows repeated cancellations to expose a new adjacent pair, while every measurement, RZ, unknown gate, and malformed instruction remains a barrier. It never commutes gates or removes operations across measurement.

The command-line driver runs this pass after generation and before simulation. It prints original and optimized unitary-gate counts separately from the original and optimized measurement counts. It does not report circuit depth.

### `src/simulator/simulator.cpp`

`Simulator::Simulator` stores the qubit count, seeds `randomGenerator` with `12345`, creates a complex state vector of size `1 << qubitCount`, fills it with zero amplitudes, and sets `state[0]` to `1 + 0i`. This represents the all-zero quantum state.

`applySingleQubitGate` applies a supplied 2x2 matrix to one qubit. It builds a bit mask for that qubit, visits every basis index whose selected bit is zero, finds the paired index with that bit set, saves the two old amplitudes, and replaces both amplitudes with the two matrix-multiplication results. Visiting only zero-bit indices ensures each pair is processed once.

`applyCNOT` builds masks for the control and target. It visits basis indices where the control bit is one and target bit is zero, then swaps that amplitude with the same state whose target bit is one. This is the standard controlled-not behavior and again processes each affected pair once.

`bitString` examines each qubit bit from index `0` through `qubitCount - 1` and appends `0` or `1` to a string. The displayed character order therefore follows QCompile's zero-based qubit-index order.

`measure` calculates `std::norm` for every amplitude and uses those values as a probability vector. It prints each non-negligible basis-state probability, constructs `std::discrete_distribution<int>`, samples a full-state outcome with the seeded generator, then prints only the bits named by the measurement instruction. The current implementation samples and displays a result; it does not modify the state vector after measurement.

`run` loops through generated instructions. Measurement instructions call `measure`. `H`, `X`, and `Z` provide their fixed 2x2 matrix entries to `applySingleQubitGate`. `RZ` computes phases `exp(-i theta/2)` and `exp(i theta/2)` from the evaluated parameter and applies the diagonal matrix. `CNOT` calls `applyCNOT` with its two indices.

### `src/apps/main.cpp`

`printExpression` recursively reconstructs an AST expression in readable infix form. It assigns precedence values to binary and unary nodes and adds parentheses only when an expression would otherwise bind less tightly than its parent. This is why printed AST expressions preserve important grouping such as `RZ((2 + 3) * 4)`.

`printProgram` walks statements and prints `QubitDecl`, `VariableDecl`, or `GateOp` lines. It uses `printExpression` for variable declarations and optional gate parameters.

`readSourceFile` opens the requested input file and reads it one line at a time with `std::getline`. Because `getline` removes newline characters and QCompile uses newlines as statement separators, it appends `"\n"` back after every line.

`countSourceLines` scans non-EOF tokens and returns the largest line number for the parsed-summary message.

`main` checks that a filename was supplied, reads the file, runs the lexer, prints tokens, runs the parser, prints the AST, and prints the statement summary. It then calls `checkSemantics`; if errors were collected, it prints all of them and returns without generating instructions. Otherwise it generates and prints instructions, constructs a simulator using the number of declared qubits, and runs it. The `try`/`catch` block prints lexer, parser, and file-reading exceptions cleanly.

## 5. How Translation Works, End to End

Consider the source line:

```text
H q0
```

1. **Lexer:** `Lexer::tokenize` reads `H` as a word. `classifyWord` recognizes it as a gate and emits `Gate("H", line 2)`. It reads `q0` as `Identifier("q0", line 2)` and then emits the line-ending `Newline("\\n", line 2)` token.
2. **Parser:** `parseStatement` sees the `Gate` token, saves `"H"` in `Statement::gate`, records its line, finds no `(` parameter, and calls `parseNameList`. The result is a `GateOperation` statement with `gate = "H"`, `qubits = ["q0"]`, and `hasParameter = false`.
3. **Semantic Checker:** When it reaches this statement, `checkSemantics` finds `q0` in `SemanticResult::qubitNames`, which was filled by an earlier `qubit q0` declaration. It checks that `H` has exactly one operand and that it has no parameter. If either check fails, it adds a line-numbered error and stops the later pipeline stages after all statements have been checked.
4. **Code Generator:** `generateInstructions` looks up `q0` in declaration order and obtains index `0`. It creates `Instruction { ApplyGate, "H", [0], 0.0 }`, printed as `APPLY H [0]`.
5. **Simulator:** `Simulator::run` recognizes `H` and calls `applySingleQubitGate` with the Hadamard matrix entries: `1/sqrt(2)`, `1/sqrt(2)`, `1/sqrt(2)`, and `-1/sqrt(2)`. For qubit zero, the simulator updates the amplitude pair for `|0>` and `|1>`. Starting from amplitude `1` for `|0>` and `0` for `|1>`, this produces equal amplitudes for `|0>` and `|1>`, giving probabilities `0.5` and `0.5` when measured.

This sequence demonstrates the central compiler idea: each stage changes representation while preserving meaning. Text becomes tokens, tokens become structured syntax, syntax is checked, checked syntax becomes simple machine-like instructions, and those instructions change the simulator state.

## 6. Design Decisions and Scope

QCompile deliberately keeps the DSL small so every implemented rule can be demonstrated and understood end to end.

- **No control flow:** There are no conditionals, loops, or runtime branches. This avoids adding execution-state rules beyond the current straight-line circuit model.
- **No registers or arrays:** Qubits are declared as individual names in a comma-separated list. Array indexing and register-wide operations would require additional grammar, semantic rules, and code-generation behavior.
- **No subroutines:** The language has no user-defined gate definitions, parameterized function calls, or reusable circuit blocks. This keeps source ordering and gate expansion simple.
- **Limited circuit optimization:** The pass cancels only adjacent identical H, X, Z, or same-ordered CNOT pairs. It does not reorder instructions, optimize RZ, cross measurements, or compute circuit depth.
- **Small state-vector simulation:** The simulator stores every basis-state amplitude, so its memory requirement grows as `2^n`. This is a textbook-correct and transparent approach for the small examples in this project, not a large-scale simulation strategy.
- **Numeric parameter evaluation only:** Variables and arithmetic expressions are evaluated to `double` values during semantic checking. Symbolic algebra is intentionally not implemented because the current simulator requires concrete values for `RZ` angles.
- **Measurement sampling without collapse:** Measurement calculates and samples probabilities but does not update the state vector after sampling. The provided sample programs measure at the end, which fits the current scope.

These are engineering choices, not incomplete parser cases: the implemented language is intentionally a compact circuit description language with a complete lexer-to-simulator path.

## 7. Fixed H₂ Energy Benchmark

`qcompile_h2` evaluates neutral singlet H₂ at a 0.735 Å bond distance in the STO-3G minimal basis with two electrons. It uses Qiskit Nature's two-qubit parity mapping reduced with particle counts `(1, 1)`. The fixed electronic Hamiltonian, in Hartree, is:

```text
-1.05237325 I + 0.39793742 Z0 - 0.39793742 Z1
-0.01128010 Z0 Z1 + 0.18093120 X0 X1
```

In QCompile's Pauli representation, an operator vector is indexed by qubit number and qubit zero is the least-significant state-vector bit. A Pauli expectation is calculated directly from the normalized complex state vector by applying its bit flips and phases; this path does not call or change the circuit language's `measure` operation.

`QubitHamiltonian` stores weighted Pauli terms and evaluates their sum against a `Simulator`. `pauliExpectation` also accepts a state vector directly, validates its dimension, finite amplitudes, and unit norm, and supports I, X, Y, and Z products. The exact reference path independently constructs the fixed 4-by-4 matrix from the Pauli terms and uses Jacobi rotations to find its lowest eigenvalue.

The one-parameter trial circuit prepares the encoded reference state, applies CNOT–RY–CNOT within the odd-parity two-state sector, and then returns to the encoded basis. RY is synthesized with the existing H and RZ gates. The golden-section search uses the bounded interval `[0, π/2]`, stops when its width is at most `1e-12` radians, and has a 128-iteration limit. The reported state is checked against the independent matrix minimum within `1e-10` Hartree.

The rounded Hamiltonian's exact minimum is approximately `-1.857275029288` Hartree; the Qiskit Nature electronic-structure reference is `-1.857275030202` Hartree. Nuclear repulsion (`+0.719968994449` Hartree) is reported separately and added once for total Born–Oppenheimer energy. The coefficient-rounding difference is less than `1e-9` Hartree.

Authoritative model and reference data: [Qiskit Nature electronic-structure tutorial](https://qiskit-community.github.io/qiskit-nature/tutorials/01_electronic_structure.html) and [qubit-mapping tutorial](https://qiskit-community.github.io/qiskit-nature/tutorials/06_qubit_mappers.html).

## 8. Geometry-Generated H₂ Input

The offline generator in `benchmarks/h2_geometry/generate_h2.py` accepts a bond length in Ångström and writes a versioned `QCOMPILE-H2-PAULI 1` Hamiltonian plus provenance JSON. It generates neutral singlet H₂ in STO-3G with RHF integrals, two spatial orbitals, no frozen orbitals, and Qiskit Nature `ParityMapper(num_particles=(1,1))`. The format records both Ångström and bohr distances, spin-orbital order, parity reduction, label endianness, Hartree units, the electronic-only Hamiltonian, and nuclear repulsion separately.

```bash
benchmarks/heh_plus/.venv-generator/bin/python -B \
  benchmarks/h2_geometry/generate_h2.py \
  --distance-angstrom 0.735 \
  --output /tmp/h2-0735.qham

make all
./qcompile_benchmarks --generated-h2 /tmp/h2-0735.qham
```

`loadGeneratedH2` accepts only this H₂ schema and rejects unsupported model metadata, invalid dimensions, non-finite values, missing or duplicate terms, and trailing content. The two-qubit parity basis uses q0 as its least-significant bit. The generated reference solver diagonalizes the two-dimensional odd-parity block at indices 1 and 2, checks decoupling from the complementary block, and rejects imaginary or asymmetric projected matrices when the residual exceeds `1e-11 Ha`. `evaluateGeneratedH2Energy` runs the existing one-parameter H₂ circuit against the loaded `QubitHamiltonian`; it verifies that the state remains normalized and inside the physical block. The bounded golden-section optimizer uses `[0, π/2]`, a `1e-12` radian interval threshold, and a 128-iteration cap. It verifies its exact-sector result against the generated singlet-FCI reference and rejects a variational state that fails the `1e-8 Ha` energy criterion or undercuts the exact result by more than `1e-10 Ha`.

The separate `qcompile_benchmarks --generated-h2 <file.qham>` path reports calculated electronic energy, physical-sector exact energy, generated FCI reference, variational error, nuclear repulsion, total energies, evaluation count, and convergence. Nuclear repulsion is added once and only when reporting total Born–Oppenheimer energy. The fixed catalog IDs and existing fixed H₂/HeH⁺ runners are unchanged. The offline generator and Psi4 validator are not C++ runtime dependencies; arbitrary molecules, bases, active spaces, and dashboard access are not supported.

Fresh generator outputs at 0.735 Å and 1.400 Å were validated with a determinant matrix assembled from MO integrals and an independent Psi4 FCI run. Test fixtures exercise the C++ parser, while the acceptance workflow separately runs the C++ executable on freshly generated files. Provenance and full generator details are documented in `benchmarks/h2_geometry/README.md`.

## 9. Fixed HeH⁺ Energy Benchmark

`qcompile_heh` evaluates only the frozen HeH⁺ benchmark generated in `benchmarks/heh_plus/`: He at `(0,0,0)` and H at `(0,0,1.463283)` bohr, STO-3G, charge `+1`, two electrons, singlet `S=0`, two retained spatial orbitals and four Jordan–Wigner spin orbitals. The mapper has no symmetry reduction. Spin-orbital/qubit ordering is `q0=MO0 α`, `q1=MO1 α`, `q2=MO0 β`, `q3=MO1 β`; bit `q0` is least significant. Frozen Pauli labels in the reference JSON are printed `q3` through `q0`, so the C++ loader reverses the label into its per-qubit operator vector. The stored electronic Hamiltonian excludes nuclear repulsion.

`hehHamiltonian` contains the 27 full-precision Pauli terms from `reference_data.json`. `hehSingletTrialCircuit(theta, phi)` starts in determinant index 5 and prepares

```text
cos(theta) |D0> + sin(theta) cos(phi) |S> + sin(theta) sin(phi) |D1>
|S> = (|6> + |9>) / sqrt(2)
```

where `|D0>` is index 5 and `|D1>` is index 10 in q0-least-significant occupation-bit order. The circuit compiles two CSF-plane Givens rotations to existing X, CNOT, H, and RZ instructions. Each RY used in the decomposition is synthesized as `RZ(-π/2), H, RZ(angle), H, RZ(+π/2)`. No new DSL gate or measurement behavior is introduced. `Simulator::stateVector()` is a read-only reference used by validation code.

`exactHeHSingletEnergy` constructs the 16-by-16 complex matrix from the Pauli terms, projects it onto the three singlet CSFs, and uses the shared dimension-checked Jacobi eigensolver. This excludes the Ms=0 triplet determinant combination. The variational optimizer uses a deterministic 3-by-4 grid of starting points, bounded coordinate-wise golden-section minimization, at most eight sweeps per start, 48 iterations per coordinate, `1e-8` radian interval stopping, and a 10,000 energy-evaluation ceiling. The reported state must lie within `1e-8` Ha of the exact singlet minimum and cannot undercut it by more than `1e-10` Ha. The exact energy is also required to agree with the independently validated Psi4 FCI reference within `1e-10` Ha.

The external reference reports electronic FCI energy `−4.218264154490385` Ha, nuclear repulsion `+1.3667896093920315` Ha, and total Born–Oppenheimer energy `−2.8514745450983536` Ha. Nuclear repulsion is added only when producing total energy. Frozen generation settings, integrals, coefficient provenance, and matrix-equivalence checks are in `benchmarks/heh_plus/reference_data.json`; the independent Psi4 result and calculation log are retained alongside it. Python/Qiskit Nature/PySCF/Psi4 are offline reference-generation tools, not runtime or C++ build dependencies.

This is an educational implementation of established Pauli-Hamiltonian, state-vector VQE-style optimization, and exact-diagonalization methods. The generated-H₂ path is restricted to the documented H₂/STO-3G/parity model. It does not provide arbitrary-molecule coverage, quantum advantage, or a novel algorithm.

## 10. Fixed Benchmark Catalog and Matrix Input Contract

`benchmarkCatalog`, `findBenchmark`, and `runBenchmark` provide a fixed C++17 catalog with stable IDs `h2` and `heh_plus`. `qcompile_benchmarks --list` lists the registered entries; passing either ID runs its existing molecule-specific optimizer and reports a common `BenchmarkResult`. `BenchmarkMetadata` records fixed geometry, charge, electron/spin/orbital/qubit counts, basis, mapping, ordering, exact-sector description, reference provenance, and energy conventions. The result keeps optimized parameters and count, electronic energy, nuclear repulsion, total energy, exact electronic energy, external reference energy, variational error, evaluation count, and convergence status. Catalog completion is marked converged only after the underlying optimizer returns successfully; optimizer failures propagate as errors. The separate `qcompile_h2` and `qcompile_heh` entry points and specialized calculation paths remain unchanged. The generated-H₂ file path is separate from the fixed catalog; it does not register arbitrary geometries or any additional molecule.

`lowestSymmetricEigenvalue` accepts a compile-time-sized real matrix with dimension at least two. It rejects non-finite entries and rejects asymmetry when `|aᵢⱼ-aⱼᵢ| > 1e-12 × max(1, |aᵢⱼ|, |aⱼᵢ|)`. The scale floor keeps the contract well-defined for small matrix elements; the tolerance exceeds the observed double-precision symmetry residuals of the frozen matrices. The Jacobi rotation method and its `1e-14` off-diagonal convergence threshold are unchanged. A matrix that does not converge within the existing rotation limit still throws.

## 11. Local Dashboard

`node dashboard/server.js` serves a local dashboard at `http://127.0.0.1:4173`; build `qcompile_benchmarks` first with `make all`. The Node.js server uses built-in modules only, binds exclusively to loopback, serves the local frontend assets, and runs the existing C++ benchmark executable through `execFile` with an allowlist of `h2` and `heh_plus` IDs. It does not pass user-provided values to a shell.

`qcompile_benchmarks --json-list` emits fixed catalog metadata, while `qcompile_benchmarks --json <benchmark-id>` emits one structured run result. The pre-existing `--list` and human-readable `<benchmark-id>` modes retain their output. The dashboard renders current process results, clears the prior result when a new run begins or fails, and shows electronic, exact-sector, documented reference, nuclear-repulsion, and total Born–Oppenheimer energies separately. Total energy is electronic plus nuclear repulsion exactly once. Missing numeric values are treated as unavailable rather than zero.

The UI reuses the Tabler Core v1.4.0 component stylesheet under its MIT license, stored locally in `dashboard/vendor/`; it does not load a CDN or introduce an npm runtime dependency. The dashboard still exposes only the two fixed catalog models; the separate generated-H₂ CLI is not wired into the dashboard. No arbitrary-molecule generation, molecular visualization, or quantum-advantage claim is provided.

## 12. Native H₂ Bonding Laboratory

`make qcompile_lab` builds the separate native desktop target using the installed Qt 6 Widgets development package and `pkg-config`. `make all` continues to build the pre-existing compiler and benchmark executables without adding Qt to their dependency path. Its native header, dark workbench treatment, card grouping, and status hierarchy adapt the existing local Tabler Core v1.4.0 dashboard’s visual direction; the MIT-licensed CSS remains a web asset, not a native runtime dependency. The implementation uses Qt widgets and QPainter-projected 3D geometry; raylib and Dear ImGui are not installed in the inspected environment and are not required.

Run `make all`, `make qcompile_lab`, then `./qcompile_lab`. The preview bond length is adjustable from 0.400 to 3.000 Å. The viewer distinguishes preview geometry from the latest completed calculation and supports mouse-drag rotation and wheel zoom. Nuclei and the bond are schematic representations, not classical trajectories or an electron-density model.

Each calculation creates a private `QTemporaryDir`, invokes the existing offline Python generator with a `QProcess` argument list, then invokes `qcompile_benchmarks --generated-h2-json` with another argument list. The versioned JSON output is checked for schema, geometry, supported model metadata, finite values, convergence, and energy identities; the existing human-readable `--generated-h2` CLI mode is unchanged. No shell command string is constructed. Generator/solver failures, cancellation, malformed output, and stale request IDs do not become current results. Temporary Hamiltonians and sidecars are removed with the private directory.

`h2TrialInstructions(parameter)` is the single gate-list source used by both H₂ energy evaluation and the circuit panel. Playback selects a gate prefix; the inspector constructs a fresh existing `Simulator`, executes that real prefix, and reads its actual `stateVector()`. It shows computational basis labels in `q1q0` order, amplitudes, probabilities `|amplitude|²`, and norm. This does not change sampling-only, non-collapsing DSL measurement semantics.

The five-point sweep uses fresh generator/solver runs at 0.735, 0.900, 1.100, 1.250, and 1.400 Å. The plot displays only actual calculated total-energy points and corresponding physical-sector exact-reference points; it does not draw an interpolated curve. The guided prediction compares the two endpoint totals and identifies the minimum only among the sweep’s sampled geometries. The one-parameter ansatz is not claimed to remain accurate throughout dissociation.

`./qcompile_lab --smoke-test` exercises camera rotation/zoom, bond preview updates, circuit stepping/play/pause, two individual geometry runs, and the five-point sweep, then saves a screenshot to `/tmp/qcompile-lab-smoke.png` without overwriting an existing file. Optional executable overrides are `QCOMPILE_LAB_PYTHON`, `QCOMPILE_LAB_GENERATOR`, and `QCOMPILE_LAB_BENCHMARKS`. The app adds no chemistry-package runtime dependency, web frontend, or support for molecules other than H₂.
