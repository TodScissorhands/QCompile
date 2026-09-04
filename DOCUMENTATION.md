# QCompile Technical Documentation

## 1. Project Overview

QCompile is a small C++ compiler for a custom domain-specific language (DSL) that describes quantum circuits. It compiles a text program into instructions for a classical state-vector simulator. It does not communicate with real quantum hardware.

The implemented pipeline has five stages:

1. **Lexer** — Reads source characters and creates tokens such as `KeywordQubit`, `Gate`, `Identifier`, `Number`, `Comma`, `Newline`, and arithmetic operators. Tokens preserve their original text and source location.
2. **Parser** — Reads the token list and creates an abstract syntax tree (AST). The AST represents qubit declarations, variable declarations, gate operations, and arithmetic expression trees.
3. **Semantic Checker** — Validates information that cannot be checked by token shape alone: qubits must be declared before use, gates must receive the right number of qubits, CNOT must use different control and target qubits, and parameter expressions must evaluate to numbers using already-declared variables.
4. **Code Generator** — Converts valid AST statements into a flat vector of instructions. It changes qubit names such as `q0` into zero-based integer indices such as `0`.
5. **Simulator** — Maintains a vector of complex amplitudes, applies the generated gates, calculates measurement probabilities, and samples a measurement result with a fixed random seed.

The driver in `src/main.cpp` connects these stages in order. If parsing fails, later stages do not run. If semantic checking reports one or more errors, code generation and simulation do not run.

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

### `include/ast.h`

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

### `include/lexer.h`

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

### `include/parser.h`

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

### `include/semantic_checker.h`

#### `struct SemanticResult`

`SemanticResult` returns both errors and the information needed by code generation.

- `errors` collects all semantic error strings instead of stopping after the first one.
- `qubitNames` stores declared qubits in declaration order, which becomes code-generation index order.
- `variableNames` stores valid variable names in declaration order.
- `variableValues` stores the corresponding evaluated `double` values at the same indices as `variableNames`.
- `parameterValues` stores evaluated parameter values aligned with `Program::statements`; non-parameter statements keep `0.0`.

#### `checkSemantics`

- `SemanticResult checkSemantics(const Program &program)` walks the AST in order, collects all semantic errors, evaluates valid variables and parameters, and returns the result.

### `include/code_generator.h`

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

### `include/simulator.h`

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

### `src/lexer.cpp`

`classifyWord` receives a complete word and its location. It returns keyword tokens for `qubit`, `measure`, and `let`, a `Gate` token for `H`, `X`, `Z`, `CNOT`, and `RZ`, and an `Identifier` for every other valid word. This separation lets the parser distinguish statement starts without needing to compare raw strings everywhere.

`Lexer::Lexer` copies the supplied source text into the object's `source` field.

`Lexer::tokenize` is the main character scanner. It maintains `position`, `line`, and `column`. It skips spaces, tabs, and carriage returns; emits a `Newline` token for `\n`; emits individual tokens for punctuation and operators; scans a run of digits with an optional decimal point as a `Number`; and scans valid word characters as a word before passing it to `classifyWord`. An unrecognized character throws a line-numbered lexer error. At the end, it appends an `EndOfFile` token so the parser has an explicit stopping marker.

`tokenTypeName` is a `switch` that maps each `TokenType` enum member to the readable names shown in the driver output.

### `src/parser.cpp`

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

### `src/semantic_checker.cpp`

The anonymous-namespace helper `findNameIndex` searches a vector of names with a plain loop and returns its index or `-1`. Both qubit and variable lookup use this simple approach, which keeps the symbol tables beginner-readable.

The anonymous-namespace helper `evaluateExpression` recursively evaluates an `Expression` into a `double`. Number nodes use `std::stod`. Identifier nodes look up the name in the already-built variable vectors; if absent, the function adds `Error: undefined variable 'name' used at line N`. Unary nodes evaluate their one child and negate it. Binary nodes evaluate both children first, allowing errors from both sides to be collected, then apply the stored operator. Division by zero is reported as a semantic error. This stage intentionally evaluates expressions to ordinary doubles rather than keeping symbolic values because the simulator needs a concrete numeric angle for `RZ`.

`checkSemantics` creates a `SemanticResult` and sizes `parameterValues` to match the statement vector. It walks statements in source order. For qubit declarations, it rejects duplicate names and otherwise appends names to `qubitNames`. For variable declarations, it rejects duplicate variable names, evaluates the declaration using only variables already stored, and adds a valid new name/value pair afterward. This ordering enforces the rule that a variable cannot be used before its `let` declaration.

For gate operations, `checkSemantics` reports every undeclared operand, checks required arity for `H`, `X`, `Z`, `RZ`, and `CNOT`, rejects a CNOT that uses the same qubit twice, requires a parameter for `RZ`, and rejects parameters for all other operation names. If a parameter exists, it evaluates the expression and saves the result at the same index as the statement. The function returns every collected error rather than throwing after the first one.

### `src/code_generator.cpp`

The anonymous-namespace helper `findQubitIndex` searches the validated declaration-order qubit list and returns a zero-based index. Semantic checking has already guaranteed that each use exists, so generation can safely use this result.

`generateInstructions` walks the AST in order. It skips qubit and variable declarations because they are compile-time information, not simulator actions. For each gate operation, it creates an `Instruction`, copies the gate name and evaluated parameter value, chooses `Measure` when the source name is `measure`, converts every source qubit name to an integer index, and appends the instruction. Converting names here separates the user-facing source language from the simulator's fast numeric representation.

`printInstructions` prints each instruction in a demonstration-friendly form. Gate instructions use `APPLY`; `RZ` includes its evaluated numeric value; measurement instructions use `MEASURE`; and all operands print as bracketed indices.

### `src/simulator.cpp`

`Simulator::Simulator` stores the qubit count, seeds `randomGenerator` with `12345`, creates a complex state vector of size `1 << qubitCount`, fills it with zero amplitudes, and sets `state[0]` to `1 + 0i`. This represents the all-zero quantum state.

`applySingleQubitGate` applies a supplied 2x2 matrix to one qubit. It builds a bit mask for that qubit, visits every basis index whose selected bit is zero, finds the paired index with that bit set, saves the two old amplitudes, and replaces both amplitudes with the two matrix-multiplication results. Visiting only zero-bit indices ensures each pair is processed once.

`applyCNOT` builds masks for the control and target. It visits basis indices where the control bit is one and target bit is zero, then swaps that amplitude with the same state whose target bit is one. This is the standard controlled-not behavior and again processes each affected pair once.

`bitString` examines each qubit bit from index `0` through `qubitCount - 1` and appends `0` or `1` to a string. The displayed character order therefore follows QCompile's zero-based qubit-index order.

`measure` calculates `std::norm` for every amplitude and uses those values as a probability vector. It prints each non-negligible basis-state probability, constructs `std::discrete_distribution<int>`, samples a full-state outcome with the seeded generator, then prints only the bits named by the measurement instruction. The current implementation samples and displays a result; it does not modify the state vector after measurement.

`run` loops through generated instructions. Measurement instructions call `measure`. `H`, `X`, and `Z` provide their fixed 2x2 matrix entries to `applySingleQubitGate`. `RZ` computes phases `exp(-i theta/2)` and `exp(i theta/2)` from the evaluated parameter and applies the diagonal matrix. `CNOT` calls `applyCNOT` with its two indices.

### `src/main.cpp`

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
- **No circuit optimization:** The code generator preserves source operation order. It does not cancel gates, merge rotations, reduce gate count, or map gates to hardware connectivity.
- **Small state-vector simulation:** The simulator stores every basis-state amplitude, so its memory requirement grows as `2^n`. This is a textbook-correct and transparent approach for the small examples in this project, not a large-scale simulation strategy.
- **Numeric parameter evaluation only:** Variables and arithmetic expressions are evaluated to `double` values during semantic checking. Symbolic algebra is intentionally not implemented because the current simulator requires concrete values for `RZ` angles.
- **Measurement sampling without collapse:** Measurement calculates and samples probabilities but does not update the state vector after sampling. The provided sample programs measure at the end, which fits the current scope.

These are engineering choices, not incomplete parser cases: the implemented language is intentionally a compact circuit description language with a complete lexer-to-simulator path.
