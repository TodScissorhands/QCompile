# Compiler and flags used to build QCompile.
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -pedantic

# Header location, source files, and executable targets.
INCLUDES = -Iinclude
SOURCES = src/apps/main.cpp src/compiler/lexer.cpp src/compiler/parser.cpp src/compiler/semantic_checker.cpp src/compiler/code_generator.cpp src/simulator/simulator.cpp src/chemistry/hamiltonian.cpp src/chemistry/h2_energy.cpp
HEADERS = $(wildcard include/*/*.h)
TEST_SOURCES = tests/qcompile_tests.cpp src/compiler/lexer.cpp src/compiler/parser.cpp src/compiler/semantic_checker.cpp src/compiler/code_generator.cpp src/simulator/simulator.cpp src/chemistry/hamiltonian.cpp src/chemistry/h2_energy.cpp src/chemistry/h2_geometry.cpp src/chemistry/heh_energy.cpp src/chemistry/benchmark_catalog.cpp
TARGET = qcompile
H2_TARGET = qcompile_h2
HEH_TARGET = qcompile_heh
BENCHMARK_TARGET = qcompile_benchmarks
LAB_TARGET = qcompile_lab
TEST_TARGET = qcompile_tests
H2_SOURCES = src/apps/h2_main.cpp src/simulator/simulator.cpp src/chemistry/hamiltonian.cpp src/chemistry/h2_energy.cpp
HEH_SOURCES = src/apps/heh_main.cpp src/simulator/simulator.cpp src/chemistry/hamiltonian.cpp src/chemistry/heh_energy.cpp
BENCHMARK_SOURCES = src/apps/benchmark_main.cpp src/chemistry/benchmark_catalog.cpp src/chemistry/h2_energy.cpp src/chemistry/h2_geometry.cpp src/chemistry/heh_energy.cpp src/simulator/simulator.cpp src/chemistry/hamiltonian.cpp
LAB_SOURCES = src/gui/h2_lab.cpp src/simulator/simulator.cpp src/chemistry/hamiltonian.cpp src/chemistry/h2_energy.cpp src/chemistry/h2_geometry.cpp

# Running "make" builds the executable.
all: $(TARGET) $(H2_TARGET) $(HEH_TARGET) $(BENCHMARK_TARGET)

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCES) -o $(TARGET)

$(H2_TARGET): $(H2_SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(H2_SOURCES) -o $(H2_TARGET)

$(HEH_TARGET): $(HEH_SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(HEH_SOURCES) -o $(HEH_TARGET)

$(BENCHMARK_TARGET): $(BENCHMARK_SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(BENCHMARK_SOURCES) -o $(BENCHMARK_TARGET)

$(LAB_TARGET): $(LAB_SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) -fPIC -no-pie $(INCLUDES) $$(pkg-config --cflags Qt6Widgets) $(LAB_SOURCES) -pthread -o $(LAB_TARGET) $$(pkg-config --libs Qt6Widgets)

$(TEST_TARGET): $(TEST_SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(TEST_SOURCES) -o $(TEST_TARGET)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

# Running "make clean" removes generated executables.
clean:
	rm -f $(TARGET) $(H2_TARGET) $(HEH_TARGET) $(BENCHMARK_TARGET) $(TEST_TARGET) $(LAB_TARGET)

.PHONY: all clean test
