# Compiler and flags used to build QCompile.
CXX ?= g++
BASE_CXXFLAGS = -std=c++17 -Wall -Wextra -pedantic

# Build configuration: release (default, with -O2) or debug (with -g -O0)
BUILD ?= release
ifeq ($(BUILD),debug)
  OPTFLAGS = -g -O0
  CFG = debug
else
  OPTFLAGS = -O2
  CFG = release
endif

CXXFLAGS = $(BASE_CXXFLAGS) $(OPTFLAGS)

# Header location and dependency generation flags.
INCLUDES = -Iinclude
DEPFLAGS = -MMD -MP

# Build output directory isolated by configuration.
BUILD_BASE = build
BUILD_DIR = $(BUILD_BASE)/$(CFG)
BIN_DIR = $(BUILD_DIR)/bin

# Targets.
TARGET = qcompile
H2_TARGET = qcompile_h2
HEH_TARGET = qcompile_heh
BENCHMARK_TARGET = qcompile_benchmarks
LAB_TARGET = qcompile_lab
TEST_TARGET = qcompile_tests

# Source lists.
SOURCES = src/apps/main.cpp src/compiler/lexer.cpp src/compiler/parser.cpp src/compiler/semantic_checker.cpp src/compiler/code_generator.cpp src/compiler/compiler_facade.cpp src/simulator/simulator.cpp src/simulator/simulation_session.cpp src/chemistry/hamiltonian.cpp src/chemistry/h2_energy.cpp src/experiment/experiment_bundle.cpp
H2_SOURCES = src/apps/h2_main.cpp src/simulator/simulator.cpp src/chemistry/hamiltonian.cpp src/chemistry/h2_energy.cpp
HEH_SOURCES = src/apps/heh_main.cpp src/simulator/simulator.cpp src/chemistry/hamiltonian.cpp src/chemistry/heh_energy.cpp
BENCHMARK_SOURCES = src/apps/benchmark_main.cpp src/chemistry/benchmark_catalog.cpp src/chemistry/h2_energy.cpp src/chemistry/h2_geometry.cpp src/chemistry/heh_energy.cpp src/simulator/simulator.cpp src/chemistry/hamiltonian.cpp
LAB_SOURCES = src/gui/h2_lab.cpp src/gui/circuit_view.cpp src/compiler/compiler_facade.cpp src/compiler/lexer.cpp src/compiler/parser.cpp src/compiler/semantic_checker.cpp src/compiler/code_generator.cpp src/simulator/simulator.cpp src/simulator/simulation_session.cpp src/chemistry/hamiltonian.cpp src/chemistry/h2_energy.cpp src/chemistry/h2_geometry.cpp src/chemistry/heh_energy.cpp src/experiment/experiment_bundle.cpp
TEST_SOURCES = tests/qcompile_tests.cpp tests/test_support.cpp tests/test_lexer.cpp tests/test_parser.cpp tests/test_semantics.cpp tests/test_optimizer.cpp tests/test_simulator.cpp tests/test_chemistry.cpp tests/test_compiler_facade.cpp tests/test_simulation_session.cpp tests/test_workbench_integration.cpp tests/test_experiment_bundle.cpp src/compiler/lexer.cpp src/compiler/parser.cpp src/compiler/semantic_checker.cpp src/compiler/code_generator.cpp src/compiler/compiler_facade.cpp src/simulator/simulator.cpp src/simulator/simulation_session.cpp src/chemistry/hamiltonian.cpp src/chemistry/h2_energy.cpp src/chemistry/h2_geometry.cpp src/chemistry/heh_energy.cpp src/chemistry/benchmark_catalog.cpp src/experiment/experiment_bundle.cpp

# Map source files to object files in BUILD_DIR.
OBJECTS = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
H2_OBJECTS = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(H2_SOURCES))
HEH_OBJECTS = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(HEH_SOURCES))
BENCHMARK_OBJECTS = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(BENCHMARK_SOURCES))
LAB_OBJECTS = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(LAB_SOURCES))
TEST_OBJECTS = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_SOURCES))

# All object files for current configuration dependency inclusion.
ALL_OBJECTS = $(sort $(OBJECTS) $(H2_OBJECTS) $(HEH_OBJECTS) $(BENCHMARK_OBJECTS) $(LAB_OBJECTS) $(TEST_OBJECTS))
ALL_DEPS = $(ALL_OBJECTS:.o=.d)

# Configuration tracker file to force relinking/resymlinking top-level binaries when switching BUILD.
CONFIG_FILE = $(BUILD_BASE)/.current_config
CONFIG_SWITCH = $(BUILD_BASE)/.config_switched_to_$(CFG)

$(CONFIG_SWITCH):
	@mkdir -p $(BUILD_BASE)
	@echo "Switching build configuration to $(CFG)..."
	@echo $(CFG) > $(CONFIG_FILE)
	@rm -f $(BUILD_BASE)/.config_switched_to_*
	@touch $(CONFIG_SWITCH)
	@rm -f $(TARGET) $(H2_TARGET) $(HEH_TARGET) $(BENCHMARK_TARGET) $(TEST_TARGET) $(LAB_TARGET)

# Default rule builds standard production CLI targets.
all: $(TARGET) $(H2_TARGET) $(HEH_TARGET) $(BENCHMARK_TARGET)

$(TARGET): $(BIN_DIR)/$(TARGET) $(CONFIG_SWITCH)
	ln -sf $(BIN_DIR)/$(TARGET) $(TARGET)

$(H2_TARGET): $(BIN_DIR)/$(H2_TARGET) $(CONFIG_SWITCH)
	ln -sf $(BIN_DIR)/$(H2_TARGET) $(H2_TARGET)

$(HEH_TARGET): $(BIN_DIR)/$(HEH_TARGET) $(CONFIG_SWITCH)
	ln -sf $(BIN_DIR)/$(HEH_TARGET) $(HEH_TARGET)

$(BENCHMARK_TARGET): $(BIN_DIR)/$(BENCHMARK_TARGET) $(CONFIG_SWITCH)
	ln -sf $(BIN_DIR)/$(BENCHMARK_TARGET) $(BENCHMARK_TARGET)

$(LAB_TARGET): $(BIN_DIR)/$(LAB_TARGET) $(CONFIG_SWITCH)
	ln -sf $(BIN_DIR)/$(LAB_TARGET) $(LAB_TARGET)

$(TEST_TARGET): $(BIN_DIR)/$(TEST_TARGET) $(CONFIG_SWITCH)
	ln -sf $(BIN_DIR)/$(TEST_TARGET) $(TEST_TARGET)
# Configuration-specific binaries in BIN_DIR
$(BIN_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

$(BIN_DIR)/$(H2_TARGET): $(H2_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(H2_OBJECTS) -o $@

$(BIN_DIR)/$(HEH_TARGET): $(HEH_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(HEH_OBJECTS) -o $@

$(BIN_DIR)/$(BENCHMARK_TARGET): $(BENCHMARK_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(BENCHMARK_OBJECTS) -o $@

$(BIN_DIR)/$(LAB_TARGET): $(LAB_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -fPIC -no-pie -pthread $(LAB_OBJECTS) -o $@ $$(pkg-config --libs Qt6Widgets)

$(BIN_DIR)/$(TEST_TARGET): $(TEST_OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(TEST_OBJECTS) -o $@
# Pattern rule for compiling C++ translation units with dependency tracking.
$(BUILD_DIR)/src/gui/%.o: src/gui/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -fPIC $(INCLUDES) $$(pkg-config --cflags Qt6Widgets) $(DEPFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEPFLAGS) -c $< -o $@

test: $(TEST_TARGET)
	./$(TEST_TARGET)

# Include generated dependency files if present.
-include $(ALL_DEPS)

# Running "make clean" removes generated executables and build objects.
clean:
	rm -f $(TARGET) $(H2_TARGET) $(HEH_TARGET) $(BENCHMARK_TARGET) $(TEST_TARGET) $(LAB_TARGET)
	rm -rf $(BUILD_BASE)

.PHONY: all clean test
