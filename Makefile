# Compiler and flags used to build QCompile.
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -pedantic

# Header location, source files, and the output executable.
INCLUDES = -Iinclude
SOURCES = src/main.cpp src/lexer.cpp src/parser.cpp src/semantic_checker.cpp src/code_generator.cpp src/simulator.cpp
TARGET = qcompile

# Running "make" builds the executable.
all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SOURCES) -o $(TARGET)

# Running "make clean" removes the generated executable.
clean:
	rm -f $(TARGET)
