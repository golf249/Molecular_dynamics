# Compiler
CXX = g++-10

# Compiler flags
CXXFLAGS = -std=c++11 -g -Wall -Iinclude -O0

# Libraries for your main executable
LDLIBS = -lboost_program_options

# Source files for production code
SRC = $(wildcard src/*.cpp)

# Object files directory
OBJDIR = build

# Production object files and executable
OBJ = $(patsubst src/%.cpp, $(OBJDIR)/%.o, $(SRC))
EXEC = $(OBJDIR)/md

# Test source file(s)
TEST_SRC = tests/unittests.cpp
TEST_OBJ = $(patsubst tests/%.cpp, $(OBJDIR)/%.o, $(TEST_SRC))
TEST_EXEC = $(OBJDIR)/unittests

# Default target builds the main executable
all: $(EXEC)

# Link production object files to create the executable
$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# Build test executable by linking test and production object files.
# Exclude main.o from the production objects so the Boost.Test main is used.
TEST_PROD_OBJ = $(filter-out $(OBJDIR)/main.o, $(OBJ))
$(TEST_EXEC): $(TEST_OBJ) $(TEST_PROD_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS) -lboost_unit_test_framework

# Compile production source files into object files
$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile test source files into object files
$(OBJDIR)/%.o: tests/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create the build directory if it doesn't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Target to run tests.
.PHONY: test
test: $(TEST_EXEC)
	./$(TEST_EXEC)

# Clean up build files
clean:
	rm -rf $(OBJDIR)

# Phony targets for all, test, and clean
.PHONY: all test clean
