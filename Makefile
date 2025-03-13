# Compiler
CXX = g++-10

# Compiler flags for both builds (serial and parallel)
CXXFLAGS = -std=c++11 -g -Wall -Iinclude -O0

# Libraries for your main executable
LDLIBS = -lboost_program_options

# Source files for production code
SRC = $(wildcard src/*.cpp)

# Object files directory
OBJDIR = build

# Production object files and executable for the serial build
OBJ = $(patsubst src/%.cpp, $(OBJDIR)/%.o, $(SRC))
TARGET_SERIAL = $(OBJDIR)/md

# Production object files and executable for the parallel build
PAR_CXXFLAGS = $(CXXFLAGS) -fopenmp
PAR_OBJ = $(patsubst src/%.cpp, $(OBJDIR)/%.par.o, $(SRC))
TARGET_PAR = $(OBJDIR)/mdpar

# Test source file(s)
TEST_SRC = tests/unittests.cpp
TEST_OBJ = $(patsubst tests/%.cpp, $(OBJDIR)/%.o, $(TEST_SRC))
TEST_EXEC = $(OBJDIR)/unittests

# Default target builds the serial main executable
all: $(TARGET_SERIAL)

# Link production object files to create the serial executable
$(TARGET_SERIAL): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# Build test executable by linking test and production object files.
# Exclude main.o from the production objects so the Boost.Test main is used.
TEST_PROD_OBJ = $(filter-out $(OBJDIR)/main.o, $(OBJ))
$(TARGET_TEST): $(TEST_OBJ) $(TEST_PROD_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS) -lboost_unit_test_framework

# Build the parallel executable (mdpar) from production sources compiled with OpenMP flags.
$(TARGET_PAR): $(PAR_OBJ)
	$(CXX) $(PAR_CXXFLAGS) -o $@ $^ $(LDLIBS)

# Compile production source files into object files (serial)
$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile production source files into parallel object files (with -fopenmp)
$(OBJDIR)/%.par.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(PAR_CXXFLAGS) -c $< -o $@

# Compile test source files into object files
$(OBJDIR)/%.o: tests/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create the build directory if it doesn't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Target to run tests.
.PHONY: test
test: $(TARGET_TEST)
	./$(TARGET_TEST)

# Target to run the parallel executable.
.PHONY: mdpar
mdpar: $(TARGET_PAR)
	./$(TARGET_PAR)

# Clean up build files
clean:
	rm -rf $(OBJDIR)

# Phony targets for all, test, and clean
.PHONY: all test clean
