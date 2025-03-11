# Compiler
CXX = g++-10

# Compiler flags
CXXFLAGS = -std=c++11 -g -Wall -Iinclude -O0

# Libraries
LDLIBS = -lboost_program_options

# Boost Test library
BOOST_LIBS = -lboost_unit_test_framework

# Source files
SRC = $(wildcard src/*.cpp)

# Test source files
TEST_SRC = $(wildcard tests/*.cpp)

# Object files directory
OBJDIR = build

# Object files
OBJ = $(patsubst src/%.cpp, $(OBJDIR)/%.o, $(SRC))

# Test object files
TEST_OBJ = $(patsubst tests/%.cpp, $(OBJDIR)/%.o, $(TEST_SRC))

# Executable name
EXEC = $(OBJDIR)/md

# Test executable name
TEST_EXEC = $(OBJDIR)/unittests

# Default target
all: $(EXEC)

# Link object files to create the executable
$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# Compile source files into object files
$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile test source files into object files
$(OBJDIR)/%.o: tests/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link test object files to create the test executable
$(TEST_EXEC): $(OBJ) $(TEST_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS) $(BOOST_LIBS)

# Run unit tests
unittests: $(TEST_EXEC)
	./$(TEST_EXEC)

# Create the build directory if it doesn't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Clean up build files
clean:
	rm -rf $(OBJDIR)

# Phony targets
.PHONY: all clean unittests
