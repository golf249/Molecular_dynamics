# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Iinclude -lboost_program_options

# Source files
SRC = $(wildcard src/*.cpp)

# Object files
OBJ = $(SRC:.cpp=.o)

# Executable name
EXEC = molecular_dynamics

# Default target
all: $(EXEC)

# Link object files to create the executable
$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJ) $(EXEC)

# Phony targets
.PHONY: all clean
