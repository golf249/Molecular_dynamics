# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -O0 -Iinclude

# Libraries
LDLIBS = -lboost_program_options

# Source files
SRC = $(wildcard src/*.cpp)

# Object files
OBJ = $(SRC:.cpp=.o)

# Executable name
EXEC = md

# Default target
all: $(EXEC)

# Link object files to create the executable
$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJ) $(EXEC) *.txt

# Phony targets
.PHONY: all clean
