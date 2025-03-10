# Compiler
CXX = g++-10

# Compiler flags
CXXFLAGS = -std=c++11 -g -Wall -Iinclude -O0

# Libraries
LDLIBS = -lboost_program_options

# Source files
SRC = $(wildcard src/*.cpp)

# Object files directory
OBJDIR = build

# Object files
OBJ = $(patsubst src/%.cpp, $(OBJDIR)/%.o, $(SRC))

# Executable name
EXEC = $(OBJDIR)/md

# Default target
all: $(EXEC)

# Link object files to create the executable
$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# Compile source files into object files
$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create the build directory if it doesn't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Clean up build files
clean:
	rm -rf $(OBJDIR)

# Phony targets
.PHONY: all clean
