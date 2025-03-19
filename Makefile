# Compiler
CXX = g++-10
NVCC = nvcc

# Compiler flags for both builds (serial and parallel)
CXXFLAGS = -std=c++11 -g -Wall -Iinclude -O3
NVCCFLAGS = -std=c++11 -g -Xcompiler -Wall -Iinclude_cuda -O3

# Libraries for your main executable
LDLIBS = -lboost_program_options

# Source files for production code
SRC = $(wildcard src/*.cpp)
CUDA_SRC = $(wildcard src_cuda/*.cu)
COMMON_SRC = src/Particle.cpp src/writeFile.cpp src/main.cpp

# Object files directory
OBJDIR = build

# Production object files and executable for the serial build
OBJ = $(patsubst src/%.cpp, $(OBJDIR)/%.o, $(SRC))
TARGET_SERIAL = $(OBJDIR)/md

# Production object files and executable for the parallel build
PAR_CXXFLAGS = $(CXXFLAGS) -fopenmp
PAR_OBJ = $(patsubst src/%.cpp, $(OBJDIR)/%.par.o, $(SRC))
TARGET_PAR = $(OBJDIR)/mdpar

# Production object files and executable for the CUDA build
CUDA_OBJ = $(patsubst src_cuda/%.cu, $(OBJDIR)/%.cu.o, $(CUDA_SRC)) $(patsubst src/%.cpp, $(OBJDIR)/%.o, $(COMMON_SRC))
TARGET_CUDA = $(OBJDIR)/mdcuda

# Test source file
TEST_SRC = tests/unittests.cpp
TEST_OBJ = $(patsubst tests/%.cpp, $(OBJDIR)/%.o, $(TEST_SRC))
TARGET_TEST = $(OBJDIR)/unittests

# Link production object files to create the serial executable
$(TARGET_SERIAL): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# Build test executable by linking test and production object files.
# Exclude main.o from the production objects so the Boost.Test main is used.
TEST_PROD_OBJ = $(filter-out $(OBJDIR)/main.o, $(OBJ))
$(TARGET_TEST): $(TEST_OBJ) $(TEST_PROD_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS) -lboost_unit_test_framework

# Build the parallel executable (mdpar) from production sources compiled with OpenMP flags
# and the macro to enable parallel forces computation.
$(TARGET_PAR): $(PAR_OBJ)
	$(CXX) $(PAR_CXXFLAGS) -DPARALLEL_FORCES -o $@ $^ $(LDLIBS)

# Build the CUDA executable (mdcuda) from production sources compiled with NVCC
$(TARGET_CUDA): $(CUDA_OBJ)
	$(NVCC) $(NVCCFLAGS) -DCUDA_ENABLED -o $@ $^ $(LDLIBS)

# Compile production source files into object files (serial)
$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile production source files into parallel object files (with -fopenmp and -DPARALLEL_FORCES)
$(OBJDIR)/%.par.o: src/%.cpp | $(OBJDIR)
	$(CXX) $(PAR_CXXFLAGS) -DPARALLEL_FORCES -c $< -o $@

# Compile CUDA source files into object files
$(OBJDIR)/%.cu.o: src_cuda/%.cu | $(OBJDIR)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Compile .cpp files with NVCC when CUDA is enabled
$(OBJDIR)/%.o: src/%.cpp | $(OBJDIR)
	$(NVCC) $(NVCCFLAGS) -DCUDA_ENABLED -c $< -o $@

# Compile test source files into object files
$(OBJDIR)/%.o: tests/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create the build directory if it doesn't exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Default target builds the serial main executable
all: $(TARGET_SERIAL)

# Target to run tests 
.PHONY: test
test: $(TARGET_TEST)
	./$(TARGET_TEST)

# Target to run the parallel executable.
.PHONY: mdpar
mdpar: $(TARGET_PAR)

# Target to run the CUDA executable.
.PHONY: mdcuda
mdcuda: $(TARGET_CUDA)

# Clean up build files
clean:
	rm -rf $(OBJDIR)

.PHONY: all test clean
