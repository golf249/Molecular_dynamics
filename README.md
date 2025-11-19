# Molecular Dynamics Simulation

A high-performance molecular dynamics simulation program written in C++ that simulates particle interactions using the Lennard-Jones potential. The project supports serial execution, parallel execution with OpenMP, and GPU acceleration with CUDA.

## Features

- **Multiple Execution Modes**:
  - Serial execution (CPU, single-threaded)
  - Parallel execution with OpenMP (CPU, multi-threaded)
  - GPU acceleration with CUDA
  
- **Flexible Initial Conditions**:
  - Predefined test cases for validation
  - Random particle initialization with configurable temperature
  - Support for multiple particle types with different masses

- **Output**:
  - Particle trajectories (positions and velocities)
  - Kinetic energy evolution over time
  - MATLAB plotting scripts for visualization

- **Unit Testing**: Comprehensive test suite using Boost.Test framework

## Prerequisites

### Required Dependencies

- **C++ Compiler**: g++ with C++11 support (g++-10 or newer recommended)
- **Boost Libraries**:
  - `libboost-program-options` for command-line parsing
  - `libboost-test` for unit testing
- **Make**: For building the project

### Optional Dependencies

- **OpenMP**: For parallel execution (usually included with g++)
- **NVIDIA CUDA Toolkit**: For GPU acceleration (nvcc compiler)
- **Doxygen**: For generating documentation
- **MATLAB**: For visualizing simulation results

### Installing Dependencies (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install g++ make libboost-program-options-dev libboost-test-dev
```

For CUDA support, install the NVIDIA CUDA Toolkit from [NVIDIA's website](https://developer.nvidia.com/cuda-downloads).

## Building the Project

The project uses a Makefile for building. Several build targets are available:

### Serial Build (Default)

```bash
make all
```

This creates the executable `build/md` for serial (single-threaded) execution.

### Parallel Build (OpenMP)

```bash
make mdpar
```

This creates the executable `build/mdpar` with OpenMP parallelization.

### CUDA Build (GPU)

```bash
make mdcuda
```

This creates the executable `build/mdcuda` with CUDA GPU acceleration.

### Building All Tests

```bash
make test
```

This compiles and runs the unit tests using Boost.Test framework.

### Cleaning Build Files

```bash
make clean
```

Removes all build artifacts from the `build/` directory.

## Usage

The simulation is controlled via command-line arguments. At minimum, you must specify an initial condition and final simulation time.

### Basic Usage Example

```bash
./build/md --ic-two --T 10.0
```

This runs a simulation with two bouncing particles for 10 time units.

### Command-Line Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--help` | flag | - | Display help message with all options |
| `--Lx` | double | 20.0 | Simulation box length in x-direction (Angstroms) |
| `--Ly` | double | 20.0 | Simulation box length in y-direction (Angstroms) |
| `--Lz` | double | 20.0 | Simulation box length in z-direction (Angstroms) |
| `--dt` | double | 0.001 | Time step for integration |
| `--T` | double | required | Final simulation time |
| `--temp` | double | - | Target temperature (Kelvin) for velocity rescaling |
| `--N` | int | - | Number of particles (required with `--ic-random`) |
| `--percent-type1` | double | 10.0 | Percentage of heavy particles (type 1) |

### Initial Condition Options

You must specify exactly **one** initial condition:

| Option | Description |
|--------|-------------|
| `--ic-one` | One stationary particle |
| `--ic-one-vel` | One moving particle |
| `--ic-two` | Two bouncing particles |
| `--ic-two-pass1` | Two passing particles |
| `--ic-two-pass2` | Two passing particles (close) |
| `--ic-two-pass3` | Two passing particles (close, one heavy) |
| `--ic-random` | N random particles (requires `--N` and optionally `--temp`) |

### Example Commands

**Single moving particle:**
```bash
./build/md --ic-one-vel --T 5.0 --Lx 30 --Ly 30 --Lz 30
```

**Random system with 1000 particles at 80K:**
```bash
./build/md --ic-random --N 1000 --T 1.0 --temp 80 --Lx 50 --Ly 50 --Lz 50
```

**Parallel execution with OpenMP:**
```bash
export OMP_NUM_THREADS=8
./build/mdpar --ic-random --N 10000 --T 0.5 --temp 80 --Lx 50 --Ly 50 --Lz 50
```

**GPU execution with CUDA:**
```bash
./build/mdcuda --ic-random --N 10000 --T 0.5 --temp 80 --Lx 50 --Ly 50 --Lz 50
```

## Output Files

The simulation generates two output files in the `build/` directory:

### `particle_data.txt`

Contains particle trajectories with the following format:
```
time particle_id x y z vx vy vz
```

- `time`: Simulation time
- `particle_id`: Unique particle identifier
- `x, y, z`: Particle position coordinates
- `vx, vy, vz`: Particle velocity components

### `kinetic_energy.txt`

Contains the total kinetic energy of the system over time:
```
time kinetic_energy
```

## Visualization

A MATLAB script `plotting.m` is provided to visualize simulation results:

```bash
matlab -r "run('plotting.m')"
```

The script reads the output files and generates plots of:
- Particle trajectories
- Kinetic energy evolution

## Running Tests

The project includes unit tests to verify simulation accuracy:

```bash
make test
```

This compiles and runs six test cases that validate:
- Single particle motion
- Two-particle interactions
- Energy conservation
- Boundary conditions

## Documentation

The code is documented using Doxygen-style comments. To generate HTML documentation:

```bash
make doc
```

This creates documentation in the `docs/` directory. Open `docs/html/index.html` in a web browser to view.

To clean documentation files:
```bash
make clean-doc
```

## Project Structure

```
Molecular_dynamics/
├── src/                   # Source files
│   ├── main.cpp          # Main entry point
│   ├── md.cpp            # MolecularDynamics class implementation (CPU)
│   ├── md.cu             # CUDA implementation for GPU
│   ├── Particle.cpp      # Particle class implementation
│   └── writeFile.cpp     # File I/O utilities
├── include/              # Header files
│   ├── md.h              # MolecularDynamics class declaration
│   ├── Particle.h        # Particle class declaration
│   └── writeFile.h       # File I/O utilities declaration
├── tests/                # Unit tests
│   └── unittests.cpp     # Boost.Test unit tests
├── plots/                # Output plot directory
├── build/                # Build artifacts (generated)
├── Makefile              # Build configuration
├── Doxyfile              # Doxygen configuration
├── plotting.m            # MATLAB visualization script
├── job-script.slr        # SLURM batch job script example
└── README.md             # This file
```

## High-Performance Computing

For running on HPC clusters with SLURM, an example job script is provided in `job-script.slr`:

```bash
sbatch job-script.slr
```

The script is configured to run the parallel version with OpenMP on 48 cores.

## Physics Model

The simulation uses the **Lennard-Jones potential** to model particle interactions:

```
V(r) = 4ε[(σ/r)^12 - (σ/r)^6]
```

Where:
- `r` is the distance between particles
- `ε` is the depth of the potential well
- `σ` is the finite distance at which the inter-particle potential is zero

The simulation supports two particle types:
- **Type 0**: Light particles (mass = 1.0)
- **Type 1**: Heavy particles (mass = 10.0)

Time integration is performed using the forward Euler method, and periodic boundary conditions are enforced.

## Performance Considerations

- **Serial (`md`)**: Good for small systems (< 100 particles)
- **Parallel (`mdpar`)**: Recommended for medium systems (100-10,000 particles)
- **CUDA (`mdcuda`)**: Best for large systems (> 10,000 particles)

The parallel version uses OpenMP to distribute force calculations across CPU cores. The CUDA version offloads force calculations to the GPU for maximum performance.

## Acknowledgments

The code structure was partially outlined with assistance from ChatGPT 4.0, with the majority of implementation done manually.

## License

Please refer to the repository or contact the author for licensing information.
