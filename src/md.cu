/**
 * @file md.cu
 * @brief Implementation of the MolecularDynamics class.
 * 
 * This file contains the implementation of the MolecularDynamics class, which
 * is responsible for simulating molecular dynamics. The class includes methods
 * for initializing particles, computing forces, updating positions and velocities,
 * and handling boundary conditions.
 * 
 * The implementation leverages CUDA for parallel computation to efficiently
 * handle large numbers of particles and interactions.
 * 
 * @note Ensure that the CUDA runtime and necessary libraries are properly
 * installed and configured before compiling and running this code.
 */

#include "../include/md.h"
#include <ctime>
#include <cmath>
#include <omp.h>
#include <iostream>

// I acknowledge the use of ChatGPT 4.0 to generate an outline of how the structure of the MolecularDynamics class should be implemented.
// The majority of the content inside the methods was implemented by me. Where the code that was produced by AI is marked with a comment inside
// the relevant method, 

/**
 * @brief Compute Lennard-Jones forces for a set of particles.
 *
 * This kernel computes the Lennard-Jones forces acting on each particle
 * based on their positions and types. The forces are stored in the provided
 * force array.
 *
 * @param position Pointer to the array of particle positions. The array should
 *                 be of size 3 * n, where n is the number of particles.
 * @param force    Pointer to the array where the computed forces will be stored.
 *                 The array should be of size 3 * n.
 * @param type     Pointer to the array of particle types. The array should be
 *                 of size n.
 * @param n        The number of particles.
 */
__global__ void calLJForces(double* position, double* force, int* type, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    
    double force_accum_x = 0.0, force_accum_y = 0.0, force_accum_z = 0.0;
    
    // Cache particle i position and type to reduce repeated lookups
    double posx_i = position[3 * i];
    double posy_i = position[3 * i + 1];
    double posz_i = position[3 * i + 2];
    int type_i = type[i];
    
    // Loop over all particles to compute interactions.
    for (int j = 0; j < n; ++j) {
        if (i == j) continue;

        // Cache particle j position and type to reduce repeated lookups
        double posx_j = position[3 * j];
        double posy_j = position[3 * j + 1];
        double posz_j = position[3 * j + 2];
        int type_j = type[j];
        
        // Compute the vector between particles and squared distance
        double dx = posx_i - posx_j;
        double dy = posy_i - posy_j;
        double dz = posz_i - posz_j;
        double r2 = dx * dx + dy * dy + dz * dz;
        
        // Find correct parameters for sigma and epsilon based on particle types
        double sigma2, epsilon;
        if (type_i == 0 && type_j == 0) {
            sigma2 = 1.0;
            epsilon = 3.0;
        } else if (type_i == 1 && type_j == 1) {
            sigma2 = 9.0;
            epsilon = 60.0;
        } else {
            sigma2 = 4.0;
            epsilon = 15.0;
        }
        
        // Here AI was used to find the most efficient way to represent the variables that is to be
            // used in the force calculation, specifically the use of sigma6, sigma12, inv_r2, inv_r6, inv_r12
            // Pre-compute powers and inverses
        double sigma6 = sigma2 * sigma2 * sigma2;
        double sigma12 = sigma6 * sigma6;
        double inv_r2 = 1.0 / r2;
        double inv_r6 = inv_r2 * inv_r2 * inv_r2;
        double inv_r12 = inv_r6 * inv_r6;
        
        double f = 24.0 * epsilon * inv_r2 * (2.0 * sigma12 * inv_r12 - sigma6 * inv_r6);
        
        force_accum_x += f * dx;
        force_accum_y += f * dy;
        force_accum_z += f * dz;

        //printf("i,j %d,%d Force: (%f, %f, %f)\n", i, j, f * dx, dy, dz);
    }
    
    // Write the net force for particle i to global memory.
    force[3 * i]     = force_accum_x;
    force[3 * i + 1] = force_accum_y;
    force[3 * i + 2] = force_accum_z;
}


/**
 * @brief Constructs a MolecularDynamics object with the specified parameters.
 * 
 * @param numParticles The number of particles in the simulation.
 * @param dt The time step for the simulation.
 * @param Lx The length of the simulation box in the x-dimension.
 * @param Ly The length of the simulation box in the y-dimension.
 * @param Lz The length of the simulation box in the z-dimension.
 * @param testCase The test case identifier for initializing the simulation.
 * @param temp The initial temperature of the system.
 * @param percent_type1 The percentage of particles of type 1.
 * @param finalTime The final time for the simulation.
 */
MolecularDynamics::MolecularDynamics(int numParticles, double dt, double Lx, double Ly, double Lz, int testCase, 
    double temp, double percent_type1, double finalTime)
    : N(numParticles), dt(dt), Lx(Lx), Ly(Ly), Lz(Lz), testCase(testCase), temp(temp), percent_type1(percent_type1), finalTime(finalTime), 
    writeFile("particle_data.txt", "kinetic_energy.txt"), position_d(nullptr), force_d(nullptr), type_d(nullptr) {
    std::srand(static_cast<unsigned int>(std::time(nullptr))); 
    initialiseParticles();
}

/**
 * @brief Destructor for the MolecularDynamics class.
 *
 * This destructor is responsible for cleaning up any resources
 * allocated to the GPU. It ensures that
 * all memory are properly released when
 * an instance of the class is destroyed.
 */
MolecularDynamics::~MolecularDynamics() {
    // Free GPU memory
    cudaFree(position_d);
    cudaFree(force_d);
    cudaFree(type_d);
}

/**
 * @brief Initialises the position, velocity and mass of particles for the molecular dynamics simulation based on the specified test case.
 * 
 * This function first clears any existing particles and initialises new particles based on the value of the `testCase` member variable.
 * For each test case, the memory for the position, force and type arrays are allocated on the GPU using `cudaMallocManaged`.
 * 
 * Test cases:
 * - 1: Initialies a single particle at position (10.0, 10.0, 10.0) with zero velocity.
 * - 2: Initialises a single particle at position (10.0, 10.0, 10.0) with velocity (5.0, 2.0, 1.0).
 * - 3: Initialises two particles at positions (8.5, 10.0, 10.0) and (11.5, 10.0, 10.0) with zero velocity.
 * - 4: Initialises two particles at positions (8.5, 11.5, 10.0) and (11.5, 8.5, 10.0) with velocities (0.5, 0.0, 0.0) and (-0.5, 0.0, 0.0) respectively.
 * - 5: Initialises two particles at positions (8.5, 11.3, 10.0) and (11.5, 8.7, 10.0) with velocities (0.5, 0.0, 0.0) and (-0.5, 0.0, 0.0) respectively.
 * - 6: Initialises two particles at positions (8.5, 11.3, 10.0) and (11.5, 8.7, 10.0) with velocities (0.5, 0.0, 0.0) and (-0.5, 0.0, 0.0) respectively, and type 1 mass.
 * - Random: Initialises a specified number of particles of type 0 and type 1 with random positions and velocities, ensuring the positions between
 *  particles isn't <= 5 unit length. This is done via the `stabilityCheck` function.
 * 
 * The number of particles and its type is determined by the `percent_type1` and `N` member variables.
 * 
 * @note The random positions and velocities are generated using the `rand()` function and to ensure true randomness at every generation, a seed using the current time was used.
 */
// Here AI was used to generate the code to set the positions and velocity since 
// it is a repetitive task, but the structure of the conditional statement was implemented by me.
void MolecularDynamics::initialiseParticles() {
    particles.clear(); // Clear any existing particles

    if (testCase == 1) {
        particles = {
            Particle({10.0, 10.0, 10.0}, {0.0, 0.0, 0.0}, 0)
        };
        int N = particles.size();
        // Allocate GPU memory on the device
        cudaMallocManaged(&position_d, 3 * N * sizeof(double));
        cudaMallocManaged(&force_d, 3 * N * sizeof(double));
        cudaMallocManaged(&type_d, N * sizeof(int));
    } else if (testCase == 2) {
        particles = {
            Particle({10.0, 10.0, 10.0}, {5.0, 2.0, 1.0}, 0)
        };
        int N = particles.size();
        // Allocate GPU memory on the device
        cudaMallocManaged(&position_d, 3 * N * sizeof(double));
        cudaMallocManaged(&force_d, 3 * N * sizeof(double));
        cudaMallocManaged(&type_d, N * sizeof(int));
    } else if (testCase == 3) {
        particles = {
            Particle({8.5, 10.0, 10.0}, {0.0, 0.0, 0.0}, 0),
            Particle({11.5, 10.0, 10.0}, {0.0, 0.0, 0.0}, 0)
        };
        int N = particles.size();
        // Allocate GPU memory on the device
        cudaMallocManaged(&position_d, 3 * N * sizeof(double));
        cudaMallocManaged(&force_d, 3 * N * sizeof(double));
        cudaMallocManaged(&type_d, N * sizeof(int));
    } else if (testCase == 4) {
        particles = {
            Particle({8.5, 11.5, 10.0}, {0.5, 0.0, 0.0}, 0),
            Particle({11.5, 8.5, 10.0}, {-0.5, 0.0, 0.0}, 0)
        };
        int N = particles.size();
        // Allocate GPU memory on the device
        cudaMallocManaged(&position_d, 3 * N * sizeof(double));
        cudaMallocManaged(&force_d, 3 * N * sizeof(double));
        cudaMallocManaged(&type_d, N * sizeof(int));
    } else if (testCase == 5) {
        particles = {
            Particle({8.5, 11.3, 10.0}, {0.5, 0.0, 0.0}, 0),
            Particle({11.5, 8.7, 10.0}, {-0.5, 0.0, 0.0}, 0)
        };
        int N = particles.size();
        // Allocate GPU memory on the device
        cudaMallocManaged(&position_d, 3 * N * sizeof(double));
        cudaMallocManaged(&force_d, 3 * N * sizeof(double));
        cudaMallocManaged(&type_d, N * sizeof(int));
    } else if (testCase == 6) {
        particles = {
            Particle({8.5, 11.3, 10.0}, {0.5, 0.0, 0.0}, 1),
            Particle({11.5, 8.7, 10.0}, {-0.5, 0.0, 0.0}, 1)
        };
        int N = particles.size();
        // Allocate GPU memory on the device
        cudaMallocManaged(&position_d, 3 * N * sizeof(double));
        cudaMallocManaged(&force_d, 3 * N * sizeof(double));
        cudaMallocManaged(&type_d, N * sizeof(int));
    } else {
        int numType1 = static_cast<int>(std::ceil(percent_type1 / 100.0 * N));
        int numType0 = N - numType1;
        // Allocate GPU memory on the device
        cudaMallocManaged(&position_d, 3 * N * sizeof(double));
        cudaMallocManaged(&force_d, 3 * N * sizeof(double));
        cudaMallocManaged(&type_d, N * sizeof(int));

        assignRandomStates(numType1, 1);
        assignRandomStates(numType0, 0);
}
}

/**
 * @brief Assigns random positions and velocities to a specified number of particles of a given type.
 *
 * This function generates random positions and velocities for a specified number of particles
 * of a given type. The positions are generated within the bounds of the simulation box defined
 * by Lx, Ly, and Lz. The velocities are generated randomly with components in the range [-0.5, 0.5].
 * The function ensures that the generated positions are valid by checking them with the stabilityCheck function.
 *
 * @param numType The number of particles to assign random states to.
 * @param type The type of particles to assign random states to.
 */
void MolecularDynamics::assignRandomStates(const int numType, const int type) {
    const double invRand = 1.0 / RAND_MAX;
    for (int i = 0; i < numType; ++i) {
        bool validPosition = false;
        std::array<double, 3> position, velocity;

        while (!validPosition) {
            position = { Lx * ((double)rand() * invRand), Ly * ((double)rand() * invRand), Lz * ((double)rand() * invRand) };
            velocity = { ((double)rand() * invRand - 0.5), ((double)rand() * invRand - 0.5), ((double)rand() * invRand - 0.5) };

            validPosition = stabilityCheck(position);
        }

        particles.emplace_back(position, velocity, type);
    }
}

// Here AI was used to outline implementation of how to detect whether particles are too close 
// to each other, specifically the norm of the difference between the positions of two particles.
/**
 * @brief Checks the stability of a given position within the molecular dynamics system.
 *
 * This function determines whether a given position is stable by ensuring that it is not
 * too close to any existing particles in the system. The minimum allowed distance between
 * particles is defined by the constant R2 (0.25, which is 0.5 squared).
 *
 * @param position The position to check, represented as an array of three doubles.
 * @return true if the position is stable (i.e., not too close to any other particle), false otherwise.
 */
bool MolecularDynamics::stabilityCheck(const std::array<double, 3>& position) {
    constexpr double R2 = 0.25; // Minimum allowed distance squared (0.5^2) {
    for (const Particle& p : particles) {
        const std::array<double, 3> otherPosition = p.getPosition();
        double dx = position[0] - otherPosition[0];
        double dy = position[1] - otherPosition[1];
        double dz = position[2] - otherPosition[2];
        if ((dx*dx + dy*dy + dz*dz) < R2) {
            return false;
        }
    }
    return true;
}


/**
 * @brief Calculates the forces acting on particles in the molecular dynamics simulation using CUDA.
 * 
 * This function utilizes CUDA to perform parallel computations of the forces between particles.
 * First, the function copies the current particle positions and types into unified memory arrays,
 * resets the force array to zero, and then launches the CUDA kernel "calLJForces". In this kernel,
 * each thread is assigned to a single particle and computes its net force by iterating over all
 * other particles using the Lennard-Jones potential. Predefined lookup tables are used for the
 * interaction parameters. After synchronizing the device, the computed forces are copied back to
 * the particle data structure via setForce().
 */
void MolecularDynamics::calForcesCUDA() {
    const int n = particles.size();
    
    // Populate managed arrays from particle data.
    for (int i = 0; i < n; i++) {
        const std::array<double, 3>& pos = particles[i].getPosition();

        // Get the position of the particle
        position_d[3 * i] = pos[0];
        position_d[3 * i + 1] = pos[1];
        position_d[3 * i + 2] = pos[2];

        // Set the forces for each particle to 0
        force_d[3 * i] = 0.0;
        force_d[3 * i + 1] = 0.0;
        force_d[3 * i + 2] = 0.0;

        // Get particle type.
        type_d[i] = particles[i].getType();
    }
    
    // Launch the kernel with one thread per particle.
    int threadsPerBlock = 256;
    int numBlocks = (n + threadsPerBlock - 1) / threadsPerBlock;
    calLJForces<<<numBlocks, threadsPerBlock>>>(position_d, force_d, type_d, n);
    cudaDeviceSynchronize();

    for (int i = 0; i < n; ++i) {
        std::array<double, 3> accum_force = {0.0,0.0,0.0};
        accum_force[0] = force_d[3 * i];
        accum_force[1] = force_d[3 * i + 1];
        accum_force[2] = force_d[3 * i + 2];
        particles[i].setForce(accum_force);
    }    
}

/**
 * @brief Advances the state of the molecular dynamics simulation using the forward Euler method.
 *
 * This method updates the velocity and position of each particle in the simulation based on the
 * current forces acting on them. The forward Euler method is a simple numerical integration technique
 * that approximates the new state of the system over a small time step (dt).
 *
 * The velocity of each particle is updated using the formula:
 *     v_new = v_old + (dt * force / mass)
 *
 * The position of each particle is then updated using the formula:
 *     x_new = x_old + (dt * v_new)
 *
 * After updating the velocities and positions, the method checks for boundary conditions.
 *
 * @note This method assumes that the particles have already been initialized with their respective
 *       velocities, positions, forces, and masses.
 */
void MolecularDynamics::forwardEuler() {
    for (Particle& p : particles) {
        std::array<double, 3> velocity = p.getVelocity();
        std::array<double, 3> position = p.getPosition();
        const std::array<double, 3> force = p.getForce();
        const double mass = p.getMass();
        for (int k = 0; k < 3; ++k) {
            velocity[k] += dt * force[k] / mass;
            position[k] += dt * velocity[k];
        }
        p.setVelocity(velocity);
        p.setPosition(position);
    }
    bcCheck();
}

/**
 * @brief Checks and applies boundary conditions to particles in the simulation.
 *
 * This function iterates over all particles in the simulation and ensures that
 * they remain within the defined simulation box boundaries. If a particle goes
 * out of bounds, it is reflected back into the simulation box, and its velocity
 * is adjusted accordingly to simulate a reflective boundary condition.
 *
 * The simulation box is defined by the dimensions Lx, Ly, and Lz along the x, y,
 * and z axes, respectively. For each particle, the function checks its position
 * along each axis:
 * - If the position is less than 0, the particle is reflected back into the box
 *   by setting its position to the negative of its current position and its
 *   velocity to the absolute value of its current velocity.
 * - If the position is greater than the box dimension, the particle is reflected
 *   back into the box by setting its position to twice the box dimension minus
 *   its current position and its velocity to the negative absolute value of its
 *   current velocity.
 *
 * After adjusting the position and velocity of a particle, the function updates
 * the particle's state.
 */
void MolecularDynamics::bcCheck() {
    // Reflect particles if they go out of bounds
    for (Particle& p : particles) {
        std::array<double, 3> position = p.getPosition();
        std::array<double, 3> velocity = p.getVelocity();
        // For each coordinate: 0->Lx, 1->Ly, 2->Lz
        for (int k = 0; k < 3; ++k) {
            // Here AI was used only to help define L in a very compact and efficient form.
            double L = (k == 0) ? Lx : ((k == 1) ? Ly : Lz); // Select the correct bound size for each axis
            if (position[k] < 0) {
                position[k] = -position[k];
                velocity[k] = std::abs(velocity[k]);
            } else if (position[k] > L) {
                position[k] = 2 * L - position[k];
                velocity[k] = -std::abs(velocity[k]);
            }
        }
        p.setPosition(position);
        p.setVelocity(velocity);
    }
}

/**
 * @brief Calculates the total kinetic energy of the system.
 *
 * This function iterates over all particles in the system at current 
 * time, retrieves their velocities, and computes the kinetic energy 
 * using the formula:
 * 
 *     KE = 0.5 * mass * (velocity_x^2 + velocity_y^2 + velocity_z^2)
 * 
 * The total kinetic energy is then accumulated and stored in the member
 * variable `kineticEnergy`.
 */
void MolecularDynamics::calKE() {
    kineticEnergy = 0.0;
    for (const Particle& p : particles) {
        const std::array<double, 3>& velocity = p.getVelocity();
        double speedSquared = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]; 
        kineticEnergy += 0.5 * p.getMass() * speedSquared;
    }
}

/**
 * @brief Rescales the velocities of particles to match the desired temperature.
 *
 * This function rescales the velocities of all particles in the system to ensure
 * that the kinetic temperature matches the desired temperature (`temp`). If the
 * temperature (`temp`) is not set (i.e., it is -1.0), the function returns immediately.
 *
 * The rescaling factor (`lambda`) is calculated based on the ratio of the desired
 * temperature to the current kinetic temperature. Each component of the velocity
 * of every particle is then multiplied by this factor.
 *
 * @note The Boltzmann constant (`kb`) is defined as 0.8314459920816467.
 */
void MolecularDynamics::velRescale() {
    // check if temp is not set else continue
    if (temp == -1.0) {
        return;
    }

    constexpr double kb = 0.8314459920816467; // Boltzmann constant
    double tempKE = (2.0 / (3.0 * kb * N)) * kineticEnergy;
    const double lambda = std::sqrt(temp / tempKE);

    for (Particle& p : particles) {
        std::array<double, 3> velocity = p.getVelocity();
        for (int k = 0; k < 3; ++k) {
            velocity[k] *= lambda;
        }
        p.setVelocity(velocity);
    }
}

/**
 * @brief Runs the molecular dynamics simulation.
 * 
 * This function initializes the kinetic energy and outputs the initial conditions.
 * It then iterates over the simulation time steps, updating the system state using
 * the forward Euler method and computing forces using CUDA. The particle data and 
 * kinetic energy are output at specified time intervals.
 * 
 * @details
 * - Computes initial kinetic energy and outputs initial conditions.
 * - Iterates over the simulation time steps:
 *   - Updates the current time.
 *   - Computes kinetic energy.
 *   - Computes forces using CUDA.
 *   - Updates particle positions and velocities using forward Euler method.
 *   - Outputs particle data and kinetic energy at specified intervals.
 * 
 * @note If `testCase` is -1, velocity rescaling is performed.
 * 
 * @pre The simulation parameters such as `finalTime`, `dt`, and `testCase` must be set.
 * @pre The output functions `outputParticleData` and `outputKineticEnergy` must be defined.
 * 
 * @param None
 * @return None
 */
void MolecularDynamics::runSimulation() {
    // Compute initial kinetic energy
    calKE();
    
    // Output initial conditions before starting integration such that
    // the initial state is also recorded at time 0.
    if (testCase != -1) {
        outputParticleData(0);
    }
    outputKineticEnergy(0);

    const double outputTime = 0.1; // output time for the txt files
    const int steps = static_cast<int>(finalTime / dt); // simulation time steps
    const int outputStep = static_cast<int>(outputTime / dt); // file output time steps

    // Start simulation from step 1
    for (int step = 1; step <= steps; ++step) {
        const double currentTime = step * dt;
        calKE();

        if (testCase == -1) {
            velRescale();
        }
        
        calForcesCUDA();

        forwardEuler();

        // Output particle data and kinetic energy at 0.1 unit time intervals for the output files
        if (step % outputStep == 0) {
            if (testCase != -1) {
                outputParticleData(currentTime);
            }
            outputKineticEnergy(currentTime);
        }
    }
}

/**
 * @brief Outputs the data of all particles at a given time.
 *
 * This function iterates over all particles in the system and writes their
 * data (position and velocity) to a file.
 *
 * @param time The current simulation time.
 */
void MolecularDynamics::outputParticleData(double time) {
    for (int i = 0; i < particles.size(); ++i) {
        const Particle& p = particles[i];
        writeFile.writeParticleData(time, i, p.getPosition(), p.getVelocity());
    }
}


/**
 * @brief Outputs the kinetic energy of the system at a given time.
 * 
 * This function writes the current kinetic energy of the molecular dynamics
 * system to a file, associating it with the specified time.
 * 
 * @param time The current time at which the kinetic energy is being recorded.
 */
void MolecularDynamics::outputKineticEnergy(double time) {
    writeFile.writeKineticEnergy(time, kineticEnergy);
}
