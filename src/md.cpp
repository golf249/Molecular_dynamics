/**
 * @file md.cpp
 * @brief Implementation of the MolecularDynamics class.
 *
 * This file contains the implementation of the MolecularDynamics class, which
 * simulates the behavior of particles in a molecular dynamics system. The class
 * provides methods for initializing particles, calculating forces, updating
 * particle positions and velocities, and running the simulation.
 *
 * The simulation supports various test cases for initializing particles with
 * specific positions and velocities, as well as random initialization. The
 * forces between particles are calculated using the Lennard-Jones potential,
 * and the simulation can be run in parallel using OpenMP.
 */

#include "../include/md.h"
#include <ctime>
#include <cmath>
#include <omp.h>

// I acknowledge the use of ChatGPT 4.0 to generate an outline of how the structure of the MolecularDynamics class should be implemented.
// The majority of the content inside the methods was implemented by me. Where the code that was produced by AI is marked with a comment inside
// the relevant method, 


/**
 * @brief Constructs a MolecularDynamics object with the specified parameters.
 * 
 * @param numParticles The number of particles in the simulation.
 * @param dt The time step for the simulation.
 * @param Lx The length of the simulation box in the x-direction.
 * @param Ly The length of the simulation box in the y-direction.
 * @param Lz The length of the simulation box in the z-direction.
 * @param testCase The test case identifier for the simulation.
 * @param temp The initial temperature of the system.
 * @param percent_type1 The percentage of particles of type 1.
 * @param finalTime The final time for the simulation.
 */
MolecularDynamics::MolecularDynamics(int numParticles, double dt, double Lx, double Ly, double Lz, int testCase, double temp, double percent_type1, double finalTime)
    : N(numParticles), dt(dt), Lx(Lx), Ly(Ly), Lz(Lz), testCase(testCase), temp(temp), percent_type1(percent_type1), finalTime(finalTime), writeFile("particle_data.txt", "kinetic_energy.txt") {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    initialiseParticles();
}

/**
 * @brief Destructor for the MolecularDynamics class.
 *
 * This destructor is responsible for cleaning up any resources
 * allocated by the MolecularDynamics class. Nothing to clean up in this case.
 */
MolecularDynamics::~MolecularDynamics() {
}

/**
 * @brief Initialises the position, velocity and mass of particles for the molecular dynamics simulation based on the specified test case.
 * 
 * This function first clears any existing particles and initialises new particles based on the value of the `testCase` member variable.
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
void MolecularDynamics::initialiseParticles() {
    particles.clear(); // Clear any existing particles
    // Here AI was used to generate the code to set the positions and velocity since 
    // it is a repetitive task.
    if (testCase == 1) {
        particles = {
            Particle({10.0, 10.0, 10.0}, {0.0, 0.0, 0.0}, 0)
        };
    } else if (testCase == 2) {
        particles = {
            Particle({10.0, 10.0, 10.0}, {5.0, 2.0, 1.0}, 0)
        };
    } else if (testCase == 3) {
        particles = {
            Particle({8.5, 10.0, 10.0}, {0.0, 0.0, 0.0}, 0),
            Particle({11.5, 10.0, 10.0}, {0.0, 0.0, 0.0}, 0)
        };
    } else if (testCase == 4) {
        particles = {
            Particle({8.5, 11.5, 10.0}, {0.5, 0.0, 0.0}, 0),
            Particle({11.5, 8.5, 10.0}, {-0.5, 0.0, 0.0}, 0)
        };
    } else if (testCase == 5) {
        particles = {
            Particle({8.5, 11.3, 10.0}, {0.5, 0.0, 0.0}, 0),
            Particle({11.5, 8.7, 10.0}, {-0.5, 0.0, 0.0}, 0)
        };
    } else if (testCase == 6) {
        particles = {
            Particle({8.5, 11.3, 10.0}, {0.5, 0.0, 0.0}, 1),
            Particle({11.5, 8.7, 10.0}, {-0.5, 0.0, 0.0}, 1)
        };
    } else {
        int numType1 = static_cast<int>(std::ceil(percent_type1 / 100.0 * N));
        int numType0 = N - numType1;
        
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
    int count = 0;

    for (int i = 0; i < numType; ++i) {
        bool validPosition = false;
        std::array<double, 3> position, velocity;

        while (!validPosition) {
            position = { Lx * ((double)rand() * invRand), Ly * ((double)rand() * invRand), Lz * ((double)rand() * invRand) };
            velocity = { ((double)rand() * invRand - 0.5), ((double)rand() * invRand - 0.5), ((double)rand() * invRand - 0.5) };

            validPosition = stabilityCheck(position);
            if (count > 1000) {
                throw std::runtime_error("Failed to generate stable initial configuration.");
            }
            count++;
        }

        particles.emplace_back(position, velocity, type);
    }
}

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
    // Here AI was used to outline implementation of how to detect whether particles are too close 
    // to each other, specifically the norm of the difference between the positions of two particles.
    constexpr double R2 = 0.25; // Minimum allowed distance squared (0.5^2) 
    for (const Particle& p : particles) {
        const std::array<double, 3> otherPosition = p.getPosition();
        double dx = position[0] - otherPosition[0];
        double dy = position[1] - otherPosition[1];
        double dz = position[2] - otherPosition[2];
        if (dx*dx + dy*dy + dz*dz < R2) {
            return false;
        }
    }
    return true;
}

// Here AI was used to find the best way to optimize the performance of the forces calculation between 
// particles, specifically using a lookup table to store the interaction parameters, the use of 'static'
// and 'constexpr', and the use of sigma2, sigma6, sigma12, inv_r2, inv_r6, inv_r12 to reduce the number 
// of calculations and floating point division for a better performance overall.

/**
 * @brief Calculate the forces acting on each particle in the system.
 *
 * This function computes the pairwise forces between particles using the 
 * Lennard-Jones potential. It first clears the existing forces on each 
 * particle, then iterates over all unique pairs of particles to compute 
 * the interaction forces. The forces are updated based on the interaction 
 * parameters (sigma and epsilon) which are looked up from predefined tables.
 *
 * The Lennard-Jones potential is given by:
 *     V(r) = 4 * epsilon * [(sigma / r)^12 - (sigma / r)^6]
 * where r is the distance between particles.
 *
 * The force is derived from the potential as:
 *     F(r) = 24 * epsilon * (2 * (sigma / r)^12 - (sigma / r)^6) / r^2
 *
 * @note This function assumes that the particles are stored in a vector 
 *       and that each particle has methods to get and set its position, 
 *       type, and force.
 */
void MolecularDynamics::calForces() {
    // Clear forces first
    for (Particle& p : particles) {
        p.setForce({0.0, 0.0, 0.0});
    }

    // Use lookup tables for interaction parameters
    static constexpr double sigma2_table[2][2] = {
        {1.0, 4.0},
        {4.0, 9.0}
    };
    static constexpr double epsilon_table[2][2] = {
        {3.0, 15.0},
        {15.0, 60.0}
    };

    const int n = particles.size();
    for (int i = 0; i < n; ++i) {
        // Cache particle i position and type to reduce repeated lookups
        const auto& pos_i = particles[i].getPosition();
        const int type_i = particles[i].getType();
        auto& force_i = particles[i].getForce();

        for (int j = i + 1; j < n; ++j) {
            // Cache particle j position and type
            const auto& pos_j = particles[j].getPosition();
            const int type_j = particles[j].getType();
            auto& force_j = particles[j].getForce();

            // Compute the vector between particles and squared distance
            const double dx = pos_i[0] - pos_j[0];
            const double dy = pos_i[1] - pos_j[1];
            const double dz = pos_i[2] - pos_j[2];
            const double r2 = dx*dx + dy*dy + dz*dz;

            // Lookup parameters for sigma and epsilon
            const double sigma2 = sigma2_table[type_i][type_j];
            const double epsilon = epsilon_table[type_i][type_j];

            // Here AI was used to find the most efficient way to represent the variables that is to be
            // used in the force calculation, specifically the use of sigma6, sigma12, inv_r2, inv_r6, inv_r12
            // Pre-compute powers and inverses
            const double sigma6 = sigma2 * sigma2 * sigma2;
            const double sigma12 = sigma6 * sigma6;
            const double inv_r2 = 1.0 / r2;
            const double inv_r6 = inv_r2 * inv_r2 * inv_r2;
            const double inv_r12 = inv_r6 * inv_r6;

            // Compute magnitude of the force
            const double f = 24.0 * epsilon * inv_r2 * (2.0 * sigma12 * inv_r12 - sigma6 * inv_r6);

            // Update forces along each component
            force_i[0] += f * dx;
            force_i[1] += f * dy;
            force_i[2] += f * dz;
            force_j[0] -= f * dx;
            force_j[1] -= f * dy;
            force_j[2] -= f * dz;
        }
    }
};

/**
 * @brief Calculate forces on particles in parallel using OpenMP.
 *
 * This function calculates the forces on particles in a molecular dynamics
 * simulation using a parallel approach with OpenMP. It first clears the forces
 * on all particles, then computes the forces in parallel, and finally combines
 * the thread-local force accumulations into the particles' forces.
 *
 * The force calculation uses a Lennard-Jones potential with precomputed
 * interaction parameters for different particle types.
 *
 * @note If OpenMP is not enabled, the function will run in a single thread.
 */
void MolecularDynamics::calForcesParallel() {
    const int n = particles.size();

    // Clear forces globally.
    for (int i = 0; i < n; ++i) {
        particles[i].setForce({0.0, 0.0, 0.0});
    }

    // Determine number of threads. If OpenMP is enabled, use the maximum number of threads.
    int num_threads = 1;
    #ifdef _OPENMP
        num_threads = omp_get_max_threads();
    #endif
    
    // thread_forces[thread_id][i] will store local force for particle i computed by thread threadID in a 3D array.
    std::vector<std::vector<std::array<double, 3>>> thread_forces(num_threads, 
        std::vector<std::array<double, 3>>(n, {0.0, 0.0, 0.0}));

    // Use lookup tables for interaction parameters.
    static constexpr double sigma2_table[2][2] = {
        {1.0, 4.0},
        {4.0, 9.0}
    };
    static constexpr double epsilon_table[2][2] = {
        {3.0, 15.0},
        {15.0, 60.0}
    };

    #ifdef _OPENMP
        #pragma omp parallel
    #endif
    { // Start parallel region
        int thread_id = 0;

        #ifdef _OPENMP
            thread_id = omp_get_thread_num();
        #endif
        #ifdef _OPENMP
            #pragma omp for schedule(dynamic)
        #endif
        for (int i = 0; i < n; ++i) {
            const auto& pos_i = particles[i].getPosition();
            const int type_i = particles[i].getType();

            // Loop from i+1 to take advantage of symmetry in the force calculation.
            // This also reduces the number of calculations by half.
            for (int j = i + 1; j < n; ++j) {
                const auto& pos_j = particles[j].getPosition();
                const int type_j = particles[j].getType();
                
                // Compute the vector between particles and squared distance
                const double dx = pos_i[0] - pos_j[0];
                const double dy = pos_i[1] - pos_j[1];
                const double dz = pos_i[2] - pos_j[2];
                const double r2 = dx * dx + dy * dy + dz * dz;
                
                // Lookup parameters for sigma and epsilon
                const double sigma2 = sigma2_table[type_i][type_j];
                const double epsilon = epsilon_table[type_i][type_j];
                
                // Pre-compute powers and inverses
                const double sigma6 = sigma2 * sigma2 * sigma2;
                const double sigma12 = sigma6 * sigma6;
                const double inv_r2 = 1.0 / r2;
                const double inv_r6 = inv_r2 * inv_r2 * inv_r2;
                const double inv_r12 = inv_r6 * inv_r6;
                
                // Compute magnitude of force
                const double f = 24.0 * epsilon * inv_r2 * (2.0 * sigma12 * inv_r12 - sigma6 * inv_r6);
                
                // Accumulate contributions into the thread-local arrays.
                thread_forces[thread_id][i][0] += f * dx;
                thread_forces[thread_id][i][1] += f * dy;
                thread_forces[thread_id][i][2] += f * dz;
                
                thread_forces[thread_id][j][0] -= f * dx;
                thread_forces[thread_id][j][1] -= f * dy;
                thread_forces[thread_id][j][2] -= f * dz;
            }
        }
    } // End parallel region

    // Combine the thread-local force accumulations into the particles’ force
    #ifdef _OPENMP
        #pragma omp parallel for
    #endif
    for (int i = 0; i < n; ++i) {
        std::array<double, 3> accum_force = {0.0, 0.0, 0.0};
        for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
            accum_force[0] += thread_forces[thread_id][i][0];
            accum_force[1] += thread_forces[thread_id][i][1];
            accum_force[2] += thread_forces[thread_id][i][2];
        }
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
 * the forward Euler method and computing forces either in parallel or serially based
 * on the preprocessor macro. The particle data and kinetic energy are output at
 * specified time intervals.
 * 
 * @details
 * - Computes initial kinetic energy and outputs initial conditions.
 * - Iterates over the simulation time steps:
 *   - Updates the current time.
 *   - Computes kinetic energy.
 *   - Computes forces using either parallel or serial implementation.
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

        // Use a preprocessor macro to determine which implementation of calForces to use (serial/ OpenMP)
        #ifdef PARALLEL_FORCES
            calForcesParallel();
        #else
            calForces();
        #endif
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
    int n = particles.size();
    for (int i = 0; i < n; ++i) {
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