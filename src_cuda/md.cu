#include "../include_cuda/md.h"
#include <ctime>
#include <cmath>
#include <omp.h>
#include <iostream>

// New kernel: each thread computes the net force on one particle.
__global__ void computeLJForces(double* position, double* force, int* type, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    
    double force_accum_x = 0.0, force_accum_y = 0.0, force_accum_z = 0.0;
    
    double posx_i = position[3 * i];
    double posy_i = position[3 * i + 1];
    double posz_i = position[3 * i + 2];
    int type_i = type[i];
    
    // Loop over all particles to compute interactions.
    for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        
        double posx_j = position[3 * j];
        double posy_j = position[3 * j + 1];
        double posz_j = position[3 * j + 2];
        int type_j = type[j];
        
        double dx = posx_i - posx_j;
        double dy = posy_i - posy_j;
        double dz = posz_i - posz_j;
        double r2 = dx * dx + dy * dy + dz * dz;
        
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

// I acknowledge the use of ChatGPT 4.0 to generate an outline of how the structure of the MolecularDynamics class should be implemented.
// The majority of the content inside the methods was implemented by me. Where the code that was produced by AI is marked with a comment inside
// the relevant method, 

// Here AI was used to find the way to generate random number using the current time.
MolecularDynamics::MolecularDynamics(int numParticles, double dt, double Lx, double Ly, double Lz, int testCase, 
    double temp, double percent_type1, double finalTime)
    : N(numParticles), dt(dt), Lx(Lx), Ly(Ly), Lz(Lz), testCase(testCase), temp(temp), percent_type1(percent_type1), finalTime(finalTime), 
    writeFile("particle_data.txt", "kinetic_energy.txt"), position_d(nullptr), force_d(nullptr), type_d(nullptr) {
    std::srand(static_cast<unsigned int>(std::time(nullptr))); 
    initialiseParticles();
}

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
}
}


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
    computeLJForces<<<numBlocks, threadsPerBlock>>>(position_d, force_d, type_d, n);
    cudaDeviceSynchronize();
    MolecularDynamics::setParticleForces();
}

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

void MolecularDynamics::setParticleForces() {
    int n = particles.size();
    std::array<double, 3> particleForce;
    double fx, fy, fz;

    for (int i = 0; i < n; ++i) {
        fx = force_d[3 * i];
        fy = force_d[3 * i + 1];
        fz = force_d[3 * i + 2];
        particleForce = {fx, fy, fz};
        particles[i].setForce(particleForce);
    }
}

void MolecularDynamics::bcCheck() {
    // Reflect particles if they go out of bounds
    for (Particle& p : particles) {
        std::array<double, 3> position = p.getPosition();
        std::array<double, 3> velocity = p.getVelocity();
        // For each coordinate: 0->Lx, 1->Ly, 2->Lz
        for (int k = 0; k < 3; ++k) {
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

void MolecularDynamics::runSimulation() {
    // Compute initial kinetic energy
    calKE();
    
    // Output initial conditions before starting integration such that
    // the initial state is also recorded at time 0.
    outputParticleData(0);
    outputKineticEnergy(0);

    const double outputTime = 0.1; // output time for the txt files
    const int steps = static_cast<int>(finalTime / dt); // simulation time steps
    const int outputStep = static_cast<int>(outputTime / dt); // file output time steps

    // Start simulation from step 1
    for (int step = 1; step <= steps; ++step) {
        const double currentTime = step * dt;
        calKE();

        if (testCase == -1) 
            velRescale();
        
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

void MolecularDynamics::outputParticleData(double time) {
    for (size_t i = 0; i < particles.size(); ++i) {
        const Particle& p = particles[i];
        writeFile.writeParticleData(time, i, p.getPosition(), p.getVelocity());
    }
}

void MolecularDynamics::outputKineticEnergy(double time) {
    writeFile.writeKineticEnergy(time, kineticEnergy);
}
