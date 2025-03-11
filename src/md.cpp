#include "../include/md.h"
#include <ctime>
#include <cmath>

MolecularDynamics::MolecularDynamics(int numParticles, double dt, double Lx, double Ly, double Lz, int testCase, double temp, double percent_type1, double finalTime)
    : N(numParticles), dt(dt), Lx(Lx), Ly(Ly), Lz(Lz), testCase(testCase), temp(temp), percent_type1(percent_type1), finalTime(finalTime), writeFile("particle_data.txt", "kinetic_energy.txt") {
    std::srand(static_cast<unsigned int>(std::time(nullptr))); // Seed the random number generator with the current time
    initializeParticles();

    // Set final time based on test case if not provided
    if (this->finalTime == -1.0) {
        if (testCase == 1) {
            this->finalTime = 1.0;
        } else if (testCase == 2) {
            this->finalTime = 20.0;
        } else if (testCase >= 3 && testCase <= 6) {
            this->finalTime = 50.0;
        }
    }
}

void MolecularDynamics::initializeParticles() {
    particles.clear(); // Clear any existing particles

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
        const double invRand = 1.0 / RAND_MAX;

        for (int i = 0; i < numType0; ++i) {
            bool validPosition = false;
            std::array<double, 3> position, velocity;

            while (!validPosition) {
                position = { Lx * ((double)rand() * invRand), Ly * ((double)rand() * invRand), Lz * ((double)rand() * invRand) };
                velocity = { ((double)rand() * invRand - 0.5), ((double)rand() * invRand - 0.5), ((double)rand() * invRand - 0.5) };

                validPosition = stabilityCheck(position);
            }

            particles.emplace_back(position, velocity, 0); // append a new type 0 particle to the particles vector with the given position, velocity, and type
        }

        for (int i = 0; i < numType1; ++i) {
            bool validPosition = false;
            std::array<double, 3> position, velocity;

            while (!validPosition) {
                position = {Lx * (double)(rand()) / RAND_MAX, Ly * (double)(rand()) / RAND_MAX, Lz * (double)(rand()) / RAND_MAX};
                velocity = {(double)(rand()) / RAND_MAX - 0.5, (double)(rand()) / RAND_MAX - 0.5, (double)(rand()) / RAND_MAX - 0.5};

                validPosition = stabilityCheck(position);
            }

            particles.emplace_back(position, velocity, 1); // append a new type 1 particle to the particles vector with the given position, velocity, and type
        }
    }
}

bool MolecularDynamics::stabilityCheck(const std::array<double, 3>& position) {
    constexpr double minDist2 = 0.25; // Minimum allowed distance squared (0.5^2) {
    for (const Particle& p : particles) {
        const std::array<double, 3> otherPosition = p.getPosition();
        double dx = position[0] - otherPosition[0];
        double dy = position[1] - otherPosition[1];
        double dz = position[2] - otherPosition[2];
        if ((dx*dx + dy*dy + dz*dz) < minDist2) {
            return false;
        }
    }
    return true;
}

void MolecularDynamics::calForces() {
    // Clear forces first
    for (auto& p : particles) {
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

    const size_t n = particles.size();
    for (size_t i = 0; i < n; ++i) {
        // Cache particle i position and type to reduce repeated lookups
        const auto& pos_i = particles[i].getPosition();
        const int type_i = particles[i].getType();
        auto& force_i = particles[i].getForce();

        for (size_t j = i + 1; j < n; ++j) {
            // Cache particle j position and type
            const auto& pos_j = particles[j].getPosition();
            const int type_j = particles[j].getType();
            auto& force_j = particles[j].getForce();

            // Get sigma2 and epsilon from lookup table
            const double sigma2 = sigma2_table[type_i][type_j];
            const double epsilon = epsilon_table[type_i][type_j];

            // Compute the vector between particles and squared distance
            const double dx = pos_i[0] - pos_j[0];
            const double dy = pos_i[1] - pos_j[1];
            const double dz = pos_i[2] - pos_j[2];
            const double r2 = dx*dx + dy*dy + dz*dz;

            // Calculate required powers
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
 * This function iterates over all particles in the system, retrieves their
 * velocities, and computes the kinetic energy using the formula:
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
    double tempKE = (2.0 / (3.0 * kb)) * kineticEnergy;
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
    
    // Output initial conditions before starting integration
    outputParticleData(0);
    outputKineticEnergy(0);

    const double outputTime = 0.1;
    const int steps = static_cast<int>(finalTime / dt);
    const int outputStep = static_cast<int>(outputTime / dt);

    // Start simulation from step 1
    for (int step = 1; step <= steps; ++step) {
        const double currentTime = step * dt;
        calKE();
        if (testCase == -1) 
            velRescale();
        calForces();
        forwardEuler();

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