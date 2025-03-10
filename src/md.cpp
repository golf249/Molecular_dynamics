#include "../include/md.h"

MolecularDynamics::MolecularDynamics(int numParticles, double dt, double Lx, double Ly, double Lz, int testCase, double temp, double percent_type1, double finalTime)
    : N(numParticles), dt(dt), Lx(Lx), Ly(Ly), Lz(Lz), testCase(testCase), temp(temp), percent_type1(percent_type1), finalTime(finalTime) {
    initializeParticles();
    if (testCase != -1) {
        particleDataFile.open("particle_data.txt");
    }
    kineticEnergyFile.open("kinetic_energy.txt");

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
        for (int i = 0; i < N; ++i) {
            std::array<double, 3> position = {Lx * (double)(rand()) / RAND_MAX, Ly * (double)(rand()) / RAND_MAX, Lz * (double)(rand()) / RAND_MAX};
            std::array<double, 3> velocity = {(double)(rand()) / RAND_MAX - 0.5, (double)(rand()) / RAND_MAX - 0.5, (double)(rand()) / RAND_MAX - 0.5};
            int type = ((double)(rand()) / RAND_MAX < percent_type1 / 100.0) ? 1 : 0; // percent_type1 probability for type 1 particles
            particles.emplace_back(position, velocity, type); // append a new particle to the particles vector with the given position, velocity, and type
        }
    }
}

void MolecularDynamics::computeForces() {
    for (auto& p : particles) {
        p.setForce({0.0, 0.0, 0.0});
    }
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            int type_i = particles[i].getType();
            int type_j = particles[j].getType();
            double sigma2, epsilon;

            if (type_i == 0 && type_j == 0) {
                sigma2 = 1.0;
                epsilon = 3.0;
            } else if ((type_i == 0 && type_j == 1) || (type_i == 1 && type_j == 0)) {
                sigma2 = 4.0;
                epsilon = 15.0;
            } else if (type_i == 1 && type_j == 1) {
                sigma2 = 9.0;
                epsilon = 60.0;
            }

            std::array<double, 3> rij;
            double r2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                rij[k] = particles[i].getPosition()[k] - particles[j].getPosition()[k];
                r2 += rij[k] * rij[k];
            }
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            double f = 24.0 * epsilon * (2.0 * sigma2 * sigma2 * sigma2 * sigma2 * sigma2 * sigma2 / r12 - sigma2 * sigma2 * sigma2 / r6) / r2;
            for (int k = 0; k < 3; ++k) {
                double forceComponent = f * rij[k];
                std::array<double, 3> force_i = particles[i].getForce();
                std::array<double, 3> force_j = particles[j].getForce();
                force_i[k] += forceComponent;
                force_j[k] -= forceComponent;
                particles[i].setForce(force_i);
                particles[j].setForce(force_j);
            }
        }
    }
}

void MolecularDynamics::integrate() {
    for (Particle& p : particles) {
        std::array<double, 3> velocity = p.getVelocity();
        std::array<double, 3> position = p.getPosition();
        std::array<double, 3> force = p.getForce();
        double mass = p.getMass();
        for (int k = 0; k < 3; ++k) {
            velocity[k] += dt * force[k] / mass;
            position[k] += dt * velocity[k];
        }
        p.setVelocity(velocity);
        p.setPosition(position);
    }
    applyBoundaryConditions();
}

void MolecularDynamics::applyBoundaryConditions() {
    for (Particle& p : particles) {
        std::array<double, 3> position = p.getPosition();
        std::array<double, 3> velocity = p.getVelocity();
        for (int k = 0; k < 3; ++k) {
            if (position[k] < 0) {
                position[k] = -position[k];
                velocity[k] = std::abs(velocity[k]);
            } else if (position[k] > ((k == 0) ? Lx : (k == 1) ? Ly : Lz)) {
                position[k] = 2 * ((k == 0) ? Lx : (k == 1) ? Ly : Lz) - position[k];
                velocity[k] = -std::abs(velocity[k]);
            }
        }
        p.setPosition(position);
        p.setVelocity(velocity);
    }
}

void MolecularDynamics::runSimulation() {
    int steps = static_cast<int>(finalTime / dt);
    for (int step = 0; step <= steps; ++step) {
        double currentTime = step * dt;
        computeForces();
        integrate();
        if (step % static_cast<int>(0.1 / dt) == 0) {
            if (testCase != -1) {
                outputParticleData(currentTime);
            }
            outputKineticEnergy(currentTime);
        }
    }
    // Close the files after the simulation
    if (particleDataFile.is_open()) {
        particleDataFile.close();
    }
    if (kineticEnergyFile.is_open()) {
        kineticEnergyFile.close();
    }
}

void MolecularDynamics::outputParticleData(double time) {
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        const auto& position = p.getPosition();
        const auto& velocity = p.getVelocity();
        particleDataFile << time << " " << i << " "
                         << position[0] << " " << position[1] << " " << position[2] << " "
                         << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\n";
    }
}

void MolecularDynamics::outputKineticEnergy(double time) {
    double kineticEnergy = 0.0;
    for (const auto& p : particles) {
        const auto& velocity = p.getVelocity();
        double speedSquared = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];
        kineticEnergy += 0.5 * p.getMass() * speedSquared;
    }
    kineticEnergyFile << time << " " << kineticEnergy << "\n";
}