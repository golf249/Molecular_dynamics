// md.cpp - Molecular Dynamics Implementation

#include "../include/md.h"

MolecularDynamics::MolecularDynamics(int numParticles, double dt, double Lx, double Ly, double Lz)
    : N(numParticles), dt(dt), Lx(Lx), Ly(Ly), Lz(Lz), uniform_dist(0.0, 1.0) {
    initializeRandomEngine();
    initializeParticles();
}

void MolecularDynamics::initializeRandomEngine() {
    std::random_device rd;
    rng = std::mt19937(rd());
}

void MolecularDynamics::initializeParticles() {
    particles.resize(N);
    for (auto& p : particles) {
        p.position = {Lx * uniform_dist(rng), Ly * uniform_dist(rng), Lz * uniform_dist(rng)};
        p.velocity = {uniform_dist(rng) - 0.5, uniform_dist(rng) - 0.5, uniform_dist(rng) - 0.5};
        p.force = {0.0, 0.0, 0.0};
        p.mass = (uniform_dist(rng) < 0.1) ? 10.0 : 1.0; // 10% probability for type 1 particles
    }
}

void MolecularDynamics::computeForces() {
    const double sigma = 1.0;
    const double epsilon = 3.0;
    for (auto& p : particles) {
        p.force = {0.0, 0.0, 0.0};
    }
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            std::array<double, 3> rij;
            double r2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                rij[k] = particles[j].position[k] - particles[i].position[k];
                r2 += rij[k] * rij[k];
            }
            double r6 = r2 * r2 * r2;
            double r12 = r6 * r6;
            double f = 24 * epsilon * (2 * sigma * sigma * sigma * sigma * sigma * sigma / r12 - sigma * sigma * sigma / r6) / r2;
            for (int k = 0; k < 3; ++k) {
                double forceComponent = f * rij[k];
                particles[i].force[k] += forceComponent;
                particles[j].force[k] -= forceComponent;
            }
        }
    }
}

void MolecularDynamics::integrate() {
    for (auto& p : particles) {
        for (int k = 0; k < 3; ++k) {
            p.velocity[k] += dt * p.force[k] / p.mass;
            p.position[k] += dt * p.velocity[k];
        }
    }
    applyBoundaryConditions();
}

void MolecularDynamics::applyBoundaryConditions() {
    for (auto& p : particles) {
        for (int k = 0; k < 3; ++k) {
            if (p.position[k] < 0) {
                p.position[k] = -p.position[k];
                p.velocity[k] = std::abs(p.velocity[k]);
            } else if (p.position[k] > ((k == 0) ? Lx : (k == 1) ? Ly : Lz)) {
                p.position[k] = 2 * ((k == 0) ? Lx : (k == 1) ? Ly : Lz) - p.position[k];
                p.velocity[k] = -std::abs(p.velocity[k]);
            }
        }
    }
}

void MolecularDynamics::runSimulation(int steps) {
    for (int step = 0; step < steps; ++step) {
        computeForces();
        integrate();
        if (step % 10 == 0) printState(step);
    }
}

void MolecularDynamics::printState(int step) const {
    std::cout << "Step " << step << ":\n";
    for (const auto& p : particles) {
        std::cout << "Position: (" << p.position[0] << ", " << p.position[1] << ", " << p.position[2] << ") ";
        std::cout << "Velocity: (" << p.velocity[0] << ", " << p.velocity[1] << ", " << p.velocity[2] << ")\n";
    }
}
