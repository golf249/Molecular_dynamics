// md.h - Molecular Dynamics Header
#ifndef MD_H
#define MD_H

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <random>

struct Particle {
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> force;
    double mass;
};

class MolecularDynamics {
public:
    MolecularDynamics(int numParticles, double dt, double Lx, double Ly, double Lz);
    void initializeParticles();
    void computeForces();
    void integrate();
    void applyBoundaryConditions();
    void runSimulation(int steps);
    void printState(int step) const;

private:
    int N;
    double dt;
    double Lx, Ly, Lz;
    std::vector<Particle> particles;
    std::mt19937 rng;
};

#endif // MD_H