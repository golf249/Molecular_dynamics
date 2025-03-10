#ifndef MD_H
#define MD_H

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Particle.h"

class MolecularDynamics {
public:
    MolecularDynamics(int numParticles, double dt, double Lx, double Ly, double Lz, int testCase = -1, double temp = -1.0, double percent_type1 = 10.0, double finalTime = -1.0);
    void initializeParticles();
    void computeForces();
    void integrate();
    void applyBoundaryConditions();
    void runSimulation();
    void outputParticleData(double time);
    void outputKineticEnergy(double time);

private:
    const int N;
    const double dt;
    const double Lx, Ly, Lz;
    const int testCase;
    double temp;
    double percent_type1;
    double finalTime;
    std::vector<Particle> particles;
    std::ofstream particleDataFile;
    std::ofstream kineticEnergyFile;
};

#endif // MD_H