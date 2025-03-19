#ifndef MD_H
#define MD_H

#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Particle.h"
#include "writeFile.h"

class MolecularDynamics {
public:
    MolecularDynamics(int numParticles, double dt, double Lx, double Ly, double Lz, int testCase = -1, double temp = -1.0, double percent_type1 = 10.0, double finalTime = -1.0);
    ~MolecularDynamics();
    void initialiseParticles();
    void runSimulation();
    void outputParticleData(double time);
    void outputKineticEnergy(double time);
    std::vector<Particle> getParticles() const { return particles; }
    double getKineticEnergy() const { return kineticEnergy; }

private:
    const int N;
    const double dt;
    const double Lx, Ly, Lz;
    const int testCase;
    double temp;
    double percent_type1;
    double finalTime;
    double kineticEnergy;
    std::vector<Particle> particles;
    Particle* particlesCUDA;
    WriteFile writeFile;
    double* position_d; 
    double* force_d;  
    int* type_d; 

    void assignRandomStates(const int numType, const int type) ;
    bool stabilityCheck(const std::array<double, 3>& position);
    void calForcesCUDA();
    void forwardEuler();
    void bcCheck();
    void calKE();
    void velRescale();
};
#endif