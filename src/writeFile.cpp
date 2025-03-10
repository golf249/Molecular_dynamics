#include "../include/writeFile.h"
#include <iostream>

WriteFile::WriteFile(const std::string& particleDataFileName, const std::string& kineticEnergyFileName) {
    particleDataFile.open(particleDataFileName);
    if (!particleDataFile.is_open()) {
        std::cerr << "Error opening particle data file: " << particleDataFileName << std::endl;
    }

    kineticEnergyFile.open(kineticEnergyFileName);
    if (!kineticEnergyFile.is_open()) {
        std::cerr << "Error opening kinetic energy file: " << kineticEnergyFileName << std::endl;
    }
}

WriteFile::~WriteFile() {
    if (particleDataFile.is_open()) {
        particleDataFile.close();
    }
    if (kineticEnergyFile.is_open()) {
        kineticEnergyFile.close();
    }
}

void WriteFile::writeParticleData(double time, int particleIndex, const std::array<double, 3>& position, const std::array<double, 3>& velocity) {
    if (particleDataFile.is_open()) {
        particleDataFile << time << " " << particleIndex << " "
                         << position[0] << " " << position[1] << " " << position[2] << " "
                         << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\n";
        particleDataFile.flush(); // Ensure data is written to the file
    } else {
        std::cerr << "Particle data file is not open." << std::endl;
    }
}

void WriteFile::writeKineticEnergy(double time, double kineticEnergy) {
    if (kineticEnergyFile.is_open()) {
        kineticEnergyFile << time << " " << kineticEnergy << "\n";
        kineticEnergyFile.flush(); // Ensure data is written to the file
    } else {
        std::cerr << "Kinetic energy file is not open." << std::endl;
    }
}