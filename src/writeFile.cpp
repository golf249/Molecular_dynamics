#include "../include/writeFile.h"

WriteFile::WriteFile(const std::string& particleDataFileName, const std::string& kineticEnergyFileName) {
    particleDataFile.open(particleDataFileName);
    kineticEnergyFile.open(kineticEnergyFileName);
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
    particleDataFile << time << " " << particleIndex << " "
                     << position[0] << " " << position[1] << " " << position[2] << " "
                     << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\n";
}

void WriteFile::writeKineticEnergy(double time, double kineticEnergy) {
    kineticEnergyFile << time << " " << kineticEnergy << "\n";
}