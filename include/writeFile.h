#ifndef WRITEFILE_H
#define WRITEFILE_H

#include <fstream>
#include <string>
#include <array>

class WriteFile {
public:
    WriteFile(const std::string& particleDataFileName, const std::string& kineticEnergyFileName);
    ~WriteFile();

    void writeParticleData(double time, int particleIndex, const std::array<double, 3>& position, const std::array<double, 3>& velocity);
    void writeKineticEnergy(double time, double kineticEnergy);

private:
    std::ofstream particleDataFile;
    std::ofstream kineticEnergyFile;
};

#endif // WRITEFILE_H