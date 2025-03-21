/**
 * @file writeFile.cpp
 * @brief Implements the WriteFile class for outputting simulation data.
 *
 * This file contains the definitions for the WriteFile class methods. The class handles
 * opening and closing output file streams, writing particle data (positions and velocities),
 * and writing kinetic energy values to separate output files.
 */

#include "../include/writeFile.h"
#include <iostream>

/**
 * @brief Constructs a WriteFile object.
 *
 * Opens the output files for particle data and kinetic energy.
 *
 * @param particleDataFileName Name of the file in which to write particle data.
 * @param kineticEnergyFileName Name of the file in which to write kinetic energy data.
 */
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

/**
 * @brief Destroys the WriteFile object.
 *
 * Closes any open file streams for particle data and kinetic energy.
 */
WriteFile::~WriteFile() {
    if (particleDataFile.is_open()) {
        particleDataFile.close();
    }
    if (kineticEnergyFile.is_open()) {
        kineticEnergyFile.close();
    }
}

/**
 * @brief Writes particle data to the designated output file.
 *
 * The function writes the simulation time, particle index, position and velocity data.
 * It flushes the file stream to ensure that data is immediately written to the file.
 *
 * @param time The current simulation time.
 * @param particleIndex The index of the particle.
 * @param position The particle's position as a 3-element array.
 * @param velocity The particle's velocity as a 3-element array.
 */
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

/**
 * @brief Writes kinetic energy data to the designated output file.
 *
 * The function writes the simulation time and associated kinetic energy value.
 * It flushes the file stream to ensure immediate output.
 *
 * @param time The current simulation time.
 * @param kineticEnergy The total kinetic energy of the system.
 */
void WriteFile::writeKineticEnergy(double time, double kineticEnergy) {
    if (kineticEnergyFile.is_open()) {
        kineticEnergyFile << time << " " << kineticEnergy << "\n";
        kineticEnergyFile.flush(); // Ensure data is written to the file
    } else {
        std::cerr << "Kinetic energy file is not open." << std::endl;
    }
}