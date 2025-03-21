/**
 * @file writeFile.h
 * @brief Declaration of the WriteFile class.
 *
 * This file declares the WriteFile class, which provides functionalities to write
 * particle data and kinetic energy values to output files.
 */

#ifndef WRITEFILE_H
#define WRITEFILE_H

#include <fstream>
#include <string>
#include <array>

/**
 * @class WriteFile
 * @brief Handles output file operations for particle data and kinetic energy.
 *
 * The WriteFile class opens two output files—one for particle data (positions and velocities)
 * and one for kinetic energy data. It offers methods to write formatted data to these files.
 */
class WriteFile {
public:
    /**
     * @brief Constructs a WriteFile object.
     *
     * Opens the specified files for writing particle data and kinetic energy.
     *
     * @param particleDataFileName Name of the file to write particle data.
     * @param kineticEnergyFileName Name of the file to write kinetic energy data.
     */
    WriteFile(const std::string& particleDataFileName, const std::string& kineticEnergyFileName);

    /**
     * @brief Destructor for WriteFile.
     *
     * Closes the opened output files.
     */
    ~WriteFile();

    /**
     * @brief Writes particle data to the particle data file.
     *
     * Writes the simulation time, particle index, position, and velocity of a particle
     * to the particle data file in a formatted way.
     *
     * @param time The current simulation time.
     * @param particleIndex The index of the particle.
     * @param position The particle's position as an array of 3 doubles.
     * @param velocity The particle's velocity as an array of 3 doubles.
     */
    void writeParticleData(double time, int particleIndex, const std::array<double, 3>& position, const std::array<double, 3>& velocity);

    /**
     * @brief Writes kinetic energy data to the kinetic energy file.
     *
     * Writes the simulation time and the system's kinetic energy to the kinetic energy file.
     *
     * @param time The current simulation time.
     * @param kineticEnergy The total kinetic energy of the system.
     */
    void writeKineticEnergy(double time, double kineticEnergy);

private:
    std::ofstream particleDataFile; ///< Output file stream for particle data.
    std::ofstream kineticEnergyFile; ///< Output file stream for kinetic energy data.
};

#endif // WRITEFILE_H