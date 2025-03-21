/**
 * @file md.h
 * @brief Declaration of the MolecularDynamics class.
 *
 * This header declares the MolecularDynamics class which performs
 * molecular dynamics simulations. It includes functionalities for
 * particle initialization, simulation execution, force and kinetic
 * energy calculations, and writing output data.
 */

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

/**
 * @class MolecularDynamics
 * @brief Class to simulate molecular dynamics.
 *
 * The MolecularDynamics class encapsulates the simulation of particles
 * interacting under the Lennard-Jones potential. It provides methods to
 * initialize particles, compute forces (both with CUDA and CPU via OpenMP),
 * update positions and velocities, and output simulation data.
 */
class MolecularDynamics {
public:
    /**
     * @brief Constructor for MolecularDynamics.
     * 
     * Initializes the simulation parameters and allocates resources.
     *
     * @param numParticles Number of particles involved in the simulation.
     * @param dt Time-step used in the simulation.
     * @param Lx Length of the simulation domain in the x-direction.
     * @param Ly Length of the simulation domain in the y-direction.
     * @param Lz Length of the simulation domain in the z-direction.
     * @param testCase Test case identifier to select initial conditions. Defaults to -1.
     * @param temp Desired temperature of the system. Defaults to -1.0.
     * @param percent_type1 Percentage of particles of type 1. Defaults to 10.0.
     * @param finalTime Final time for the simulation. Defaults to -1.0.
     */
    MolecularDynamics(int numParticles, double dt, double Lx, double Ly, double Lz,
         int testCase = -1, double temp = -1.0, double percent_type1 = 10.0, double finalTime = -1.0);

    /**
     * @brief Destructor for MolecularDynamics.
     *
     * Frees allocated GPU memory and other resources.
     */
    ~MolecularDynamics();

    /**
     * @brief Initializes the particles for the simulation.
     *
     * Clears any existing particles and creates new particles based on the
     * specified test case, ensuring proper initialization of position, velocity,
     * mass, and type.
     */
    void initialiseParticles();

    /**
     * @brief Runs the molecular dynamics simulation.
     *
     * This function performs the time integration, including calculating forces, 
     * updating positions and velocities, applying boundary conditions, and 
     * writing output data.
     */
    void runSimulation();

    /**
     * @brief Outputs particle data for each simulation step.
     *
     * Writes particle positions and velocities to a file for a given simulation time.
     * 
     * @param time The simulation time at which the data is recorded.
     */
    void outputParticleData(double time);

    /**
     * @brief Outputs the total kinetic energy of the system.
     *
     * Writes the kinetic energy to a file for a given simulation time.
     *
     * @param time The simulation time at which the kinetic energy is recorded.
     */
    void outputKineticEnergy(double time);

    /**
     * @brief Retrieves a list of particles.
     * 
     * This function returns a vector containing all the particles.
     * 
     * @return std::vector<Particle> A vector of Particle objects.
     */
    std::vector<Particle> getParticles() const { return particles; };

    /**
     * @brief Returns the current total kinetic energy.
     *
     * @return The kinetic energy as a double.
     */
    double getKineticEnergy() const { return kineticEnergy; }

private:
    const int N;              /**< Number of particles. */
    const double dt;          /**< Time-step for the simulation. */
    const double Lx, Ly, Lz;  /**< Dimensions of the simulation domain. */
    const int testCase;       /**< Test case identifier for initialization. */
    double temp;              /**< Desired temperature of the system. */
    double percent_type1;     /**< Percentage of particles of type 1. */
    double finalTime;         /**< Final time for the simulation. */
    double kineticEnergy;     /**< Current kinetic energy of the system. */
    std::vector<Particle> particles;  /**< Container for particles. */
    WriteFile writeFile;      /**< Utility class for writing output files. */
    double* position_d;       /**< Pointer to GPU managed memory for particle positions. */
    double* force_d;          /**< Pointer to GPU managed memory for forces. */
    int* type_d;              /**< Pointer to GPU managed memory for particle types. */

    /**
     * @brief Assigns random states to particles of given type.
     *
     * This method generates random positions and velocities for particles,
     * ensuring that particle positions maintain a minimum distance via a
     * stability check.
     *
     * @param numType Number of particles to assign a state.
     * @param type The type (identifier) of particles.
     */
    void assignRandomStates(const int numType, const int type);

    /**
     * @brief Checks the stability of a particle position.
     *
     * Determines if the provided position is valid by ensuring particles are 
     * not too close to each other.
     * 
     * @param position The position to be checked.
     * @return True if the position is stable (valid); false otherwise.
     */
    bool stabilityCheck(const std::array<double, 3>& position);

    /**
     * @brief Calculates the forces acting on particles in the system.
     *
     * This function computes the forces based on the current positions and 
     * interactions of the particles. It updates the force vectors for each 
     * particle accordingly.
     */
    void calForces() ;

    /**
     * @brief Calculates the forces in parallel using openMP.
     *
     * This function is responsible for computing the forces acting on particles
     * using parallelisation with OpenMP.
     */
    void calForcesParallel();

    /**
     * @brief Calculates forces on particles using CUDA.
     *
     * Copies particle data to the GPU, launches a CUDA kernel to compute 
     * interactions based on the Lennard-Jones potential, and updates the force array.
     */
    void calForcesCUDA();

    /**
     * @brief Updates particle positions and velocities using the forward Euler method.
     */
    void forwardEuler();

    /**
     * @brief Enforces boundary conditions on particle positions.
     *
     * Checks and corrects any particles that leave the simulation domain.
     */
    void bcCheck();

    /**
     * @brief Calculates the total kinetic energy of the system.
     */
    void calKE();

    /**
     * @brief Rescales particle velocities to maintain the desired temperature.
     *
     * If a target temperature is set, the velocities of all particles are adjusted 
     * by a scaling factor.
     */
    void velRescale();

    /**
     * @brief Updates particle forces from the GPU computed forces.
     *
     * Copies the forces computed on the GPU back to the particle container.
     */
    void setParticleForces();
};

#endif