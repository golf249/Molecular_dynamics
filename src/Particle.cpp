/**
 * @file Particle.cpp
 * @brief Implementation of the Particle class for the MolecularDynamics simulation.
 *
 * This file contains the implementation of the Particle class, which represents
 * a particle in the MolecularDynamics simulation. The class includes methods for
 * setting and getting the particle's position, velocity, force, mass, and type.
 */

#include "../include/Particle.h"

/**
 * @brief Array of mass values for different particle types.
 * 
 * This constexpr array holds the mass values for two different types of particles.
 */
constexpr double Properties::mass[2];

/**
 * @brief Constructs a new Particle object.
 * 
 * @param position The initial position of the particle as a 3-element array.
 * @param velocity The initial velocity of the particle as a 3-element array.
 * @param type The type of the particle, used to determine its mass.
 */
Particle::Particle(const std::array<double, 3>& position, const std::array<double, 3>& velocity, int type)
    : position(position), velocity(velocity), force({0.0, 0.0, 0.0}), mass(Properties::mass[type]), type(type) {}

/**
 * @brief Get the position of the particle.
 * 
 * This function returns a constant reference to an array containing 
 * the position components of the particle.
 * 
 * @return const std::array<double, 3>& A constant reference to the 
 *         array of the particle's position.
 */
const std::array<double, 3>& Particle::getPosition() const {
    return position;
}

/**
 * @brief Sets the position of the particle.
 * 
 * @param position A std::array containing the x, y, and z coordinates of the particle's position.
 */
void Particle::setPosition(const std::array<double, 3>& position) {
    this->position = position;
}

/**
 * @brief Get the velocity of the particle.
 * 
 * This function returns a constant reference to an array containing 
 * the velocity components of the particle.
 * 
 * @return const std::array<double, 3>& A constant reference to the velocity array.
 */
const std::array<double, 3>& Particle::getVelocity() const {
    return velocity;
}

/**
 * @brief Sets the velocity of the particle.
 *
 * @param velocity A std::array containing the velocity components in 
 * the x, y, and z directions.
 */
void Particle::setVelocity(const std::array<double, 3>& velocity) {
    this->velocity = velocity;
}

/**
 * @brief Retrieves the force vector of the particle.
 * 
 * This function returns a reference to the force vector, which is represented
 * as a std::array of three double values.
 * 
 * @return A reference to the force vector of the particle.
 */
std::array<double, 3>& Particle::getForce() {
    return force;
}

/**
 * @brief Sets the force vector for the particle.
 * 
 * @param force A 3-element array representing the force vector to be applied to the particle.
 */
void Particle::setForce(const std::array<double, 3>& force) {
    this->force = force;
}

/**
 * @brief Get the mass of the particle.
 * 
 * @return double The mass of the particle.
 */
double Particle::getMass() const {
    return mass;
}

/**
 * @brief Get the type of the particle.
 * 
 * @return int The type of the particle.
 */
int Particle::getType() const {
    return type;
}

/**
 * @brief Sets the type of the particle.
 * 
 * @param type An integer representing the type of the particle.
 */
void Particle::setType(int type) {
    this->type = type;
}