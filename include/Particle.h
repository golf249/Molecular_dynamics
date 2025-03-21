/**
 * @file Particle.h
 * @brief Declaration of the Particle class and Properties structure.
 *
 * This header declares the Particle class which represents a particle in a molecular dynamics simulation.
 * It includes functionalities for particle initialization, retrieving and updating particle properties such as
 * position, velocity, force, mass, and type. Additionally, it declares the Properties structure which provides
 * constants such as the masses for different particle types.
 */

#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>

/**
 * @brief Contains properties for particles.
 *
 * This structure provides constants such as the masses for different particle types.
 */
struct Properties {
    /// Mass values for each particle type. For example, type 0 has mass 1.0 and type 1 mass 10.0.
    static constexpr double mass[2] = {1.0, 10.0};
};

/**
 * @class Particle
 * @brief Represents a particle in a molecular dynamics simulation.
 *
 * The Particle class encapsulates properties of a particle, such as position, velocity,
 * force, mass, and type. It provides methods to retrieve and update these properties.
 */
class Particle {
public:
    /**
     * @brief Constructs a Particle.
     *
     * Initializes the particle with the specified position, velocity, and type.
     * The mass is determined based on the provided type using the Properties structure.
     *
     * @param position The initial position as an array of 3 doubles.
     * @param velocity The initial velocity as an array of 3 doubles.
     * @param type The particle type identifier.
     */
    Particle(const std::array<double, 3>& position, const std::array<double, 3>& velocity, int type);

    /**
     * @brief Retrieves the particle's position.
     *
     * @return A constant reference to the position array.
     */
    const std::array<double, 3>& getPosition() const;

    /**
     * @brief Sets the particle's position.
     *
     * @param position The new position as an array of 3 doubles.
     */
    void setPosition(const std::array<double, 3>& position);

    /**
     * @brief Retrieves the particle's velocity.
     *
     * @return A constant reference to the velocity array.
     */
    const std::array<double, 3>& getVelocity() const;

    /**
     * @brief Sets the particle's velocity.
     *
     * @param velocity The new velocity as an array of 3 doubles.
     */
    void setVelocity(const std::array<double, 3>& velocity);

    /**
     * @brief Retrieves the particle's force.
     *
     * @return A reference to the force array.
     */
    std::array<double, 3>& getForce();

    /**
     * @brief Sets the particle's force.
     *
     * @param force The new force as an array of 3 doubles.
     */
    void setForce(const std::array<double, 3>& force);

    /**
     * @brief Retrieves the particle's mass.
     *
     * @return The mass of the particle.
     */
    double getMass() const;
    
    /**
     * @brief Retrieves the particle's type.
     *
     * @return The particle type identifier.
     */
    int getType() const;

    /**
     * @brief Sets the particle's type.
     *
     * @param type The new type identifier.
     */
    void setType(int type);

private:
    std::array<double, 3> position; /**< The position of the particle. */
    std::array<double, 3> velocity; /**< The velocity of the particle. */
    std::array<double, 3> force;    /**< The force acting on the particle. */
    double mass;                    /**< The mass of the particle. */
    int type;                       /**< The type of particle */
};

#endif // PARTICLE_H