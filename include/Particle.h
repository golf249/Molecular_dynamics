#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>

struct Properties {
    static constexpr double mass[2] = {1.0, 10.0};
};

class Particle {
public:
    Particle(const std::array<double, 3>& position, const std::array<double, 3>& velocity, int type);

    const std::array<double, 3>& getPosition() const;
    void setPosition(const std::array<double, 3>& position);

    const std::array<double, 3>& getVelocity() const;
    void setVelocity(const std::array<double, 3>& velocity);

    const std::array<double, 3>& getForce() const;
    void setForce(const std::array<double, 3>& force);

    double getMass() const;
    
    int getType() const;
    void setType(int type);

private:
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> force;
    double mass;
    int type;
};

#endif // PARTICLE_H