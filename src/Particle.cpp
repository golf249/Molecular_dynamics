#include "../include/Particle.h"

constexpr double Properties::mass[2];

Particle::Particle(const std::array<double, 3>& position, const std::array<double, 3>& velocity, int type)
    : position(position), velocity(velocity), force({0.0, 0.0, 0.0}), mass(Properties::mass[type]), type(type) {}

const std::array<double, 3>& Particle::getPosition() const {
    return position;
}

void Particle::setPosition(const std::array<double, 3>& position) {
    this->position = position;
}

const std::array<double, 3>& Particle::getVelocity() const {
    return velocity;
}

void Particle::setVelocity(const std::array<double, 3>& velocity) {
    this->velocity = velocity;
}

std::array<double, 3>& Particle::getForce() {
    return force;
}

void Particle::setForce(const std::array<double, 3>& force) {
    this->force = force;
}

double Particle::getMass() const {
    return mass;
}

int Particle::getType() const {
    return type;
}

void Particle::setType(int type) {
    this->type = type;
}