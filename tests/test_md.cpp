#include "../include/md.h"

#define BOOST_TEST_MODULE MolecularDynamicsTest
#include <boost/test/included/unit_test.hpp>

// Helper function to run a simulation and return the final state
void getFinalStates(int numParticles, double dt, double Lx, double Ly, double Lz, int testCase, double temp, double percent_type1, double finalTime, std::vector<Particle>& finalParticles, double& finalKineticEnergy) {
    MolecularDynamics sims(numParticles, dt, Lx, Ly, Lz, testCase, temp, percent_type1, finalTime);
    sims.runSimulation();
    finalParticles = sims.getParticles();
    finalKineticEnergy = sims.getKineticEnergy();
}

// Test case for one stationary particle
BOOST_AUTO_TEST_CASE(OneStationaryParticle) {
    std::vector<Particle> finalParticles;
    double finalKineticEnergy;
    runSimulationAndGetFinalState(1, 0.001, 20.0, 20.0, 20.0, 1, -1.0, 10.0, 1.0, finalParticles, finalKineticEnergy);

    BOOST_REQUIRE_EQUAL(finalParticles.size(), 1);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[0], 10.0, 1e-5);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[1], 10.0, 1e-5);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[2], 10.0, 1e-5);
    BOOST_CHECK_CLOSE(finalKineticEnergy, 0.0, 1e-5);
}

// Add more test cases for other initial conditions...
