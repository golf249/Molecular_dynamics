#include "../include/md.h"

#define BOOST_TEST_MODULE MolecularDynamicsTest
#include <boost/test/included/unit_test.hpp>

// Define tolerance variable (relative tolerance in percent)
const double tol = 0.01;

// Helper function to run a simulation and return the final state
void testRun(int numParticles, double dt, double Lx, double Ly, double Lz,
             int testCase, double temp, double percent_type1, double finalTime,
             std::vector<Particle>& finalParticles, double& finalKineticEnergy) {
    MolecularDynamics sims(numParticles, dt, Lx, Ly, Lz, testCase, temp, percent_type1, finalTime);
    sims.runSimulation();
    finalParticles = sims.getParticles();
    finalKineticEnergy = sims.getKineticEnergy();
}

// Test case for one stationary particle
BOOST_AUTO_TEST_CASE(ic_one) {
    std::vector<Particle> finalParticles;
    double finalKineticEnergy;
    testRun(1, 0.001, 20.0, 20.0, 20.0, 1, -1.0, 10.0, 1.0, finalParticles, finalKineticEnergy);

    BOOST_REQUIRE_EQUAL(finalParticles.size(), 1);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[0], 10.0, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[1], 10.0, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[2], 10.0, tol);
    BOOST_CHECK_CLOSE(finalKineticEnergy, 0.0, tol);
}

// Test case for one particle with velocity
BOOST_AUTO_TEST_CASE(ic_one_vel) {
    std::vector<Particle> finalParticles;
    double finalKineticEnergy;
    testRun(1, 0.001, 20.0, 20.0, 20.0, 2, -1.0, 10.0, 20.0, finalParticles, finalKineticEnergy);

    BOOST_REQUIRE_EQUAL(finalParticles.size(), 1);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[0], 10.0, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[1], 10.0, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[2], 10.0, tol);
    BOOST_CHECK_CLOSE(finalKineticEnergy, 15, tol);
}

// Test case for two particles (ic_two)
BOOST_AUTO_TEST_CASE(ic_two) {
    std::vector<Particle> finalParticles;
    double finalKineticEnergy;
    testRun(2, 0.001, 20.0, 20.0, 20.0, 3, -1.0, 10.0, 50.0, finalParticles, finalKineticEnergy);

    BOOST_REQUIRE_EQUAL(finalParticles.size(), 2);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[0], 8.50767, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[1], 10.0, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[2], 10.0, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[0], 11.4923, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[1], 10.0, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[2], 10.0, tol);
    BOOST_CHECK_CLOSE(finalKineticEnergy, 0.000515153, tol);
}

// Test case for two particles (ic_two_pass1)
BOOST_AUTO_TEST_CASE(ic_two_pass1) {
    std::vector<Particle> finalParticles;
    double finalKineticEnergy;
    testRun(2, 0.001, 20.0, 20.0, 20.0, 4, -1.0, 10.0, 50.0, finalParticles, finalKineticEnergy);

    BOOST_REQUIRE_EQUAL(finalParticles.size(), 2);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[0], 7.249, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[1], 5.71562, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[2], 10.0, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[0], 12.751, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[1], 14.2844, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[2], 10.0, tol);
    BOOST_CHECK_CLOSE(finalKineticEnergy, 0.247955, tol);
}

// Test case for two particles (ic_two_pass2)
BOOST_AUTO_TEST_CASE(ic_two_pass2) {
    std::vector<Particle> finalParticles;
    double finalKineticEnergy;
    testRun(2, 0.001, 20.0, 20.0, 20.0, 5, -1.0, 10.0, 50.0, finalParticles, finalKineticEnergy);

    BOOST_REQUIRE_EQUAL(finalParticles.size(), 2);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[0], 12.7544, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[1], 18.1728, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[2], 10.0, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[0], 7.24565, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[1], 1.82722, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[2], 10.0, tol);
    BOOST_CHECK_CLOSE(finalKineticEnergy, 0.246937, tol);
}

// Test case for two particles (ic_two_pass3)
BOOST_AUTO_TEST_CASE(ic_two_pass3) {
    std::vector<Particle> finalParticles;
    double finalKineticEnergy;
    testRun(2, 0.001, 20.0, 20.0, 20.0, 6, -1.0, 10.0, 50.0, finalParticles, finalKineticEnergy);

    BOOST_REQUIRE_EQUAL(finalParticles.size(), 2);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[0], 10.5965, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[1], 8.10445, tol);
    BOOST_CHECK_CLOSE(finalParticles[0].getPosition()[2], 10.0, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[0], 9.40349, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[1], 11.8956, tol);
    BOOST_CHECK_CLOSE(finalParticles[1].getPosition()[2], 10.0, tol);
    BOOST_CHECK_CLOSE(finalKineticEnergy, 2.28109, tol);
}