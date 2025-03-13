#include "../include/md.h"
#include <iostream>
#include <cstdlib>
#include <boost/program_options.hpp>
#include <fstream>
#include <chrono> 

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
    // Parse command line options
    po::options_description opts("Available options");
    opts.add_options()
        ("help", "Print available options.")
        ("Lx", po::value<double>()->default_value(20.0), "x length (Angstrom)")
        ("Ly", po::value<double>()->default_value(20.0), "y length (Angstrom)")
        ("Lz", po::value<double>()->default_value(20.0), "z length (Angstrom)")
        ("dt", po::value<double>()->default_value(0.001), "Time-step")
        ("T", po::value<double>(), "Final time")
        ("ic-one", "Initial condition: one stationary particle")
        ("ic-one-vel", "Initial condition: one moving particle")
        ("ic-two", "Initial condition: two bouncing particles")
        ("ic-two-pass1", "Initial condition: two passing particles")
        ("ic-two-pass2", "Initial condition: two passing particles close")
        ("ic-two-pass3", "Initial condition: two passing particles close, heavy")
        ("ic-random", "Initial condition: N random particles")
        ("percent-type1", po::value<double>()->default_value(10.0), "Percentage of type 1 particles with random IC")
        ("N", po::value<int>(), "Number of particles to spawn with random IC")
        ("temp", po::value<double>(), "Temperature (degree Kelvin)")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << opts << std::endl;
        return 0;
    }

    // Check if the user provided any IC option:
    if (!(vm.count("ic-one") || vm.count("ic-one-vel") || vm.count("ic-two") ||
    vm.count("ic-two-pass1") || vm.count("ic-two-pass2") || vm.count("ic-two-pass3") ||
    vm.count("ic-random"))) {
    std::cerr << "Error: No initial condition (--ic option) provided." << std::endl;
    return 1;
    }

    // Check that only one --ic option is provided by incrementing ic_count
    int ic_count = 0;
    int testCase = -1;
    if (vm.count("ic-one")) { ic_count++; testCase = 1; }
    if (vm.count("ic-one-vel")) { ic_count++; testCase = 2; }
    if (vm.count("ic-two")) { ic_count++; testCase = 3; }
    if (vm.count("ic-two-pass1")) { ic_count++; testCase = 4; }
    if (vm.count("ic-two-pass2")) { ic_count++; testCase = 5; }
    if (vm.count("ic-two-pass3")) { ic_count++; testCase = 6; }
    if (vm.count("ic-random")) { ic_count++; testCase = -1; }

    // Check that exactly one --ic option is provided
    if ((ic_count != 1) && (!vm.count("T"))) {
        std::cerr << "Error: Exactly one --ic option must be provided." << std::endl;
        return 1;
    }

    // Check that --T is provided for the "ic-" test cases
    if (testCase == 1 || testCase == 2 || testCase == 3 || testCase == 4 || testCase == 5 || testCase == 6) {
        if (!vm.count("T")) {
            std::cerr << "Error: --T (Final time) must be provided for the \"ic-\" test case." << std::endl;            return 1;
        }
    }

    // Check that --N and --percent-type1 are provided if --ic-random is specified
    if (vm.count("ic-random")) {
        if ((!vm.count("N")) || (!vm.count("T"))) {
            std::cerr << "Error: Both --N and --T must be provided when --ic-random is specified." << std::endl;
            return 1;
        }
    }
    
    const double Lx = vm["Lx"].as<double>();
    const double Ly = vm["Ly"].as<double>();
    const double Lz = vm["Lz"].as<double>();
    const double dt = vm["dt"].as<double>();
    double T = vm.count("T") ? vm["T"].as<double>() : -1.0; // -1 means not set
    double temp = vm.count("temp") ? vm["temp"].as<double>() : -1.0; // -1 means not set
    int N = vm.count("N") ? vm["N"].as<int>() : -1; // -1 means not set
    double percent_type1 = vm["percent-type1"].as<double>();

    // Print parsed values
    std::cout << "Simulation parameters:\n"
              << "Box Size: (" << Lx << ", " << Ly << ", " << Lz << ")\n"
              << "Time Step: " << dt << "\n"
              << "Final Time: " << (T >= 0 ? std::to_string(T) : "Not Fixed") << "\n"
              << "Temperature: " << (temp >= 0 ? std::to_string(temp) : "Not Fixed") << "\n"
              << "Number of Particles: " << (N >= 0 ? std::to_string(N) : "Not Fixed") << "\n"
              << "Type 1 Particle Percentage: " << percent_type1 << "%\n"
              << "Test Case: " << (testCase == -1 ? ("Random") : std::to_string(testCase))<< "\n";
    std::cout << "Starting simulation" << std::endl;

    if (vm.count("ic-one") || vm.count("ic-one-vel") || vm.count("ic-two") || vm.count("ic-two-pass1") || vm.count("ic-two-pass2") || vm.count("ic-two-pass3")) {
        // Start time measurement.
        auto start = std::chrono::high_resolution_clock::now();
        
        // Initialize and run the simulation with the specified test case
        MolecularDynamics simulation(N, dt, Lx, Ly, Lz, testCase, temp, percent_type1, T);
        simulation.runSimulation();

        // End time measurement.
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        double duration_seconds = duration / 1000.0;
        std::cout << "Simulation took " << duration_seconds << " seconds." << std::endl;
        std::cout << "Relevant files stored!." << std::endl;
        std::cout << "Simulation complete." << std::endl;
        
    } else if (vm.count("ic-random")) {
        // Start time measurement.
        auto start = std::chrono::high_resolution_clock::now();

        // Initialize and run the simulation with random initial conditions
        MolecularDynamics simulation(N, dt, Lx, Ly, Lz, -1, temp, percent_type1, T);
        simulation.runSimulation();

        // End time measurement.
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        double duration_seconds = duration / 1000.0;
        std::cout << "Simulation took " << duration_seconds << " seconds." << std::endl;
        std::cout << "Relevant files stored!." << std::endl;
        std::cout << "Simulation complete." << std::endl;
    }

    return 0;
}