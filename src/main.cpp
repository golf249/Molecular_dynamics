#include "../include/md.h"
#include <iostream>
#include <cstdlib>
#include <boost/program_options.hpp>

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

    // Check that only one --ic option is provided
    int ic_count = 0;
    if (vm.count("ic-one")) ic_count++;
    if (vm.count("ic-one-vel")) ic_count++;
    if (vm.count("ic-two")) ic_count++;
    if (vm.count("ic-two-pass1")) ic_count++;
    if (vm.count("ic-two-pass2")) ic_count++;
    if (vm.count("ic-two-pass3")) ic_count++;
    if (vm.count("ic-random")) ic_count++;

    if (ic_count != 1) {
        std::cerr << "Error: Exactly one --ic option must be provided." << std::endl;
        return 1;
    }

    // Check that --N and --percent-type1 are provided if --ic-random is specified
    if (vm.count("ic-random")) {
        if (!vm.count("N")) {
            std::cerr << "Error: --N must be provided when --ic-random is specified." << std::endl;
            return 1;
        }
    }
    
    const double Lx = vm["Lx"].as<double>();
    const double Ly = vm["Ly"].as<double>();
    const double Lz = vm["Lz"].as<double>();
    const double dt = vm["dt"].as<double>();
    double T = vm.count("T") ? vm["T"].as<double>() : -1.0; // -1 means not set
    double temp = vm.count("temp") ? vm["temp"].as<double>() : -1.0; // -1 means not set
    int N = vm.count("N") ? vm["N"].as<int>() : 0;
    double percent_type1 = vm["percent-type1"].as<double>(); // Default is fine



    // Print parsed values (for debugging)
    std::cout << "Simulation parameters:\n"
              << "Box Size: (" << Lx << ", " << Ly << ", " << Lz << ")\n"
              << "Time Step: " << dt << "\n"
              << "Final Time: " << (T >= 0 ? std::to_string(T) : "Not Fixed") << "\n"
              << "Temperature: " << (temp >= 0 ? std::to_string(temp) : "Not Fixed") << "\n"
              << "Number of Particles: " << N << "\n"
              << "Type 1 Particle Percentage: " << percent_type1 << "%\n";
    


    return 0;
}

    // Initialize and run the simulation
    // MolecularDynamics simulation(numParticles, dt, Lx, Ly, Lz);
    // simulation.runSimulation(steps);
