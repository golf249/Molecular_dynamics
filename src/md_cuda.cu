#include "../include/md.h"

__global__ void calculateLJForces(double* position, double* force, double* type, int n) {
    // Get the index of the particle
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        for (int j = i+1; j < n; ++j) {
            
        }
    }
}

void MolecularDynamics::calForcesCUDA() {
    const int n = particles.size();

    // Allocate GPU memory on the device
    double* position;
    double* force;
    double* type;
    cudaMallocManaged(&position, 3 * n * sizeof(double));
    cudaMallocManaged(&force, 3 * n * sizeof(double));
    cudaMallocManaged(&type, n * sizeof(double));

    // Populate the arrays in the allocated memory
    for (int i = 0; i < n; i++) {
        const std::array<double, 3>& pos = particles[i].getPosition();
        // Get the position of the particle
        position[3 * i] = pos[0];
        position[3 * i + 1] = pos[1];
        position[3 * i + 2] = pos[2];
        // Set the forces for each particle to 0
        force[3 * i] = 0.0;
        force[3 * i + 1] = 0.0;
        force[3 * i + 2] = 0.0;
        // Get the type of the particle
        type[i] = particles[i].getType();
    }

    // Calculate num thread and thread blocks
    int threads = std::min(256, n);
    int blocks = std::max(n/256, 1);

    // Launch the kernal to compute the Leonard-Jones forces
    calculateLJForces<<<blocks, threads>>>(position, force, type, n);
}