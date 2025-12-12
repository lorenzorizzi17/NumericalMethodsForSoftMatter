#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "configuration.hpp"

int main(int argc, char* argv[]) {
    int MAX_TIME;
    int BLUR;

    double T;
    double H;
    int N;

    if (argc == 6){ //Blur provided
        // Simulation parameters
        MAX_TIME = std::stoi(argv[4]);   
        BLUR = std::stoi(argv[5]);
        // Define sizes, temperature and magnetic field
        N = std::stoi(argv[1]);
        H = std::stod(argv[3]);
        T = std::stod(argv[2]);
    } else if (argc == 5){ // No blur provided
        // Simulation parameters
        MAX_TIME = std::stoi(argv[4]);   
        BLUR = 1; // Default blur
        // Define sizes, temperature and magnetic field
        N = std::stoi(argv[1]);
        H = std::stod(argv[3]);
        T = std::stod(argv[2]);
    } else {
        std::cerr << "Usage: " << argv[0] << " <Lattice Size N> <Temperature T> <Magnetic Field H> <Max steps> (eventually <Blur>, default 1)\n";
        return EXIT_FAILURE;
    }
    // Initialization probability
    double p_init = 0.5;
    // Build the lattice
    SpinConfiguration config(N, T, p_init, H);
    // Attach a MCMC engine
    config.mountMCMCengine(MCMCType::MetropolisHastings);
    // Set periodic boundary conditions
    config.setPBC(true);  
    // In the graphics mode we do not keep track of magnetization and energy
    // Now run the MCMC chain
    config.runGraphics(MAX_TIME, BLUR);
}