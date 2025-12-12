#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "configuration.hpp"

int main(int argc, char* argv[]) {

    if (argc != 7){
        std::cerr << "Usage: " << argv[0] << " <Lattice Size N> <Temperature T> <Magnetic Field H> <Max steps> <Burn in> <Thinning>\n";
        return EXIT_FAILURE;
    }

    // Simulation parameters
    int MAX_TIME = std::stoi(argv[4]);   
    int THINNING = std::stoi(argv[6]);
    int BURNIN = std::stoi(argv[5]);
    // Define sizes, temperature and magnetic field
    int N = std::stoi(argv[1]);
    double H = std::stod(argv[3]);
    double T = std::stod(argv[2]);
    // Initialization probability
    double p_init = 0.5;
    // Build the lattice
    SpinConfiguration config(N, T, p_init, H);
    // Attach a MCMC engine
    config.mountMCMCengine(MCMCType::MetropolisHastings);
    // Set periodic boundary conditions
    config.setPBC(true);  
    // Save both magnetization and energy
    config.keepTrack(true, true); 
    // Now run the MCMC chain
    config.run(MAX_TIME, BURNIN, THINNING);
    // Get the results
    std::vector<double> energies = config.getResults().getEnergy();
    std::vector<double> magnetizations = config.getResults().getMagnetizations();
}
