#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <omp.h>

#include "configuration.hpp"

int main() {
    // Simulation parameters
    double H = 0;
    int MAX_TIME = 10000;   
    int THINNING = 1;
    int BURNIN = 0;
    // Define sizes and temperatures
    int numTemperatures = 32;
    int numSizes = 1;

    int* Sizes = new int[numSizes]{100*100}; 
    double* Temperatures = new double[numTemperatures];
    for (int i = 0; i < numTemperatures; i++) {
        Temperatures[i] = 0.5 + i * 0.125; 
    }

    omp_set_num_threads(16);
    for(int i = 0; i < numSizes; i++){
        int N = *(Sizes+i);
        #pragma omp parallel for
        for (int j = 0; j < numTemperatures; j++){
                double T = *(Temperatures+j);
                // Build the lattice
                double p_init;
                if (T < 2.28){
                    p_init = 0.80;
                } else {
                    p_init = 0.5;
                }
                SpinConfiguration config(N, T, p_init, H);
                // Attach a MCMC engine
                config.mountMCMCengine(MCMCType::MetropolisHastings);
                config.setPBC(true); // Set periodic boundary conditions
                config.keepTrack(true, true); // Magnetization only
                // Now run the MCMC chain
                config.run(MAX_TIME, BURNIN, THINNING);
                // Get the results
                std::vector<double> energies = config.getResults().getEnergy();
                std::vector<double> magnetizations = config.getResults().getMagnetizations();

                // Write to file according to the specified format:
                if (!magnetizations.empty()){
                    std::ofstream outFileMagnetization("./data/magnetization/magnetization_N" + std::to_string(N) + "_T" + std::to_string(T).substr(0,5) + "_H" + std::to_string(H).substr(0,4) + "_t" + std::to_string(MAX_TIME) + "_th" + std::to_string(THINNING) + ".txt");
                    for (const auto& mag : magnetizations) {
                        outFileMagnetization << std::setprecision(10) << mag << "\n"; 
                    }
                    std::cout << "Wrote to file: " << "./data/magnetization/magnetization_N" + std::to_string(N) + "_T" + std::to_string(T).substr(0,5) + "_H" + std::to_string(H).substr(0,4) + "_t" + std::to_string(MAX_TIME) + "_th" + std::to_string(THINNING) + ".txt" << std::endl;
                    outFileMagnetization.close();
                }
                if (!energies.empty()){
                    std::ofstream outFileEnergy("./data/energy/energy_N" + std::to_string(N) + "_T" + std::to_string(T).substr(0,5) + "_H" + std::to_string(H).substr(0,4) + "_t" + std::to_string(MAX_TIME) + "_th" + std::to_string(THINNING) + ".txt");
                    for (const auto& energy : energies) {
                        outFileEnergy << std::setprecision(10) << energy << "\n";
                    }
                    std::cout << "Wrote to file: " << "./data/energy/energy_N" + std::to_string(N) + "_T" + std::to_string(T).substr(0,5) + "_H" + std::to_string(H).substr(0,4) + "_t" + std::to_string(MAX_TIME) + "_th" + std::to_string(THINNING) + ".txt" << std::endl;
                    outFileEnergy.close();
                }

                std::cout << "Simulation completed." << std::endl;
                std::cout << "Average magnetization: " << std::accumulate(magnetizations.begin(), magnetizations.end(), 0.0) / magnetizations.size() << std::endl;
            }
    }
}
