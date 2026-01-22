#include <omp.h>
#include <array>
#include "NDMolDyn.hpp"

int main(){
    omp_set_num_threads(8);
    std::random_device rd;
    unsigned int master_seed = rd();
    std::cout << "Master seed: " << master_seed << std::endl;
    std::mt19937 gen(master_seed + 12345);
    std::uniform_int_distribution<int> dist(0, 7);
    std::array<double, 8> beta = {3,4,5,6,7,8,9,10};  // beta_c = 1.666
    std::array<NDMolDyn<2>*, 8> thermalMap; //Each element of this array will point to the correct chain with the corresponding temperature
    std::vector<NDMolDyn<2>> chains;  // Backup chains for swap attempts
    std::vector<int> swaps(8,0);
    chains.reserve(8);
    // I/O Files
    std::vector<std::ofstream> traj_files(8);
    std::vector<std::ofstream> ene_files(8);
    /*
Total swap acceptance ratio: 0.307
Swap counts per chain: 
Chain 0 (beta=1.05): 42 swaps.
Chain 1 (beta=1.25): 95 swaps.
Chain 2 (beta=1.45): 88 swaps.
Chain 3 (beta=1.65): 62 swaps.
Chain 4 (beta=1.85): 66 swaps.
Chain 5 (beta=2.05): 91 swaps.
Chain 6 (beta=2.25): 111 swaps.
Chain 7 (beta=2.45): 59 swaps.


Total swap acceptance ratio: 0.438
Swap counts per chain: 
Chain 0 (beta=11): 101 swaps.
Chain 1 (beta=10): 190 swaps.
Chain 2 (beta=9): 182 swaps.
Chain 3 (beta=8): 158 swaps.
Chain 4 (beta=7): 118 swaps.
Chain 5 (beta=6): 85 swaps.
Chain 6 (beta=5): 37 swaps.
Chain 7 (beta=4): 5 swaps.

Total swap acceptance ratio: 0.245
Swap counts per chain: 
Chain 0 (beta=11.5): 32 swaps.
Chain 1 (beta=10.5): 73 swaps.
Chain 2 (beta=9.5): 75 swaps.
Chain 3 (beta=8.5): 95 swaps.
Chain 4 (beta=7.5): 112 swaps.
Chain 5 (beta=6.5): 72 swaps.
Chain 6 (beta=5.5): 26 swaps.
Chain 7 (beta=4.5): 5 swaps.

    */



    int N_TOT = 1000000;
    int N_SWAP = 1000;
    std::cout << "Starting Parallel Tempering with " << N_TOT << " total MC steps and swaps every " << N_SWAP << " steps. Max num of swaps:" << N_TOT/N_SWAP << std::endl;
    // Create K MC chains serially (it's ok, it's just initialization)
    for(int i = 0; i < 8; ++i) {
        Parameters params;
        params.readFromFile("../parameter.yaml");
        NDMolDyn<2> sim(params, 5267*i + master_seed);  // Seed the RNG differently for each thread
        sim.setTemperature(1./beta[i]);                 // Set the initial temperatures
        sim.setNumMCSteps(N_SWAP);                      // Set number of MC steps per block
        sim.setIndex(i);
        chains.push_back(sim);                          // Store a backup copy of the chain
        thermalMap[i] = &chains.back();                 // Store the pointer to the chain

        traj_files[i].open("../data/dump" + std::to_string(beta[i]) + "._noSwap.xyz");
        ene_files[i].open("../data/energy_dump" + std::to_string(beta[i]) + "._noSwap.dat");
        //traj_files[i].open("../data/traj_beta" + std::to_string(beta[i]) + ".xyz");
        //ene_files[i].open("../data/energy_beta" + std::to_string(beta[i]) + ".dat");
    }
    int num_swaps = 0;
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        for(int block = 0; block < N_TOT / N_SWAP; ++block) {
            
            NDMolDyn<2>* myChain = thermalMap[thread_id];
            myChain->MCRun(10, 10, traj_files[thread_id], ene_files[thread_id], block * N_SWAP);
            #pragma omp barrier
            #pragma omp master
            {
                std::uniform_int_distribution<int> dist_idx(0, 8 - 2); // Indici da 0 a 6
                int k = dist_idx(gen); // Tentiamo scambio tra k e k+1

                NDMolDyn<2>* chain_cold = thermalMap[k];
                NDMolDyn<2>* chain_hot  = thermalMap[k+1];

                double beta_cold = beta[k];
                double beta_hot  = beta[k+1];
                
                double E_cold = chain_cold->getTotalEnergy(); 
                double E_hot  = chain_hot->getTotalEnergy();

                double delta = (beta_cold - beta_hot) * (E_cold - E_hot);

                std::uniform_real_distribution<double> dist_prob(0.0, 1.0);
                if (delta > 0 || std::exp(delta) > dist_prob(gen)) {
                    chain_cold->setTemperature(1./beta[k+1]);
                    chain_hot->setTemperature(1./beta[k]);
                    std::swap(thermalMap[k], thermalMap[k+1]);

                    std::cout << "Block " << block << ": Swap accepted " << k << " <-> " << k+1 << std::endl;
                    std::cout << std::endl;
                    num_swaps++;
                    swaps[k]++;
                    swaps[k+1]++;
                }
            } 
            #pragma omp barrier
        }
    } 
    std::cout << "Total swap acceptance ratio: " << static_cast<double>(num_swaps) / (N_TOT / N_SWAP) << std::endl;
    std::cout << "Swap counts per chain: " << std::endl;
    for(int i = 0; i < 8; ++i) {
        std::cout << "Chain " << i << " (beta=" << beta[i] << "): " << swaps[i] << " swaps." << std::endl;
    }
    for(int i = 0; i < 8; ++i) {
        traj_files[i].close();
        ene_files[i].close(); 
    }
}