#include "ActiveIsingModel.hpp"
#include <iostream>
#include <omp.h>

int main(){
    // Order parameters
    double beta = 1.9;
    int Nx = 200;
    int Ny = 50;

    int totalSteps = 20000;
    int saveEvery = 40;

    double* epsilons = new double[8]{0, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75};
    double* rho = new double[8]{1, 2, 3, 4, 5, 6, 7, 8};    
    omp_set_num_threads(8);
    for (int j = 0; j < 8; j++){
        double r = rho[j];
        std::cout << "Starting simulations for r = " << r << "\n";
        #pragma omp parallel for
        for(int i = 0; i < 8; i++){
            double epsilon = epsilons[i];
            ActiveIsingModel model(Nx, Ny, r, epsilon, 1./beta, initType::UNIFORM, true);
            model.run(totalSteps, saveEvery);
            std::cout << "Done with T = " << 1./beta << "\n";
        }
        std::cout << "Done with r = " << r << "\n";
    }
}

