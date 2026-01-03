#include "NDMolDyn.hpp"
#include <omp.h>

int main(){
    double* cosThetaMax = new double[2]{0.5, 0.93937271285}; 
    double* Temperatures = new double[8]{0.25, 0.26, 0.28, 0.30, 0.35, 0.40, 0.45, 0.5}; // gas, liquid, solid 
    omp_set_num_threads(8);
    for(int i = 1; i < 2; ++i){
        #pragma omp parallel for // to speed up the simulations
        for(int j = 0; j < 8; ++j){
            NDMolDyn<3> sim(200, 0.05, 1, 0.5, 7, 0.4, Temperatures[j], 200000, cosThetaMax[i]);
            sim.MCRun(true, true);
        }
    }
}