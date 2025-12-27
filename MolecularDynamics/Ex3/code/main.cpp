#include "NDMolDyn.hpp"
#include <iostream>
#include <omp.h>

/*
int main(){
    omp_set_num_threads(1);
    #pragma omp parallel for
    for (int N = 57; N <= 357; N += 25){
        for(int _ = 0; _ < 10; ++_){
            double V = 500.0;
            double density = double(N) / V;

            double r_cut = 2.5;
            double skin = 0.3;
            int M = 2;

            double L = std::pow(static_cast<double>(N) / density, 1.0 / 3.0);
            double ax = int(L/(r_cut+skin));
            std::cout << ax << " > " << M << std::endl;

            int num_steps = 100000;
            double initialTemperature = 2; //Does not guarantee precise T = 2!

            NDMolDyn<3> sim("../fcc32.xyz", r_cut, skin, M, initialTemperature);
            sim.attachAndersenThermostat(initialTemperature, 10);

            std::string filename = "../data/Andersen/energy_omega10_N" + std::to_string(N) + "_run" + std::to_string(_) + ".dat";
            sim.MDRun(num_steps, 0.005, 10, filename);
        }
    }
*/

int main(){
    double r_cut = 3.5;
    double skin = 0.3;
    int M = 1;
    double initialTemperature = 0.0001;

    int num_steps = 100000;
    double L = 4.621845;
    NDMolDyn<3> sim("../fcc108.xyz", L, r_cut, skin, M, initialTemperature);
    sim.attachAndersenThermostat(initialTemperature, 100);

    sim.MDRun(num_steps, 0.005, 10, "dump.txt"); // Devo rifare andersen omega = 100, ricorda
}