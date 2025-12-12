#include "ActiveIsingModel.hpp"
#include <iostream>

int main(int argc, char* argv[]){
    if (argc != 6){
        std::cout << "Usage: ./sample.out <beta> <rho> <epsilon> <totalSteps> <saveEvery>\n";
        return -1;
    }
    // Read parameters from command line
    double beta = std::stod(argv[1]);
    double rho = std::stod(argv[2]);
    double epsilon = std::stod(argv[3]);
    int totalSteps = std::stoi(argv[4]);
    int saveEvery = std::stoi(argv[5]);
    // Geometric constraints
    int Nx = 200;
    int Ny = 50;
    // Instantiate the model ...
    ActiveIsingModel model(Nx, Ny, rho, epsilon, 1./beta, initType::UNIFORM, true);
    // ... and run it
    model.run(totalSteps, saveEvery);
}

