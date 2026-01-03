#include "NDMolDyn.hpp"

int main(){
    Parameters params;
    params.readFromFile("../parameter.yaml");
    NDMolDyn<2> sim(params, 1545);
    sim.MCRun(1, 1, "../data/traj.xyz", "../data/energy.dat");
}

// L = 16