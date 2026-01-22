#include "NDMolDyn.hpp"

int main(){
    Parameters params;
    params.readFromFile("../parameter.yaml");
    NDMolDyn<2> sim(params);
    sim.MCRun(10, 10, "../data/traj.xyz", "../data/energy.dat");
    //acceptance rate: solid 0.12 \pm 0.06
    //acceptance rate: liquid 0.38 \pm 0.05
    //acceptance rate: gas 0.92 \pm 0.02
}