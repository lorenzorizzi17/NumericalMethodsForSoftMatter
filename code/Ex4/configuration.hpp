#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP

#include <random>
#include <vector>
#include <array>
#include <SFML/Graphics.hpp>

#include"MCMCengine.hpp"
#include "results.hpp"


class SpinConfiguration {
    private:
    // Physical parameters
    int m_N;
    double m_T;
    double m_h = 0.0; // Magnetic field, default to 0.0
    // Data structure to hold the spins
    std::vector<bool> m_spins;

    // RNG and MCMC engine
    std::mt19937 m_rng;
    MCMCengine m_engine;


    bool m_PBC = true; // Periodic Boundary Conditions, true by default (otherwise open boundary conditions)
    bool m_keepTrackof[2];  // Magnetization / energy
    Results m_results; // Results object to store the results of the simulation


    public:
    // Standard constructor
    SpinConfiguration(int, double, double);
    SpinConfiguration(int, double, double, double);

    // getters
    double getMagnetization() const;
    double getEnergy() const;
    std::vector<bool>& getSpins() { return m_spins; }
    double getTemperature() const { return m_T; }
    double getMagneticField() const { return m_h; }
    void setMagneticField(double h) { m_h = h; }
    Results& getResults() { return m_results; }

    // MCMC methods (single one spin flip or N spins)
    void mountMCMCengine(MCMCType type) {
        m_engine = MCMCengine(this, type, m_N);
    }
    void setPBC(bool pbc) { m_PBC = pbc; }
    void keepTrack(bool magnetization, bool energy) {
        m_keepTrackof[0] = magnetization;
        m_keepTrackof[1] = energy;
    }
    
    void step();
    void sweep();
    std::array<int,4> getNeighbourhoodSpins(int, int) const;

    // Graphics / simulations
    void run(int, int, int);
    void runGraphics(int, int);

};

#endif 