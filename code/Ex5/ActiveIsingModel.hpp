#include <vector>
#include <numeric>
#include <random>
#include "SFML/Graphics.hpp"

// A simple wrapper for a vector of ints, i.e. spins in a cell
struct Cells {
    private:
        std::vector<int> m_spins;
    public:
        Cells(int N_SPINS) {m_spins.reserve(2*N_SPINS);}
        int getMagnetization() const {
            int sum = std::accumulate(m_spins.begin(), m_spins.end(), 0);
            return sum;
        }
        int getDensity() const {
            return (m_spins.size());
        }
        std::vector<int>& getData() { return m_spins; }
        int& operator[](int index){ return m_spins[index]; }
        void addSpin(int spin){ m_spins.push_back(spin); }
};

enum initType { UNIFORM, PEAKED };

// Main class
class ActiveIsingModel {
    private:
        std::mt19937 m_gen;
        std::vector<Cells> m_cells;
        std::vector<int> m_fluxes;
        bool m_saveFlux = false;
        int N_CELL_X;
        int N_CELL_Y;
        double m_epsilon = 0.9; // Activity parameter
        double m_deltaT;
        double m_T = 1.0; // Temperature
        double m_rho;
        int N_SPINS;
        double m_time = 0.0;
        double m_D = 1;

        void MCMC_SpinFlip(int, int);
        void MCMC_Hop(int, int);
        
    public:
        //Standard constructor, initializes the system with random spins and uniform distribution
        ActiveIsingModel(int, int, double, double, double, initType = UNIFORM, bool = false);
        // A getter for the raw data
        const std::vector<Cells>& getSpins() const;
        // Enable or disable flux saving
        void saveFlux(bool var) {m_saveFlux = var;};
        // Run and launch and SFML window
        void runGraphics();
        // Run the simulation and save data to file
        void run(int, int);
        void MCMCStep();
};  