#ifndef MOL_DYN_3D_HPP
#define MOL_DYN_3D_HPP

#include <vector>
#include <iostream>
#include <random>
#include "Particle.hpp"
#include "Parameter.hpp"

#define multidx std::array<int, Dim>

template<int Dim>
class NDMolDyn{
    private:
        std::vector<Particle<Dim>> m_particles;      // The vector holding all particles
        std::vector<std::vector<int>> m_verletLists; // verlet lists for each particle
        std::vector<int> m_headOfChain;              // head of chain for cell list
        std::vector<vecd> m_backupParticles; // backup of particles
        
        void buildVerletLists();             // Build the verlet lists
        double dist_sq_pbc(const vecd&, const vecd&) const;
        double computeAngle(Patch const&, vecd const&, vecd const&, double) const;

        // Internal counters
        int m_acceptCounter = 0;
        int m_numVerletCalls = 0;
        // Random number generators
        std::mt19937 m_gen; 
        std::uniform_real_distribution<double> m_unif_dist;
    public:
        double m_L;       // Box side length, measured in sigma units ( L = 1 means box of size 1 diameter)
        double m_V;
        double m_rho;
        double m_temperature;
        
        int m_numMCSteps;
        double m_rcut;
        double m_skin;
        int m_M;
        double max_displacement;
        int m_idx;
        
        double m_cosThetaMax;
        double m_gamma;
        double m_energyWell;

        NDMolDyn(int const, double const, double const, double const, int const, double const, double const, int, double);
        NDMolDyn(Parameters const& params, int seed);
        
        void MCRun(int, int, std::string const&, std::string const&);
        void MCRun(int, int, std::ofstream&, std::ofstream&, int);

        void singleMCStep();
        double getEnergy(vecd const&,Patch const&, int) const;
        double getEnergyNoPatches(vecd const&,Patch const&, int) const;
        double getTotalEnergy() const;

        std::vector<vecd> getPositions() const {
            std::vector<vecd> positions;
            positions.reserve(m_particles.size());
            for(const auto& p : m_particles) {
                positions.push_back(p.position);
            }
            return positions;
        }
        void setTemperature(double T) {
            m_temperature = T;
        }
        void setNumMCSteps(int nSteps) {
            m_numMCSteps = nSteps;
        }
        void setIndex(int idx) {
            m_idx = idx;
        }
};

#include "NDMolDyn.tpp"
#endif