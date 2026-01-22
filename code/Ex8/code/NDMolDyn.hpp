#ifndef MOL_DYN_3D_HPP
#define MOL_DYN_3D_HPP

#include <vector>
#include <iostream>
#include "Particle.hpp"

#define multidx std::array<int, Dim>
#define vecd std::array<double, Dim>


enum class ThermostatType {
    NONE,
    BERENDSEN,
    ANDERSEN
};

template<int Dim>
class NDMolDyn{
    private:
        std::vector<Particle<Dim>> m_particles;      // The vector holding all particles
        std::vector<std::vector<int>> m_verletLists; // verlet lists for each particle
        std::vector<int> m_headOfChain;              // head of chain for cell list
        std::vector<vecd> m_backupParticles; // backup of particles

        std::vector<vecd> m_forces;               // Forces on each particle, cache
        
        void buildVerletLists();             // Build the verlet lists
        double dist_sq_pbc(const vecd&, const vecd&) const;

        // Internal counters
        int m_numVerletCalls = 0;
        ThermostatType m_Thermostat = ThermostatType::NONE;

    public:
        double m_L;       // Box side length, measured in sigma units ( L = 1 means box of size 1 diameter)
        double m_V;
        double m_rho;
        double m_initialTemperature;
        
        int m_numMDSteps;
        double m_rcut;
        double m_skin;
        int m_M;

        int m_tau_berendsen;
        double m_BerendsenTemperature;
        double m_omega_Andersen;

        NDMolDyn(int const, double const, double const, double const, int const, double const);
        NDMolDyn(std::string const&, double const, double const, double const, int const, double const);
        
        void MDRun(int, double, int, std::string const&);

        void computeForces();

        double computeKinEnergy() const;
        double computePotEnergy() const;

        void attachBerendsenThermostat(double targetTemp, int tau){
            m_Thermostat = ThermostatType::BERENDSEN;
            m_BerendsenTemperature = targetTemp;
            m_tau_berendsen = tau;
        }

        void attachAndersenThermostat(double targetTemp, double omega){
            m_Thermostat = ThermostatType::ANDERSEN;
            m_omega_Andersen = omega;
            m_BerendsenTemperature = targetTemp; // Reuse variable
        }
        
        std::vector<vecd> getPositions() const {
            std::vector<vecd> positions;
            positions.reserve(m_particles.size());
            for(const auto& p : m_particles) {
                positions.push_back(p.position);
            }
            return positions;
        }

        std::vector<vecd> getVelocities() const {
            std::vector<vecd> velocities;
            velocities.reserve(m_particles.size());
            for(const auto& p : m_particles) {
                velocities.push_back(p.velocity);
            }
            return velocities;
        }
};

#include "NDMolDyn.tpp"
#endif