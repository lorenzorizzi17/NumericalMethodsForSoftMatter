#ifndef RESULTS_HPP
#define RESULTS_HPP

#include <vector>
#include <numeric>

class Results {
    private:
    std::vector<double> m_magnetizations;
    std::vector<double> m_energy;

    public:

    Results() = default;

    void storeMagnetization(double mag) {
        m_magnetizations.push_back(mag);
    }

    void storeEnergy(double energy) {
        m_energy.push_back(energy);
    }

    std::vector<double>& getMagnetizations() {
        return m_magnetizations;
    }
    
    std::vector<double>& getEnergy() {
        return m_energy;
    }
};

#endif