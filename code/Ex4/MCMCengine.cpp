#include "MCMCengine.hpp"
#include "configuration.hpp"


void MCMCengine::run() {
    if (m_type == MCMCType::MetropolisHastings) {
        this->MetropolisHastings();
    } else if (m_type == MCMCType::GibbsSampling) {
        throw std::runtime_error("Gibbs Sampling not implemented yet");
    } else if (m_type == MCMCType::HeatBath){
        this->HeatBath();
    }
}

void MCMCengine::sweep() {
    for(int i = 0; i < m_N; i++){
        this->run();
    }
}

void MCMCengine::MetropolisHastings() {
    std::uniform_int_distribution<int> dist(0, m_N - 1);
    int index = dist(m_rng);
    int L = static_cast<int>(std::sqrt(m_N));
    int row = index / L;
    int col = index % L;

    // extract the spin value
    int spin = m_spinConfig->getSpins()[index] ? 1 : -1;
    // Get the neighbourhood spins (1 matrix)
    int neighbourSum = 0;
    std::array<int,4> neighbours = m_spinConfig->getNeighbourhoodSpins(row, col);
    for (int i = 0; i < 4; ++i) {
        neighbourSum += neighbours[i];
    }
    double exponentialFactor = std::exp(- 2* double(spin) * static_cast<double>(neighbourSum) / m_spinConfig->getTemperature() + 2 * m_spinConfig->getMagneticField() * double(spin) / m_spinConfig->getTemperature());
    // Metropolis-Hastings criterion: accept if ...
    if (exponentialFactor > 1.0 || std::uniform_real_distribution<double>(0.0, 1.0)(m_rng) < exponentialFactor) {
        m_spinConfig->getSpins()[index] = !m_spinConfig->getSpins()[index];   // modify the pointed spins
    }
}


void MCMCengine::HeatBath() {
    std::uniform_int_distribution<int> dist(0, m_N - 1);
    int index = dist(m_rng);
    int L = static_cast<int>(std::sqrt(m_N));
    int row = index / L;
    int col = index % L;

    // Get the neighbourhood spins (1 matrix)
    int neighbourSum = 0;
    int spin = m_spinConfig->getSpins()[index] ? 1 : -1;
    std::array<int,4> neighbours = m_spinConfig->getNeighbourhoodSpins(row, col);
    for (int i = 0; i < 4; ++i) {
        neighbourSum += neighbours[i];
    }
    
    // Calculate the probability of being in state +1 or -1
    double exponentialFactor = std::exp(- 2* double(spin) * static_cast<double>(neighbourSum) / m_spinConfig->getTemperature() + 2 * m_spinConfig->getMagneticField() * double(spin) / m_spinConfig->getTemperature());
    double pFlip = 1. / (1.+(1./exponentialFactor));

    if (std::uniform_real_distribution<double>(0.0, 1.0)(m_rng) < pFlip) {
        m_spinConfig->getSpins()[index] = !m_spinConfig->getSpins()[index];  
    } 
}