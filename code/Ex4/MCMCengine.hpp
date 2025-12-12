#ifndef MCMCENGINE_HPP
#define MCMCENGINE_HPP

#include <random>
#include <stdexcept>


enum class MCMCType {
    MetropolisHastings,
    GibbsSampling,
    HeatBath
};
class SpinConfiguration; // Forward declaration

class MCMCengine {
    private:
    SpinConfiguration* m_spinConfig; // Pointer to SpinConfiguration
    std::mt19937 m_rng; // Random number generator
    MCMCType m_type;  // Type of MCMC engine
    int m_N;
    public:
    // Constructor
    MCMCengine(SpinConfiguration* spinConfig, MCMCType type, int N) : m_spinConfig(spinConfig), m_type(type), m_N(N) {
        // initialize the random number engine
        std::random_device rd;
        m_rng = std::mt19937(rd());
    }
    //default constructor
    MCMCengine() = default;

    void run();
    void sweep();
    //helper methods
    void MetropolisHastings();
    void GibbsSampling() {
        throw std::runtime_error("Gibbs Sampling not implemented yet");
    }
    void HeatBath();
};
















#endif