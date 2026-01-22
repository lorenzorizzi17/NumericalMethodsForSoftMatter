#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <array>
#include "Patches.hpp"

#define multidx std::array<int, Dim>

template<int Dim>
struct Particle{
    vecd position;
    int next = -1;   // -1 means null
    Patch patch;     // For now, single patch per particle
    Particle() = default;
    Particle(vecd pos = {}): position(pos) {patch.setRandomPatch();}
    Patch getPatch() const {return patch;}
};


#endif