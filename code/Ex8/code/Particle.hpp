#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <array>

#define multidx std::array<int, Dim>
#define vecd std::array<double, Dim>

template<int Dim>
struct Particle{
    vecd position;
    vecd velocity;

    int next = -1;   // -1 means null

    Particle() = default;
    Particle(const vecd& pos) : position(pos) {}
};


#endif