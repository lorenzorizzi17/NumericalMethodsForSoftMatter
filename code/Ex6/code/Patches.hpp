#ifndef PATCHES_HPP
#define PATCHES_HPP

#include<cmath>
#include<array>

#define vecd std::array<double, 3>

struct Patch {
public:
    vecd direction; 
    Patch() = default;
    Patch(const vecd& dir) {
        double norm = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
        direction[0] = dir[0] / norm;
        direction[1] = dir[1] / norm;
        direction[2] = dir[2] / norm;
    }

    void setRandomPatch() {
        double z = (static_cast<double>(rand()) / RAND_MAX) * 2.0 - 1.0;
        double phi = (static_cast<double>(rand()) / RAND_MAX) * 2.0 * M_PI;
        double r_xy = std::sqrt(1.0 - z*z); 
        direction[0] = r_xy * std::cos(phi);
        direction[1] = r_xy * std::sin(phi);
        direction[2] = z;
    }

    double getX() const { return direction[0]; };
    double getY() const { return direction[1]; };
    double getZ() const { return direction[2]; };
};











#endif