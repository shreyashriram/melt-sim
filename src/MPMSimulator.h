#pragma once
#include <vector>
#include "Particle.h"

class MPMSimulator {
public:
    MPMSimulator();
    void initParticles();
    void step();
    std::vector<Particle>& getParticles();

private:
    std::vector<Particle> particles;
    float dt;
    float gravity;
};