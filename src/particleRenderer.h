#pragma once
#include "particle.h"
#include <vector>
#include <glad/glad.h>

class ParticleRenderer {
public:
    ParticleRenderer();
    ~ParticleRenderer();
    void init(const std::vector<Particle>& particles);
    void update(const std::vector<Particle>& particles);
    void draw();

private:
    unsigned int VAO, VBO;
    size_t numParticles;
};
