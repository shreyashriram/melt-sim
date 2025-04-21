#pragma once
#include "Particle.h"
#include <vector>
#include <glad/glad.h>

class ParticleRenderer {
public:
    ParticleRenderer();
    void update(const std::vector<Particle>& particles);
    void draw();

private:
    GLuint vao, vbo;
    size_t maxParticles;
};