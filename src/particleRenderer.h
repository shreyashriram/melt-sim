// particleRenderer.h
#ifndef PARTICLE_RENDERER_H
#define PARTICLE_RENDERER_H

#include <vector>
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <string>
#include "particleSampler.h"

class ParticleRenderer {
public:
    ParticleRenderer();
    ~ParticleRenderer();
    
    void init(const std::string& shaderPath);
    void render(const std::vector<Particle>& particles, const glm::mat4& view, const glm::mat4& projection);
    
private:
    unsigned int shaderProgram;
    unsigned int VAO, VBO;
    
    void setupParticles(const std::vector<Particle>& particles);
};

#endif // PARTICLE_RENDERER_H