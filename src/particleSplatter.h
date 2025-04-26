#pragma once
#include "particle.h"
#include <vector>
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <algorithm> // For std::sort

class ParticleSplatter {
public:
    ParticleSplatter();
    ~ParticleSplatter();

    // Initialize the splatter with a particle radius and smoothing parameter
    void init(float particleRadius = 0.05f, float smoothingKernel = 0.8f);
    
    // Update the particles data
    void update(const std::vector<Particle>& particles);
    
    // Draw the particles as fluid-like splats
    void draw(const glm::mat4& model, const glm::mat4& view, const glm::mat4& projection);

private:
    unsigned int VAO, VBO, instanceVBO;
    unsigned int shaderProgram;
    
    std::vector<glm::vec4> particleData; // Position (xyz) and size (w)
    size_t numParticles;
    
    float radius;
    float smoothing;
    
    // Helper method to create embedded shader program (fallback)
    unsigned int createEmbeddedShaderProgram();
};