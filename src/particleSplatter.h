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

    // Set parameters for the fluid rendering
    void setMetaballThreshold(float threshold) { metaballThreshold = threshold; }
    void setMetaballStrength(float strength) { metaballStrength = strength; }
    void setDropletScale(float scale) { dropletScale = scale; }
    void setDropletIntensity(float intensity) { dropletIntensity = intensity; }
    void enableWaterDroplets(bool enable) { waterDropletsEnabled = enable; }

private:
    unsigned int VAO, VBO, instanceVBO;
    unsigned int shaderProgram;
    
    std::vector<glm::vec4> particleData; // Position and velocity data interleaved
    size_t numParticles;
    
    // Splatter parameters
    float radius;
    float smoothing;
    
    // Metaball parameters
    float metaballThreshold;
    float metaballStrength;
    
    // Water droplet parameters
    float dropletScale;
    float dropletIntensity;
    bool waterDropletsEnabled;
    
    // Texture for water droplets
    unsigned int waterTexture;
    bool waterTextureLoaded;
    
    // Helper method to load or generate water texture
    unsigned int createWaterDropletTexture();
};