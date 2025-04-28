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
    void init(float particleRadius = 0.1f, float smoothingKernel = 0.8f);
    
    // Update the particles data
    void update(const std::vector<Particle>& particles);
    
    // Draw the particles as fluid-like splats
    void draw(const glm::mat4& model, const glm::mat4& view, const glm::mat4& projection);

    // Set parameters for the water rendering
    void setReflectionStrength(float strength) { reflectionStrength = strength; }
    void setSpecularPower(float power) { specularPower = power; }
    void setSpecularStrength(float strength) { specularStrength = strength; }
    void setFresnelBias(float bias) { fresnelBias = bias; }
    void setFresnelScale(float scale) { fresnelScale = scale; }
    void setFresnelPower(float power) { fresnelPower = power; }
    void setRippleStrength(float strength) { rippleStrength = strength; }
    void setRippleSpeed(float speed) { rippleSpeed = speed; }
    void setWaterColor(const glm::vec3& color) { waterColor = color; }
    void setSurfaceTensionStrength(float strength) { surfaceTensionStrength = strength; }
    void setBlendDistance(float distance) { blendDistance = distance; }
    void enableReflections(bool enable) { reflectionsEnabled = enable; }
    void enableRefractions(bool enable) { refractionsEnabled = enable; }

private:
    unsigned int VAO, VBO, instanceVBO;
    unsigned int shaderProgram;
    
    std::vector<glm::vec4> particleData; // Position and velocity data interleaved
    size_t numParticles;
    
    // Surface parameters
    float radius;
    float smoothing;
    
    // Liquid blending parameters
    float surfaceTensionStrength;
    float blendDistance;
    
    // Water appearance parameters
    glm::vec3 waterColor;
    float reflectionStrength;
    float specularPower;
    float specularStrength;
    float fresnelBias;
    float fresnelScale;
    float fresnelPower;
    float rippleStrength;
    float rippleSpeed;
    bool reflectionsEnabled;
    bool refractionsEnabled;
    
    // Textures for water effects
    unsigned int normalMapTexture;
    unsigned int environmentMapTexture;
    bool texturesLoaded;
    
    // Helper methods to load or generate textures
    unsigned int createNormalMapTexture();
    unsigned int createEnvironmentMapTexture();
    
    // Helper method to update time
    float currentTime;
    void updateTime();
};