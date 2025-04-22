// particleSampler.h
#ifndef PARTICLE_SAMPLER_H
#define PARTICLE_SAMPLER_H

#include <vector>
#include <string>
#include <glm/glm.hpp>

struct Particle {
    glm::vec3 position;
    glm::vec3 velocity;
    float mass;
    float temperature;
    bool isFixed;
    
    Particle(const glm::vec3& pos, float m = 1.0f, float temp = 20.0f) 
        : position(pos), velocity(glm::vec3(0.0f)), mass(m), temperature(temp), isFixed(false) {}
};

class ParticleSampler {
public:
    ParticleSampler();
    ~ParticleSampler();
    
    // Create particles from mesh - sampling a volume
    std::vector<Particle> sampleMeshVolume(const std::string& meshPath, float particleRadius);
    
    // Create particles from just the surface
    std::vector<Particle> sampleMeshSurface(const std::string& meshPath, float particleSpacing);
    
    // Utility to visualize particles
    unsigned int createParticleVAO(const std::vector<Particle>& particles);

private:
    // Helper methods
    bool isPointInMesh(const glm::vec3& point, const std::vector<glm::vec3>& vertices, 
                       const std::vector<unsigned int>& indices);
    std::vector<glm::vec3> extractVertices(const std::string& meshPath);
    std::vector<unsigned int> extractIndices(const std::string& meshPath);
    glm::vec3 calculateBoundingBoxMin(const std::vector<glm::vec3>& vertices);
    glm::vec3 calculateBoundingBoxMax(const std::vector<glm::vec3>& vertices);
};

#endif // PARTICLE_SAMPLER_H