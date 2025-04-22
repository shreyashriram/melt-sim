#ifndef MPM_SIMULATION_H
#define MPM_SIMULATION_H

#include <vector>
#include <glm/glm.hpp>
#include <glad/glad.h>
#include "particleSampler.h" // Include the existing particle definition

struct MaterialPoint {
    glm::vec3 position;
    glm::vec3 velocity;
    glm::mat3 deformationGradient = glm::mat3(1.0f);
    glm::mat3 affineVelocity = glm::mat3(0.0f);
    float temperature = 0.0f;
    float fluidity = 0.0f;         // 0.0 = solid, 1.0 = fluid
    float plasticVolume = 1.0f;    // For plasticity tracking
    glm::vec3 color;               // Visualization color
};

struct GridNode {
    glm::vec3 velocity = glm::vec3(0.0f);
    float mass = 0.0f;
};

class MPMSimulation {
public:
    MPMSimulation();
    ~MPMSimulation();
    
    void initialize(const std::vector<Particle>& initialParticles);
    void update(float deltaTime);
    void render(const glm::mat4& view, const glm::mat4& projection);
    
private:
    // Simulation parameters
    static const int gridSize = 64;             // Number of grid cells in each direction
    const float dx;                             // Grid cell size
    const float invDx;                          // 1/dx
    const float particleMass = 1.0f;            // Mass of each particle
    const float particleVolume = 1.0f;          // Volume of each particle
    const float youngsModulus = 1e5f;           // Elastic stiffness
    const float poissonRatio = 0.2f;            // Poisson's ratio
    const float yieldStress = 1e3f;             // Stress threshold for melting
    const float decayRate = 0.1f;               // Rate at which material becomes fluid
    float mu;                                   // Lame parameter (shear modulus)
    float lambda;                               // Lame parameter (bulk modulus)
    
    // Grid and particles
    std::vector<GridNode> grid;                 // Background Eulerian grid
    std::vector<MaterialPoint> particles;       // Lagrangian material points
    
    // Simulation steps
    void resetGrid();
    void particleToGrid(float dt);              // Updated to include dt parameter
    void updateGrid(float dt);
    void gridToParticle(float dt);
    void updateMaterialState(MaterialPoint& p, float dt);
    
    // Rendering
    unsigned int particleVAO;
    unsigned int particleVBO;
    unsigned int particleShaderProgram;
};

#endif // MPM_SIMULATION_H