// File: src/mpmSimulator.h
#ifndef MPM_SIMULATOR_H
#define MPM_SIMULATOR_H

#include <vector>
#include <glm/glm.hpp>
#include "particleSampler.h"

// Grid cell structure for the MPM grid
struct GridCell {
    glm::vec3 velocity;
    float mass;
    glm::vec3 force;
    
    GridCell() : velocity(0.0f), mass(0.0f), force(0.0f) {}
    void reset() {
        velocity = glm::vec3(0.0f);
        mass = 0.0f;
        force = glm::vec3(0.0f);
    }
};

class MPMSimulator {
public:
    MPMSimulator(float gridSize = 0.1f, float timeStep = 0.01f);
    ~MPMSimulator();
    
    // Initialize the simulator with particles
    void initialize(std::vector<Particle>& particles);
    
    // Update the simulation by one time step
    void update(std::vector<Particle>& particles);
    
    // Set simulation parameters
    void setGravity(const glm::vec3& gravity) { m_gravity = gravity; }
    void setDensity(float density) { m_density = density; }
    void setViscosity(float viscosity) { m_viscosity = viscosity; }
    void setTimeStep(float timeStep) { m_timeStep = timeStep; }
    void setRestDensity(float restDensity) { m_restDensity = restDensity; }
    void setBoundaries(const glm::vec3& min, const glm::vec3& max) {
        m_boundaryMin = min;
        m_boundaryMax = max;
    }
    
private:
    // Grid resolution and dimensions
    float m_gridSize;
    int m_gridResolutionX;
    int m_gridResolutionY;
    int m_gridResolutionZ;
    std::vector<GridCell> m_grid;
    
    // Simulation parameters
    float m_timeStep;
    glm::vec3 m_gravity;
    float m_density;
    float m_viscosity;
    float m_restDensity;
    float m_stiffness;  // For equation of state
    
    // Simulation boundaries
    glm::vec3 m_boundaryMin;
    glm::vec3 m_boundaryMax;
    
    // Helper methods
    void resetGrid();
    void particlesToGrid(const std::vector<Particle>& particles);
    void updateGridVelocities();
    void gridToParticles(std::vector<Particle>& particles);
    void enforceBoundaries(std::vector<Particle>& particles);
    
    // Grid indexing helpers
    int gridIndex(int i, int j, int k) const;
    bool isValidGridIndex(int i, int j, int k) const;
    glm::ivec3 particleToGridPos(const glm::vec3& position) const;
    
    // Interpolation weights
    float weight(float x) const;
    float weightDerivative(float x) const;
};

#endif // MPM_SIMULATOR_H