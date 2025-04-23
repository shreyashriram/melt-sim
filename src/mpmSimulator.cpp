// File: src/mpmSimulator.cpp
#include "mpmSimulator.h"
#include <iostream>
#include <algorithm>
#include <cmath>

MPMSimulator::MPMSimulator(float gridSize, float timeStep)
    : m_gridSize(gridSize), m_timeStep(timeStep),
      m_gravity(0.0f, -9.81f, 0.0f),
      m_density(1000.0f), m_viscosity(0.1f),
      m_restDensity(1000.0f), m_stiffness(50.0f),
      m_boundaryMin(-2.0f, 0.05f, -2.0f), m_boundaryMax(2.0f, 4.0f, 2.0f) {
    
    // Initialize grid dimensions based on boundary
    m_gridResolutionX = static_cast<int>((m_boundaryMax.x - m_boundaryMin.x) / m_gridSize) + 1;
    m_gridResolutionY = static_cast<int>((m_boundaryMax.y - m_boundaryMin.y) / m_gridSize) + 1;
    m_gridResolutionZ = static_cast<int>((m_boundaryMax.z - m_boundaryMin.z) / m_gridSize) + 1;
    
    // Allocate grid cells
    int totalCells = m_gridResolutionX * m_gridResolutionY * m_gridResolutionZ;
    m_grid.resize(totalCells);
    
    std::cout << "MPM Grid initialized with dimensions: " 
              << m_gridResolutionX << " x " 
              << m_gridResolutionY << " x " 
              << m_gridResolutionZ << " (" 
              << totalCells << " cells)" << std::endl;
}

MPMSimulator::~MPMSimulator() {
    // Nothing to clean up
}

void MPMSimulator::initialize(std::vector<Particle>& particles) {
    // Initialize all particles with default properties
    for (auto& particle : particles) {
        particle.velocity = glm::vec3(0.0f);
        particle.mass = 1.0f;
        particle.temperature = 20.0f;
        particle.isFixed = false;
    }
    
    std::cout << "MPM Simulator initialized with " << particles.size() << " particles" << std::endl;
}

void MPMSimulator::update(std::vector<Particle>& particles) {
    // Reset the grid
    resetGrid();
    
    // Step 1: Transfer particle data to grid (P2G)
    particlesToGrid(particles);
    
    // Step 2: Update grid velocities (apply forces)
    updateGridVelocities();
    
    // Step 3: Transfer grid data back to particles (G2P)
    gridToParticles(particles);
    
    // Step 4: Update particle positions
    for (auto& particle : particles) {
        if (!particle.isFixed) {
            particle.position += particle.velocity * m_timeStep;
        }
    }
    
    // Step 5: Apply boundary conditions
    enforceBoundaries(particles);
}

void MPMSimulator::resetGrid() {
    for (auto& cell : m_grid) {
        cell.reset();
    }
}

void MPMSimulator::particlesToGrid(const std::vector<Particle>& particles) {
    // For each particle, distribute mass and momentum to surrounding grid nodes
    for (const auto& particle : particles) {
        // Get the grid cell containing the particle
        glm::ivec3 cellIdx = particleToGridPos(particle.position);
        
        // Loop over 3x3x3 neighborhood of cells around the particle
        for (int di = -1; di <= 1; di++) {
            for (int dj = -1; dj <= 1; dj++) {
                for (int dk = -1; dk <= 1; dk++) {
                    int i = cellIdx.x + di;
                    int j = cellIdx.y + dj;
                    int k = cellIdx.z + dk;
                    
                    if (!isValidGridIndex(i, j, k)) continue;
                    
                    // Calculate distance to grid node
                    glm::vec3 gridPos = m_boundaryMin + glm::vec3(i, j, k) * m_gridSize;
                    glm::vec3 dist = (particle.position - gridPos) / m_gridSize;
                    
                    // Calculate interpolation weight
                    float w = weight(dist.x) * weight(dist.y) * weight(dist.z);
                    
                    // Accumulate mass and momentum
                    int idx = gridIndex(i, j, k);
                    m_grid[idx].mass += w * particle.mass;
                    m_grid[idx].velocity += w * particle.mass * particle.velocity;
                }
            }
        }
    }
    
    // Normalize velocities
    for (auto& cell : m_grid) {
        if (cell.mass > 0.0f) {
            cell.velocity /= cell.mass;
        }
    }
}

void MPMSimulator::updateGridVelocities() {
    // Apply forces to grid velocities
    for (auto& cell : m_grid) {
        if (cell.mass > 0.0f) {
            // Apply gravity
            cell.velocity += m_gravity * m_timeStep;
            
            // Apply other forces (pressure, viscosity, etc.)
            cell.velocity += cell.force * m_timeStep / cell.mass;
        }
    }
}

void MPMSimulator::gridToParticles(std::vector<Particle>& particles) {
    // Update particle velocities based on interpolated grid velocities
    for (auto& particle : particles) {
        if (particle.isFixed) continue;
        
        // Get the grid cell containing the particle
        glm::ivec3 cellIdx = particleToGridPos(particle.position);
        
        // Reset particle velocity
        glm::vec3 newVelocity(0.0f);
        
        // Loop over 3x3x3 neighborhood of cells around the particle
        for (int di = -1; di <= 1; di++) {
            for (int dj = -1; dj <= 1; dj++) {
                for (int dk = -1; dk <= 1; dk++) {
                    int i = cellIdx.x + di;
                    int j = cellIdx.y + dj;
                    int k = cellIdx.z + dk;
                    
                    if (!isValidGridIndex(i, j, k)) continue;
                    
                    // Calculate distance to grid node
                    glm::vec3 gridPos = m_boundaryMin + glm::vec3(i, j, k) * m_gridSize;
                    glm::vec3 dist = (particle.position - gridPos) / m_gridSize;
                    
                    // Calculate interpolation weight
                    float w = weight(dist.x) * weight(dist.y) * weight(dist.z);
                    
                    // Interpolate grid velocity to particle
                    int idx = gridIndex(i, j, k);
                    if (m_grid[idx].mass > 0.0f) {
                        newVelocity += w * m_grid[idx].velocity;
                    }
                }
            }
        }
        
        // Update particle velocity
        particle.velocity = newVelocity;
    }
}

void MPMSimulator::enforceBoundaries(std::vector<Particle>& particles) {
    for (auto& particle : particles) {
        if (particle.isFixed) continue;
        
        // Apply boundary conditions (bounce off walls)
        const float damping = 0.5f;  // Velocity damping factor on collision
        
        // X boundaries
        if (particle.position.x < m_boundaryMin.x) {
            particle.position.x = m_boundaryMin.x;
            particle.velocity.x *= -damping;
        } else if (particle.position.x > m_boundaryMax.x) {
            particle.position.x = m_boundaryMax.x;
            particle.velocity.x *= -damping;
        }
        
        // Y boundaries (floor and ceiling)
        if (particle.position.y < m_boundaryMin.y) {
            particle.position.y = m_boundaryMin.y;
            particle.velocity.y *= -damping;
        } else if (particle.position.y > m_boundaryMax.y) {
            particle.position.y = m_boundaryMax.y;
            particle.velocity.y *= -damping;
        }
        
        // Z boundaries
        if (particle.position.z < m_boundaryMin.z) {
            particle.position.z = m_boundaryMin.z;
            particle.velocity.z *= -damping;
        } else if (particle.position.z > m_boundaryMax.z) {
            particle.position.z = m_boundaryMax.z;
            particle.velocity.z *= -damping;
        }
    }
}

int MPMSimulator::gridIndex(int i, int j, int k) const {
    return i + j * m_gridResolutionX + k * m_gridResolutionX * m_gridResolutionY;
}

bool MPMSimulator::isValidGridIndex(int i, int j, int k) const {
    return (i >= 0 && i < m_gridResolutionX &&
            j >= 0 && j < m_gridResolutionY &&
            k >= 0 && k < m_gridResolutionZ);
}

glm::ivec3 MPMSimulator::particleToGridPos(const glm::vec3& position) const {
    glm::ivec3 idx;
    idx.x = static_cast<int>((position.x - m_boundaryMin.x) / m_gridSize);
    idx.y = static_cast<int>((position.y - m_boundaryMin.y) / m_gridSize);
    idx.z = static_cast<int>((position.z - m_boundaryMin.z) / m_gridSize);
    
    // Clamp to valid grid indices
    idx.x = std::max(0, std::min(idx.x, m_gridResolutionX - 1));
    idx.y = std::max(0, std::min(idx.y, m_gridResolutionY - 1));
    idx.z = std::max(0, std::min(idx.z, m_gridResolutionZ - 1));
    
    return idx;
}

// Quadratic B-spline weight function
float MPMSimulator::weight(float x) const {
    x = std::abs(x);
    float result = 0.0f;
    
    if (x < 0.5f) {
        result = 0.75f - x * x;
    } else if (x < 1.5f) {
        float t = 1.5f - x;
        result = 0.5f * t * t;
    }
    
    return result;
}

// Derivative of the weight function
float MPMSimulator::weightDerivative(float x) const {
    float sign = (x >= 0.0f) ? 1.0f : -1.0f;
    x = std::abs(x);
    float result = 0.0f;
    
    if (x < 0.5f) {
        result = -2.0f * x;
    } else if (x < 1.5f) {
        result = -1.0f * (1.5f - x);
    }
    
    return sign * result;
}