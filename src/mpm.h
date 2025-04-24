#pragma once

#include <vector>
#include <glm/glm.hpp>
#include "particle.h"
#include "grid.h"

class MPMSimulation {
public:
    MPMSimulation(std::vector<Particle>& particles);
    void initialize();
    void step(float dt);
    void transferParticlesToGrid();
    void transferGridToParticles();
    void updateParticles(float dt, int gridSize, float gridSpacing);
    void computeForces();
    void updateGrid();

private:
    float youngsModulus;
    float poissonsRatio;
    int gridSize;
    float gridSpacing;
    std::vector<Particle>& particles;  // Reference to particles from main
    std::vector<GridNode> grid;

    float computeWeight(const glm::vec3& particlePos, const glm::vec3& nodePos);
    glm::vec3 computeWeightGradient(const glm::vec3& particlePos, const glm::vec3& nodePos);
}; 