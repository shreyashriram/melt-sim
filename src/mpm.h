#pragma once

#include <vector>
#include <glm/glm.hpp>
#include "particle.h"
#include "grid.h"

class MPMSimulation {
public:
    MPMSimulation();
    void initialize();
    void step(float dt);
    void transferParticlesToGrid();
    void transferGridToParticles();
    void computeForces();
    void updateGrid();

private:
    float youngsModulus;
    float poissonsRatio;
    int gridSize;
    float gridSpacing;
    std::vector<Particle> particles;
    std::vector<GridNode> grid;

    float computeWeight(const glm::vec3& particlePos, const glm::vec3& nodePos);
    glm::vec3 computeWeightGradient(const glm::vec3& particlePos, const glm::vec3& nodePos);
}; 