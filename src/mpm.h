#pragma once

#include <vector>
#include <glm/glm.hpp>
#include "particle.h"
#include "grid.h"
#include <Eigen/Dense>

using scalar = float;
using Vector3 = Eigen::Matrix<scalar, 3, 1>;

class MPMSimulation {
public:
    std::vector<Particle> particles;

    MPMSimulation();
    void addMeshParticles(std::vector<Vector3> sampledPoints);
    void initializeParticles();
    void step(float dt);
    void transferParticlesToGrid();
    void transferGridToParticles();
    void updateParticles(float dt, int gridSize, float gridSpacing);
    void computeForces();
    void updateGrid();

    float youngsModulus;
    float poissonsRatio;
    // int gridSize;
    // float gridSpacing;  
    
    Grid grid;
    // std::vector<GridNode> grid;

    float computeWeight(const glm::vec3& particlePos, const glm::vec3& nodePos);
    glm::vec3 computeWeightGradient(const glm::vec3& particlePos, const glm::vec3& nodePos);
}; 