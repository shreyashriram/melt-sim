#pragma once

#include <vector>
#include <glm/glm.hpp>
#include "particle.h"
#include "grid.h"
#include <Eigen/Dense>

using scalar = float;
using Vector3 = Eigen::Matrix<scalar, 3, 1>;

class MPMSimulation
{
public:
    Grid grid;    
    std::vector<Particle> particles;

    MPMSimulation();
    void addMeshParticles(std::vector<Vector3> sampledPoints);
    void initializeParticles();
    void step(float dt);
    void transferParticlesToGrid();
    void transferGridToParticles(float dt);
    void updateParticles(float dt);
    void computeForces();
    void updateGrid(float dt);

    float youngsModulus;
    float poissonsRatio;

    float computeWeight(const glm::vec3 &particlePos, const glm::vec3 &nodePos);
    glm::vec3 computeWeightGradient(const glm::vec3 &particlePos, const glm::vec3 &nodePos);
};