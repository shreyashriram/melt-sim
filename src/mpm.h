
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
    Grid grid;    
    std::vector<Particle> particles;

    MPMSimulation();
    void addMeshParticles(std::vector<Vector3> sampledPoints, MaterialType type);
    void spawnCube(MaterialType type, glm::vec3 center, float spacing, int countPerAxis);

    void step(float dt);
    void transferParticlesToGrid();
    void transferGridToParticles(float dt);
    void updateParticles(float dt);
    void computeForces();
    void updateGrid(float dt);
    
    void resolveParticleCollisions(float radius, float stiffness);


    float youngsModulus;
    float poissonsRatio;
    
    // Lamé parameters (derived from Young's modulus and Poisson's ratio)
    float shearModulus;       // Controls resistance to shape change (μ)
    float bulkModulus;        // Controls resistance to volume change (λ)

    // Plasticity parameters
    float yieldThreshold;     // When to transition from elastic to plastic behavior
    float meltRate;           // How quickly the material melts
    float globalMeltProgress; // Tracks overall melting progress


    float computeWeight(const glm::vec3& particlePos, const glm::vec3& nodePos);
    glm::vec3 computeWeightGradient(const glm::vec3& particlePos, const glm::vec3& nodePos);

    void polarDecomposition(const glm::mat3& F, glm::mat3& R, glm::mat3& S);
    glm::mat3 computeVelocityGradient(const Particle& p);
    void updatePlasticity(Particle& p, float dt);
    glm::mat3 computeStress(const Particle& p);
    void testComputeStress();
    void runTests();



    void updateSolidParticle(Particle& p, float dt);
    void updateLiquidParticle(Particle& p, float dt);
    void updateMeltingParticle(Particle& p, float dt);

};
