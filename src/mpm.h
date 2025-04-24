#pragma once

#include <vector>
#include <glm/glm.hpp>


class MPMSimulation {
public:
    MPMSimulation();
    void initialize();
    void step(float dt);
    void transferParticlesToGrid();
    void transferGridToParticles();
    void computeForces();

private:
    float youngsModulus;  //idk which one this should go in
    float poissonsRatio;

    float computeWeight(const glm::vec3& particlePos, const glm::vec3& nodePos);
    glm::vec3 computeWeightGradient(const glm::vec3& particlePos, const glm::vec3& nodePos);
}; 