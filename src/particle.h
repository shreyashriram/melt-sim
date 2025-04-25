#pragma once
#include <vector>
#include <glm/glm.hpp>

struct Particle {
    glm::vec3 position;
    glm::vec3 velocity;
    glm::vec3 force;

    float mass;
    float J;
    
    glm::mat3 F;
    glm::mat3 C;

    Particle() : position(0.0f), velocity(0.0f), force(0.0f), mass(1.0f), J(1.0f), F(1.0f), C(0.0f) {};
    Particle(const glm::vec3& pos) : position(pos), velocity(0.0f), force(0.0f), mass(1.0f), J(1.0f), F(1.0f), C(0.0f) {};
    Particle(const glm::vec3& pos, const glm::vec3& vel) : position(pos), velocity(vel), force(0.0f), mass(1.0f), J(1.0f), F(1.0f), C(0.0f) {};
};

