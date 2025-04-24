#pragma once
#include <glm/glm.hpp>

struct Particle {
    glm::vec3 position;
    glm::vec3 velocity;

    Particle(glm::vec3 pos) : position(pos), velocity(0.0f) {}
    void update(float dt);
};
