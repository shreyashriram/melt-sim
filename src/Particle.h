#pragma once
#include <glm/glm.hpp>

struct Particle {
    glm::vec2 pos;
    glm::vec2 vel;
    glm::vec2 acc;
    float mass;
    float density;
    float pressure;

    Particle(glm::vec2 p) : pos(p), vel(0.0f), acc(0.0f), mass(1.0f), density(0.0f), pressure(0.0f) {}
};