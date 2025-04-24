#include "Particle.h"

void Particle::update(float dt) {
    velocity += glm::vec3(0.0f, -9.8f, 0.0f) * dt;
    position += velocity * dt;
    if (position.y < 0.0f) {
        position.y = 0.0f;
        velocity.y *= -0.5f;
    }
}
