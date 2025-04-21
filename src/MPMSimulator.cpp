#include "MPMSimulator.h"
#include <cmath>

MPMSimulator::MPMSimulator() : dt(0.003f), gravity(-9.8f) {
    initParticles();
}

void MPMSimulator::initParticles() {
    // Simple grid of particles
    int countX = 20, countY = 20;
    for (int i = 0; i < countX; ++i) {
        for (int j = 0; j < countY; ++j) {
            glm::vec2 pos = glm::vec2(i * 0.05f, j * 0.05f + 1.0f); // slightly above ground
            particles.emplace_back(pos);
        }
    }
}

void MPMSimulator::step() {
    for (auto& p : particles) {
        p.acc = glm::vec2(0.0f, gravity);
        p.vel += p.acc * dt;
        p.pos += p.vel * dt;

        // Floor collision
        if (p.pos.y < 0.0f) {
            p.pos.y = 0.0f;
            p.vel.y *= -0.3f;
        }
    }
}

std::vector<Particle>& MPMSimulator::getParticles() {
    return particles;
}