#include "particle_renderer.h"
#include <glm/gtc/type_ptr.hpp>

ParticleRenderer::ParticleRenderer() : maxParticles(1000) {
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, maxParticles * sizeof(glm::vec2), nullptr, GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), (void*)0);
    glEnableVertexAttribArray(0);
}

void ParticleRenderer::update(const std::vector<Particle>& particles) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, particles.size() * sizeof(glm::vec2), particles.data());
}

void ParticleRenderer::draw() {
    glBindVertexArray(vao);
    glDrawArrays(GL_POINTS, 0, (GLsizei)maxParticles);
}