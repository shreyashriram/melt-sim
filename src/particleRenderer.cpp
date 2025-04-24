#include "particleRenderer.h"

ParticleRenderer::ParticleRenderer() : VAO(0), VBO(0), numParticles(0) {}

ParticleRenderer::~ParticleRenderer() {
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
}

void ParticleRenderer::init(const std::vector<Particle>& particles) {
    numParticles = particles.size();

    std::vector<float> data;
    for (const auto& p : particles) {
        data.push_back(p.position.x);
        data.push_back(p.position.y);
        data.push_back(p.position.z);
    }

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(float), data.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

void ParticleRenderer::update(const std::vector<Particle>& particles) {
    std::vector<float> data;
    for (const auto& p : particles) {
        data.push_back(p.position.x);
        data.push_back(p.position.y);
        data.push_back(p.position.z);
    }
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, data.size() * sizeof(float), data.data());
}

void ParticleRenderer::draw() {
    glBindVertexArray(VAO);
    glPointSize(10.0f);
    glDrawArrays(GL_POINTS, 0, numParticles);
}
