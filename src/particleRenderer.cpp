// particleRenderer.cpp
#include "particleRenderer.h"
#include "shader_utils.h"
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

ParticleRenderer::ParticleRenderer() : shaderProgram(0), VAO(0), VBO(0) {}

ParticleRenderer::~ParticleRenderer() {
    if (VAO != 0) {
        glDeleteVertexArrays(1, &VAO);
    }
    if (VBO != 0) {
        glDeleteBuffers(1, &VBO);
    }
}

void ParticleRenderer::init(const std::string& shaderPath) {
    // Assuming you have vertex and fragment shaders for particles
    std::string vertexPath = shaderPath + "/particle_vertex.glsl";
    std::string fragmentPath = shaderPath + "/particle_fragment.glsl";
    
    shaderProgram = createShaderProgram(vertexPath, fragmentPath);
}

void ParticleRenderer::setupParticles(const std::vector<Particle>& particles) {
    if (VAO != 0) {
        glDeleteVertexArrays(1, &VAO);
    }
    if (VBO != 0) {
        glDeleteBuffers(1, &VBO);
    }
    
    std::vector<float> particleData;
    particleData.reserve(particles.size() * 3); // Only position data for now
    
    for (const auto& particle : particles) {
        particleData.push_back(particle.position.x);
        particleData.push_back(particle.position.y);
        particleData.push_back(particle.position.z);
    }
    
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    
    glBindVertexArray(VAO);
    
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, particleData.size() * sizeof(float), particleData.data(), GL_STATIC_DRAW);
    
    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    glBindVertexArray(0);
}

void ParticleRenderer::render(const std::vector<Particle>& particles, const glm::mat4& view, const glm::mat4& projection) {
    if (particles.empty()) {
        return;
    }
    
    setupParticles(particles);
    
    glUseProgram(shaderProgram);
    
    // Pass matrices to shader
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    
    // Particle color
    glUniform3f(glGetUniformLocation(shaderProgram, "particleColor"), 0.2f, 0.6f, 1.0f);
    
    // Point size (you might want to adjust this based on screen resolution)
    glPointSize(3.0f);
    
    glBindVertexArray(VAO);
    glDrawArrays(GL_POINTS, 0, particles.size());
    glBindVertexArray(0);
}