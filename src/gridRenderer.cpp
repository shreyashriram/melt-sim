#include "gridRenderer.h"
#include <glad/glad.h>
#include <glm/gtc/type_ptr.hpp>

GridRenderer::GridRenderer(int gridSize, float gridSpacing)
    : gridSize(gridSize), gridSpacing(gridSpacing) {}

void GridRenderer::init() {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &massVBO);
    
    glBindVertexArray(VAO);
    
    // Position buffer
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    // reserve maximum space
    glBufferData(GL_ARRAY_BUFFER,
                 gridSize*gridSize*gridSize * sizeof(glm::vec3),
                 nullptr,
                 GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    
    // Mass buffer
    glBindBuffer(GL_ARRAY_BUFFER, massVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 gridSize*gridSize*gridSize * sizeof(float),
                 nullptr,
                 GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, (void*)0);
    
    glBindVertexArray(0);
}

void GridRenderer::update(const std::vector<GridNode>& nodes) {
    positions.clear();
    masses.clear();
    positions.reserve(nodes.size());
    masses.reserve(nodes.size());
    
    glm::vec3 offset = glm::vec3(gridSize*gridSpacing*0.5f);

    // only include nodes with nonzero mass
    for (int z = 0; z < gridSize; ++z)
    for (int y = 0; y < gridSize; ++y)
    for (int x = 0; x < gridSize; ++x) {
        int idx = x + y*gridSize + z*gridSize*gridSize;
        if (nodes[idx].mass > 0.0f) {
            glm::vec3 pos = glm::vec3(x, y, z)*gridSpacing - offset;
            positions.push_back(pos);
            masses.push_back(nodes[idx].mass);
        }
    }

    // Update position buffer
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferSubData(GL_ARRAY_BUFFER,
                    0,
                    positions.size()*sizeof(glm::vec3),
                    positions.data());
                    
    // Update mass buffer
    glBindBuffer(GL_ARRAY_BUFFER, massVBO);
    glBufferSubData(GL_ARRAY_BUFFER,
                    0,
                    masses.size()*sizeof(float),
                    masses.data());
}

void GridRenderer::draw(const glm::mat4& model,
                        const glm::mat4& view,
                        const glm::mat4& projection,
                        unsigned int shaderProgram) {
    glUseProgram(shaderProgram);
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"),  1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    
    // We no longer set a fixed color since the fragment shader will color by mass
    glPointSize(25.0f);

    glBindVertexArray(VAO);
    glDrawArrays(GL_POINTS, 0, (GLsizei)positions.size());
    glBindVertexArray(0);
}