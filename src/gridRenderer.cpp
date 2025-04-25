#include "gridRenderer.h"
#include <glad/glad.h>
#include <glm/gtc/type_ptr.hpp>

GridRenderer::GridRenderer(int gridSize, float gridSpacing)
    : gridSize(gridSize), gridSpacing(gridSpacing), velocityScale(0.5f) {}

void GridRenderer::init() {
    // Initialize grid nodes rendering
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &massVBO);
    
    glBindVertexArray(VAO);
    
    // Position buffer
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
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
    
    // Initialize velocity arrows rendering
    glGenVertexArrays(1, &arrowVAO);
    glGenBuffers(1, &arrowVBO);
    
    glBindVertexArray(arrowVAO);
    glBindBuffer(GL_ARRAY_BUFFER, arrowVBO);
    // Reserve space for arrow vertices (will update dynamically)
    glBufferData(GL_ARRAY_BUFFER,
                 gridSize*gridSize*gridSize * 6 * sizeof(glm::vec3), // For each node: main line + arrow tips
                 nullptr,
                 GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    
    glBindVertexArray(0);
}

void GridRenderer::update(const std::vector<GridNode>& nodes) {
    positions.clear();
    masses.clear();
    arrowVertices.clear();
    
    positions.reserve(nodes.size());
    masses.reserve(nodes.size());
    
    glm::vec3 offset = glm::vec3(gridSize*gridSpacing*0.5f);

    // Only include nodes with nonzero mass
    for (int z = 0; z < gridSize; ++z) {
        for (int y = 0; y < gridSize; ++y) {
            for (int x = 0; x < gridSize; ++x) {
                int idx = x + y*gridSize + z*gridSize*gridSize;
                if (nodes[idx].mass > 0.0f) {
                    glm::vec3 pos = glm::vec3(x, y, z)*gridSpacing - offset;
                    positions.push_back(pos);
                    masses.push_back(nodes[idx].mass);
                    
                    glm::vec3 velocity = nodes[idx].velocity;
                    float velocityMagnitude = glm::length(velocity);
                    
                    if (velocityMagnitude > 0.01f) {  // Only draw arrows for non-negligible velocities
                        glm::vec3 normalizedVelocity = glm::normalize(velocity);
                        glm::vec3 arrowEnd = pos + velocity * velocityScale;
                        
                        // Main arrow line (start and end points)
                        arrowVertices.push_back(pos);
                        arrowVertices.push_back(arrowEnd);
                        
                        // Create arrow tip (two additional lines)
                        glm::vec3 perpendicular;
                        
                        // If velocity is nearly parallel to (0,1,0), use (1,0,0) as reference
                        if (abs(normalizedVelocity.y) > 0.99f) {
                            perpendicular = glm::normalize(glm::cross(normalizedVelocity, glm::vec3(1.0f, 0.0f, 0.0f)));
                        } else {
                            perpendicular = glm::normalize(glm::cross(normalizedVelocity, glm::vec3(0.0f, 1.0f, 0.0f)));
                        }
                        
                        // Arrow tip size - proportional to velocity magnitude but clamped
                        float tipSize = glm::min(velocityMagnitude * 0.2f, gridSpacing * 0.25f);
                        
                        // Calculate arrow tip points
                        glm::vec3 tip1 = arrowEnd - normalizedVelocity * tipSize * 0.8f + perpendicular * tipSize;
                        glm::vec3 tip2 = arrowEnd - normalizedVelocity * tipSize * 0.8f - perpendicular * tipSize;
                        
                        // Add first side of arrow tip
                        arrowVertices.push_back(arrowEnd);
                        arrowVertices.push_back(tip1);
                        
                        // Add second side of arrow tip
                        arrowVertices.push_back(arrowEnd);
                        arrowVertices.push_back(tip2);
                    }
                }
            }
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
                    
    // Update velocity arrows buffer
    glBindBuffer(GL_ARRAY_BUFFER, arrowVBO);
    glBufferSubData(GL_ARRAY_BUFFER,
                    0,
                    arrowVertices.size()*sizeof(glm::vec3),
                    arrowVertices.data());
}

void GridRenderer::draw(const glm::mat4& model,
                        const glm::mat4& view,
                        const glm::mat4& projection,
                        unsigned int shaderProgram) {
    glUseProgram(shaderProgram);
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"),  1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    
    // Draw grid nodes
    glUniform1f(glGetUniformLocation(shaderProgram, "isArrow"), 0.0f); // Not an arrow
    glPointSize(12.0f);
    glBindVertexArray(VAO);
    glDrawArrays(GL_POINTS, 0, (GLsizei)positions.size());
    
    // Draw velocity arrows
    glUniform1f(glGetUniformLocation(shaderProgram, "isArrow"), 1.0f); // This is an arrow
    glUniform3f(glGetUniformLocation(shaderProgram, "arrowColor"), 1.0f, 1.0f, 0.0f); // Yellow arrows
    
    glBindVertexArray(arrowVAO);
    glLineWidth(2.0f); // Set line width for better visibility
    glDrawArrays(GL_LINES, 0, (GLsizei)arrowVertices.size());
    
    glBindVertexArray(0);
}