#include "grid.h"

#include <glad/glad.h>
#include <glm/gtc/type_ptr.hpp>

#include <glad/glad.h>
#include <glm/gtc/type_ptr.hpp>


Grid::Grid( int size, float spacing ) : size(size), spacing(spacing) {
    // Sets the size of the "flattened" grid to gridSize x gridSize x gridSize
    // ex. (5 x 5 x 5)
    nodes.resize(size * size * size);

    for (GridNode& n : nodes) { //initialize each node in the grid to zero

        n.velocity = glm::vec3(0.0f);
        n.force = glm::vec3(0.0f);
        n.mass = 0.0f;

    }
    
}
void Grid::setupBuffers() {
    if (initialized) return;

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    std::vector<glm::vec3> lines;

    float start = 0.0f;
    float end = size * spacing;

    for (int i = 0; i <= size; ++i) {
        float offset = i * spacing;

        // Y-Z planes (lines along X)
        for (int j = 0; j <= size; ++j) {
            float y = j * spacing;
            lines.emplace_back(start, y, offset);
            lines.emplace_back(end,   y, offset);
        }

        // X-Z planes (lines along Y)
        for (int j = 0; j <= size; ++j) {
            float x = j * spacing;
            lines.emplace_back(x, start, offset);
            lines.emplace_back(x, end,   offset);
        }

        // X-Y planes (lines along Z)
        for (int j = 0; j <= size; ++j) {
            float x = j * spacing;
            lines.emplace_back(x, offset, start);
            lines.emplace_back(x, offset, end);
        }
    }

    glBufferData(GL_ARRAY_BUFFER, lines.size() * sizeof(glm::vec3), lines.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);

    initialized = true;
}

void Grid::draw() const {
    if (!initialized) return;
    glBindVertexArray(VAO);
    glDrawArrays(GL_LINES, 0, (size + 1) * (size + 1) * 6); // each line has 2 vertices, 3 sets per offset
    glBindVertexArray(0);
}