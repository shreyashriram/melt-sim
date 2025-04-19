#include "floor.h"
#include <glad/glad.h>

void setupFloor(unsigned int& VAO, unsigned int& VBO, unsigned int& EBO) {
    float floorVertices[] = {
        // positions         // normals (facing up)
        -5.0f, -1.0f, -5.0f,  0.0f, 1.0f, 0.0f,
         5.0f, -1.0f, -5.0f,  0.0f, 1.0f, 0.0f,
         5.0f, -1.0f,  5.0f,  0.0f, 1.0f, 0.0f,
        -5.0f, -1.0f,  5.0f,  0.0f, 1.0f, 0.0f
    };

    unsigned int floorIndices[] = {
        0, 1, 2,
        2, 3, 0
    };

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(floorVertices), floorVertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(floorIndices), floorIndices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);
}