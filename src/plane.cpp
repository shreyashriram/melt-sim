#include "plane.h"
#include <glm/gtc/matrix_transform.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glad/glad.h>
#include <vector>


Plane::Plane(const glm::vec3& center, float rotation_degrees, float size):
    center(center), rotation_degrees(rotation_degrees), size(size)
{
    generatePlane();
    setupBuffers();  
}

void Plane::generatePlane() {

    float halfSize = size/2.0f;

    std::vector<glm::vec3> temp = {

        {-halfSize, 0.0f,  halfSize}, // 0
        {-halfSize, 0.0f, -halfSize}, // 1 
        { halfSize, 0.0f, -halfSize}, // 2
        { halfSize, 0.0f,  halfSize}, // 3
        
    };

    // vertices array 
    float radians = glm::radians(rotation_degrees);

    for (auto& v : temp ) {
        v = glm::rotateZ(v, radians); // rotate along z axis "roll"
        v += center;
    
    };
    vertices = temp;

    // normals array 
    glm::vec3 normal = glm::rotateZ(glm::vec3(0, 1, 0), glm::radians(rotation_degrees));

    for (int i = 0; i < 4; ++i) {
        normals.push_back(glm::vec3(0, 1, 0));
    }

    // indices array 
    indices = { 
        0, 1, 2,
        0, 2, 3 
    };
}

void Plane::setupBuffers() {
    // Generate VAO, VBO, EBO
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    glGenBuffers(1, &NBO);

    glBindVertexArray(VAO);

    glBindVertexArray(VAO);

    // Positions
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), vertices.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);

    // Normals
    glBindBuffer(GL_ARRAY_BUFFER, NBO);
    glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(glm::vec3), normals.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);

    // Indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

    glBindVertexArray(0);
}

void Plane::draw() const {
    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(indices.size()), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}
