#ifndef PLANE_H
#define PLANE_H

#include <glm/glm.hpp>
#include <vector>

class Plane {
public:
    Plane( const glm::vec3& center, float rotation_degrees = 0.0f, float size = 1.0f );
    
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> normals;
    std::vector<unsigned int> indices;


    void draw() const;

private:

    glm::vec3 center;
    float rotation_degrees;
    float size;

    unsigned int VAO, VBO, NBO, EBO;

    

    void setupBuffers();
    void generatePlane();
};

#endif