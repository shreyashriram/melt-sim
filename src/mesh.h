#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <glad/glad.h>

class Mesh {
public:
    Mesh(const std::string& path);
    ~Mesh();

    void draw() const;
    void bind() const;
    void unBind() const;

    std::vector<float> vertices;
    std::vector<unsigned int> indices;

    // Future extensibility:
    // glm::mat4 transform;
    // std::string name;
    // material info, texture, etc.

private:
    
    unsigned int VAO, VBO, EBO;
    void setupBuffers();
    
};

#endif