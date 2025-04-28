#pragma once
#include <vector>
#include <glm/glm.hpp>

struct GridNode {
    
    glm::vec3 velocity;
    glm::vec3 force;
    float mass;
    float fluidFraction;

    GridNode() : velocity(0.0f), force(0.0f), mass(0.0f), fluidFraction(0.0f) {}
};

class Grid {
public: 
    Grid(int size, float spacing);
    
    int size;
    float spacing;
    
    std::vector<GridNode> nodes;

    void setupBuffers(); 
    void draw() const;          

private:
    unsigned int VAO, VBO;
    
    bool initialized = false;

};