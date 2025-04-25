#pragma once
#include <vector>
#include <glm/glm.hpp>

struct GridNode {
    
    glm::vec3 velocity;
    glm::vec3 force;
    float mass;

};

class Grid {
public: 
    Grid(int size, float spacing);
    
    int size;
    float spacing;
    
    std::vector<GridNode> nodes;

};


