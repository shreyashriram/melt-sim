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
    const std::vector<GridNode>& getGrid() const { return grid; }
    int gridSize;
    float gridSpacing;

    private:
    std::vector<GridNode> grid;
    void initializeGrid();
};


