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
    void updateGrid();
    const std::vector<GridNode>& getGrid() const { return grid; }

    private:
    int gridSize;
    float gridSpacing;
    std::vector<GridNode> grid;
    void initializeGrid();
}


