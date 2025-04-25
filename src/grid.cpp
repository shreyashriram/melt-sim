#include "grid.h"

void Grid::initializeGrid() {
    int gridSize = 4; //4 cells in each dimension
    float cellSize = 0.25f; //width of each cell (h) 
    
    this->gridSize = gridSize; //set the grid properties
    this->gridSpacing = cellSize;
    
    grid.resize(gridSize * gridSize * gridSize);
    
    for (GridNode& node : grid) { //initialize each node in the grid to zero 
        node.velocity = glm::vec3(0.0f);
        node.force = glm::vec3(0.0f);
        node.mass = 0.0f;
    }
    
    }
