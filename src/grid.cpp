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
    
    /*std::cout << "===== Grid Initialization Debug =====" << std::endl;
    std::cout << "Grid size: " << gridSize << "x" << gridSize << "x" << gridSize 
              << " (" << grid.size() << " total nodes)" << std::endl;
    std::cout << "Cell spacing: " << gridSpacing << std::endl;
    std::cout << "Grid physical dimensions: " 
              << gridSize * gridSpacing << "x" 
              << gridSize * gridSpacing << "x"
              << gridSize * gridSpacing << std::endl;
              
    // Print a few sample grid nodes
    if (grid.size() > 0) {
        std::cout << "Sample node (0,0,0): mass=" << grid[0].mass 
                  << ", vel=(" << grid[0].velocity.x << "," 
                  << grid[0].velocity.y << "," 
                  << grid[0].velocity.z << ")" << std::endl;*/
    }

void Grid::updateGrid() { //STEP 2 IN THE MPM GUIDE
    //EVENTUALLY THIS WILL HAVE INTERNAL FORCES FROM STRESS, BOUNDARY CONDITIONS, DAMPING 

    //apply gravity
    for (auto& node : grid) {
        if (node.mass > 0.0f) {
            node.force += glm::vec3(0.0f, -9.8f, 0.0f) * node.mass;
        }
    }
}