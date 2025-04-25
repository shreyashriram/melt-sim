#include "grid.h"


Grid::Grid( int size, float spacing ) : size(size), spacing(spacing) {
    // Sets the size of the "flattened" grid to gridSize x gridSize x gridSize
    // ex. (5 x 5 x 5)
    nodes.resize(size * size * size);


    
    for (GridNode& n : nodes) { //initialize each node in the grid to zero

        n.velocity = glm::vec3(0.0f);
        n.force = glm::vec3(0.0f);
        n.mass = 0.0f;

    }
    
}
