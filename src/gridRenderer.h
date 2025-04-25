#pragma once
#include <vector>
#include <glm/glm.hpp>
#include "grid.h"

class GridRenderer {
public:
    GridRenderer(int gridSize, float gridSpacing);
    void init();
    void update(const std::vector<GridNode>& nodes);
    void draw(const glm::mat4& model,
              const glm::mat4& view,
              const glm::mat4& projection,
              unsigned int shaderProgram);

private:
    int    gridSize;
    float  gridSpacing;
    uint VAO, VBO;
    std::vector<glm::vec3> positions;
};
