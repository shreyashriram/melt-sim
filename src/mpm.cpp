#include "mpm.h"
#include "particle.h"
#include "grid.h"
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>


bool alreadyPrinted = false; //for debugging 

MPMSimulation::MPMSimulation() 
    : gridSize(4), gridSpacing(0.25f), 
      youngsModulus(1.0e6f), poissonsRatio(0.3f) {
}

void MPMSimulation::initialize() {
    // Initialize particles vector
    particles.resize(100); // or whatever initial count you need
    
    // Initialize grid vector
    grid.resize(gridSize * gridSize * gridSize);
    
    // Initialize grid nodes
    for (auto& node : grid) {
        node.mass = 0.0f;
        node.velocity = glm::vec3(0.0f);
        node.force = glm::vec3(0.0f);
    }
}


void MPMSimulation::step(float dt) {
    //reset grid
    for (auto& node : grid) {
        node.velocity = glm::vec3(0.0f);
        node.force = glm::vec3(0.0f);
        node.mass = 0.0f;
    }
    
    transferParticlesToGrid();
    /*if (alreadyPrinted == false){
    std::cout << "===== Grid After Particle Transfer =====" << std::endl;
    for (int z = 0; z < gridSize; ++z) {
        for (int y = 0; y < gridSize; ++y) {
            for (int x = 0; x < gridSize; ++x) {
                int linearIdx = x + y * gridSize + z * gridSize * gridSize;
                GridNode& node = grid[linearIdx];
                //if (node.mass > 0.0f) {  // Only print nodes with mass
                    std::cout << "Node (" << x << "," << y << "," << z << "): "
                              << "mass=" << node.mass 
                              << ", vel=(" << node.velocity.x << "," 
                              << node.velocity.y << "," 
                              << node.velocity.z << ")" << std::endl;
                //}
            }
        }
    }
    std::cout << "===================================" << std::endl;
   // alreadyPrinted = true;
        }*/
    
    updateGrid();
    transferGridToParticles();
    Particles particlesClass;
    particlesClass.updateParticles(1.0f / 60.0f, gridSize, gridSpacing);
}

void MPMSimulation::transferParticlesToGrid() { //STEP 1 IN MPM GUIDE
    for (const auto& p : particles) {
        glm::vec3 cellIdx = (p.position / gridSpacing) - 0.5f; 
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx)); //cell we will build the neighborhood around 
        
        for (int i = 0; i < 3; ++i) { //bc we can't do 1.5 in an int loop, still loop through 3 cells and just adjust the index to be the correct neighborhood
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    glm::ivec3 nodeIdx = baseNode + glm::ivec3(i-1, j-1, k-1); //current node in the neighborhood
                    
                    if (nodeIdx.x < 0 || nodeIdx.x >= gridSize || nodeIdx.y < 0 || nodeIdx.y >= gridSize || nodeIdx.z < 0 || nodeIdx.z >= gridSize) {
                        continue; //if not in the grid bounds, skip it 
                    }
                    
                    int linearIdx = nodeIdx.x + nodeIdx.y * gridSize + nodeIdx.z * gridSize * gridSize; //finding the index in the array
                    glm::vec3 nodePos = glm::vec3(nodeIdx) * gridSpacing; //converting from grid space to world space 
                    
                    float weight = computeWeight(p.position, nodePos); //where we will implement B-spline

                    grid[linearIdx].mass += p.mass * weight; //setting the mass and velocity in the array 
                    grid[linearIdx].velocity += p.mass * p.velocity * weight;
                    //WE STILL NEED TO DO EITHER APIC OR MLS FOR PRESERVED MOMENTUM (INSTEAD OF LINEAR VELOCITY LIKE WE ARE RIGHT NOW)
                }
            }
        }
    }
    
    //have to normalize velocity to convert accumulated momentum back to velocity
    for (auto& node : grid) {
        if (node.mass > 0.0f) {
            node.velocity /= node.mass;
        }
    }
}

void MPMSimulation::updateGrid() { //STEP 2 IN THE MPM GUIDE
    //EVENTUALLY THIS WILL HAVE INTERNAL FORCES FROM STRESS, BOUNDARY CONDITIONS, DAMPING 

    //apply gravity
    for (auto& node : grid) {
        if (node.mass > 0.0f) {
            node.force += glm::vec3(0.0f, -9.8f, 0.0f) * node.mass;
        }
    }
}


void MPMSimulation::transferGridToParticles() { //STEP 3 IN MPM GUIDE (basically same as P2G but in opposite direction?)
    for (auto& p : particles) {
        glm::vec3 cellIdx = (p.position / gridSpacing) - 0.5f;
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));
        
        p.velocity = glm::vec3(0.0f);
        p.force = glm::vec3(0.0f);
        
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    glm::ivec3 nodeIdx = baseNode + glm::ivec3(i-1, j-1, k-1);
                    if (nodeIdx.x < 0 || nodeIdx.x >= gridSize || nodeIdx.y < 0 || nodeIdx.y >= gridSize ||nodeIdx.z < 0 || nodeIdx.z >= gridSize) {
                        continue;
                    }
                    
                    int linearIdx = nodeIdx.x + nodeIdx.y * gridSize + nodeIdx.z * gridSize * gridSize;
                    glm::vec3 nodePos = glm::vec3(nodeIdx) * gridSpacing;
                    
                    float weight = computeWeight(p.position, nodePos);
                    p.velocity += grid[linearIdx].velocity * weight;
                    p.force += grid[linearIdx].force * weight;
                }
            }
        }
    }
}

float MPMSimulation::computeWeight(const glm::vec3& particlePos, const glm::vec3& nodePos) { 
    //RIGHT NOW WE ARE DOING LINEAR, WE NEED TO IMPLEMENT B-SPLINE HERE
    glm::vec3 d = (particlePos - nodePos) / gridSpacing;
    float wx = std::max(0.0f, 1.0f - std::abs(d.x));
    float wy = std::max(0.0f, 1.0f - std::abs(d.y));
    float wz = std::max(0.0f, 1.0f - std::abs(d.z));
    return wx * wy * wz;
}

glm::vec3 MPMSimulation::computeWeightGradient(const glm::vec3& particlePos, const glm::vec3& nodePos) {
    //NEED TO USE THIS FUNCTION LATER ON IN APIC/MLS AND B-SPLINE CALCULATIONS
    return glm::vec3(0.0f);
} 
