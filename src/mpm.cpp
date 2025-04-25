#include "mpm.h"
#include "particle.h"
#include "grid.h"
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>


bool alreadyPrinted = false; //for debugging 

MPMSimulation::MPMSimulation() 
    : youngsModulus(1.4e5f), poissonsRatio(0.2f), gridSize(64), gridSpacing(0.25f) {
}

void MPMSimulation::addMeshParticles(std::vector<Vector3> sampledPoints) {
    Particle p;
    for (auto& pt : sampledPoints) {
        p = Particle(glm::vec3(pt.x(), pt.y() + 1, pt.z()), glm::vec3(0.0f, 0.0f, 0.0f));
        particles.push_back(p);
    }
}

void MPMSimulation::initializeParticles() {
    particles.clear(); 
    
    Particle p;  //declare the Particle variable
    for (int i = 0; i < 5; i++) {
        p = Particle(glm::vec3(i*0.2f, 0.5f, 0.0f), glm::vec3(0.0f, i*0.05f, 0.0f));

        particles.push_back(p);
    }
}


void MPMSimulation::step(float dt) {
    // reset grid
    // for (auto& node : grid) {
    //     node.velocity = glm::vec3(0.0f);
    //     node.force = glm::vec3(0.0f);
    //     node.mass = 0.0f;
    // }
    
    //transferParticlesToGrid();
    //updateGrid();
    //transferGridToParticles();


    updateParticles(dt, gridSize, gridSpacing);
}

void MPMSimulation::transferParticlesToGrid() { //STEP 1 IN MPM GUIDE
    for (const auto& p : particles) {
        // Offset the position to handle negative coordinates
        glm::vec3 offsetPos = p.position + glm::vec3(gridSize * gridSpacing * 0.5f);
        glm::vec3 cellIdx = (offsetPos / gridSpacing) - 0.5f;
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));
        
        for (int i = 0; i < 3; ++i) { //bc we can't do 1.5 in an int loop, still loop through 3 cells and just adjust the index to be the correct neighborhood
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    glm::ivec3 nodeIdx = baseNode + glm::ivec3(i-1, j-1, k-1); //current node in the neighborhood
                    
                    if (nodeIdx.x < 0 || nodeIdx.x >= gridSize || nodeIdx.y < 0 || nodeIdx.y >= gridSize || nodeIdx.z < 0 || nodeIdx.z >= gridSize) {
                        continue; //if not in the grid bounds, skip it 
                    }
                    
                    int linearIdx = nodeIdx.x + nodeIdx.y * gridSize + nodeIdx.z * gridSize * gridSize; //finding the index in the array
                    glm::vec3 nodePos = (glm::vec3(nodeIdx) * gridSpacing) - glm::vec3(gridSize * gridSpacing * 0.5f);; //converting from grid space to world space 
                    
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
        glm::vec3 offsetPos = p.position + glm::vec3(gridSize * gridSpacing * 0.5f);
        glm::vec3 cellIdx = (offsetPos / gridSpacing) - 0.5f;
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
                    glm::vec3 nodePos = (glm::vec3(nodeIdx) * gridSpacing) - glm::vec3(gridSize * gridSpacing * 0.5f);
                    
                    float weight = computeWeight(p.position, nodePos);
                    p.velocity += grid[linearIdx].velocity * weight;
                    p.force += grid[linearIdx].force * weight;
                }
            }
        }
    }
}

void MPMSimulation::updateParticles(float dt, int gridSize, float gridSpacing) {
    float gridBoundary = gridSize * gridSpacing;
    
    for (auto& p : particles) {
        // Update velocity based on forces
        p.velocity += glm::vec3(0.0f, -9.81f, 0.0f) * dt;
        
        // Update position
        p.position += p.velocity * dt;
        if(p.position.y < 0.0f) {
            p.position.y = 0.0f;
            p.velocity.y *= -0.5f;  // Bounce with energy loss
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
