#include "mpm.h"
#include "particle.h"
#include "grid.h"
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

bool alreadyPrinted = false; // for debugging


MPMSimulation::MPMSimulation() 
    : youngsModulus(1.4e5f), poissonsRatio(0.2f), grid(5, 0.25f){
}


// TODO: make sure its in grid, currently manually translated
void MPMSimulation::addMeshParticles(std::vector<Vector3> sampledPoints) {
    Particle p;
    for (auto& pt : sampledPoints) {

        // ? Debugging: negative y velocity
        p = Particle(glm::vec3(pt.x()+0.75f, pt.y()+0.75f, pt.z()+0.75f), glm::vec3(0.0f, 0.0f, 0.0f));
        particles.push_back(p);#include "mpm.h"
#include "particle.h"
#include "grid.h"
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

bool alreadyPrinted = false; // for debugging


MPMSimulation::MPMSimulation() 
    : youngsModulus(1.4e5f), poissonsRatio(0.2f), grid(5, 0.25f){
}


// TODO: make sure its in grid, currently manually translated
void MPMSimulation::addMeshParticles(std::vector<Vector3> sampledPoints) {
    Particle p;
    for (auto& pt : sampledPoints) {

        // ? Debugging: negative y velocity
        p = Particle(glm::vec3(pt.x()+0.75f, pt.y()+0.75f, pt.z()+0.75f), glm::vec3(0.0f, 0.0f, 0.0f));
        particles.push_back(p);
    }
}

// ? Debugging: Manually create particles
void MPMSimulation::initializeParticles() {
    particles.clear(); 
    
    Particle p;  //declare the Particle variable
    for (int i = 0; i < 5; i++) {
        p = Particle(glm::vec3(i*0.2f, 0.5f, 0.0f), glm::vec3(0.0f, i*0.05f, 0.0f));

        particles.push_back(p);
    }
}

void MPMSimulation::step(float dt) {

    // std::cout << "Step start -- dt: " << dt << std::endl;
    // std::cout << "First particle pos: " << particles[0].position.y << ", vel: " << particles[0].velocity.y << std::endl;

    // reset grid
    for (auto& n : grid.nodes) {
        n.velocity = glm::vec3(0.0f);
        n.force = glm::vec3(0.0f);
        n.mass = 0.0f;
    }
    
    transferParticlesToGrid();
    updateGrid(dt);
    transferGridToParticles(dt);
    updateParticles(dt);
}

// STEP 1 IN MPM GUIDE
void MPMSimulation::transferParticlesToGrid() { 

    for (const auto& p : particles) {
        // std::cout << "Particle Pos: (" << p.position.x << ", " << p.position.y << ", " << p.position.z << ")" << std::endl;

        
        glm::vec3 cellIdx = (p.position / grid.spacing);
        
        // bottom left node relative to particle
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));
        
        // Loop through the 3x3x3 neighborhood using -1 to 1 range directly
        // 1 node apart.
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {

                    glm::ivec3 nodeIdx = baseNode + glm::ivec3(i, j, k);
                    
                    if (nodeIdx.x < 0 || nodeIdx.x >= grid.size || 
                        nodeIdx.y < 0 || nodeIdx.y >= grid.size || 
                        nodeIdx.z < 0 || nodeIdx.z >= grid.size) {


                        continue; // if not in the grid bounds, skip it 
                    }


                    int linearIdx = nodeIdx.x + (nodeIdx.y * grid.size) + (nodeIdx.z * grid.size * grid.size); // finding the index in the <node>
                
                    // node(0, 0, 0) is at (0, 0, 0)w
                    glm::vec3 nodePos = (glm::vec3(nodeIdx) * grid.spacing);
                    
                    float weight = computeWeight(p.position, nodePos); // where we will implement B-spline

    
                    grid.nodes[linearIdx].mass += p.mass * weight; // setting the mass and velocity in the array 
                    grid.nodes[linearIdx].velocity += p.mass * p.velocity * weight;
                    
                    // WE STILL NEED TO DO EITHER APIC OR MLS FOR PRESERVED MOMENTUM (INSTEAD OF LINEAR VELOCITY LIKE WE ARE RIGHT NOW)
                }
            }
        }
    }
    
    //have to normalize velocity to convert accumulated momentum back to velocity
    for (auto& node : grid.nodes) {
        if (node.mass > 0.0f) {
            node.velocity /= node.mass;
        }
    }

    float totalMass = 0.0f;
    for (const auto& n : grid.nodes) {
        totalMass += n.mass;
    }
    std::cout << "Total grid mass after P2G: " << totalMass << std::endl;

    // for (int i = 0; i < grid.nodes.size(); ++i) {
    //     if (grid.nodes[i].mass > 0.0f) {
    //         std::cout << "Node " << i << " vel.y: " << grid.nodes[i].velocity.y << std::endl;
    //     }
    // }
}

// STEP 2 IN THE MPM GUIDE
void MPMSimulation::updateGrid(float dt ) { 
    //EVENTUALLY THIS WILL HAVE INTERNAL FORCES FROM STRESS, BOUNDARY CONDITIONS, DAMPING 
    for (auto& node : grid.nodes) {
        if (node.mass > 0.0f) {
            // Apply gravity
            node.force += glm::vec3(0.0f, -9.8f, 0.0f) * node.mass;
            node.velocity += (node.force / node.mass) * dt;  // <-- Add this line

        }
    
  }
}

// STEP 3 IN MPM GUIDE (basically same as P2G but in opposite direction?)
void MPMSimulation::transferGridToParticles(float dt) { 
    for (auto& p : particles) {
        
        glm::vec3 cellIdx = (p.position / grid.spacing);
        
        // bottom left node relative to particle
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));


        glm::vec3 newVelocity(0.0f);
        
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    glm::ivec3 nodeIdx = baseNode + glm::ivec3(i, j, k);
                    if (nodeIdx.x < 0 || nodeIdx.x >= grid.size || nodeIdx.y < 0 || nodeIdx.y >= grid.size ||nodeIdx.z < 0 || nodeIdx.z >= grid.size) {
                        continue;
                    }

                    int linearIdx = nodeIdx.x + nodeIdx.y * grid.size + nodeIdx.z * grid.size * grid.size;
                    glm::vec3 nodePos = glm::vec3(nodeIdx) * grid.spacing;
                    
                    float weight = computeWeight(p.position, nodePos);
                    newVelocity += grid.nodes[linearIdx].velocity * weight;
                }
            }
        }
        
        p.velocity = newVelocity;
    }
}


void MPMSimulation::updateParticles(float dt) {
    float gridBoundary = grid.size * grid.spacing;
    
    for (auto& p : particles) {
        // Update velocity based on forces
        p.position += p.velocity * dt;
        
        // Floor collision
        if (p.position.y < 0.0f) {
            p.position.y = 0.0f;
            p.velocity.y *= -0.5f; // Bounce with energy loss
        }
    }
}

float MPMSimulation::computeWeight(const glm::vec3& particlePos, const glm::vec3& nodePos) { 
    //RIGHT NOW WE ARE DOING LINEAR, WE NEED TO IMPLEMENT B-SPLINE HERE
    glm::vec3 diff = (particlePos - nodePos) / grid.spacing; // Normalize distance by grid spacing
    glm::vec3 w;
    
    for (int i = 0; i < 3; ++i) {
        
        float d = fabs(diff[i]);
        
        if (d < 0.5f) {
            w[i] = 0.75f - d * d;
        }
        else if (d < 1.5f) {
            float t = 1.5f - d;
            w[i] = 0.5f * t * t;
        } 
        else {
            w[i] = 0.0f;
        }
    }

    return w.x * w.y * w.z;  // Tensor product of 1D kernels
}

// glm::vec3 MPMSimulation::computeWeightGradient(const glm::vec3& particlePos, const glm::vec3& nodePos) {
//     //NEED TO USE THIS FUNCTION LATER ON IN APIC/MLS AND B-SPLINE CALCULATIONS
//     return glm::vec3(0.0f);
// } 

    }
}

// ? Debugging: Manually create particles
void MPMSimulation::initializeParticles() {
    particles.clear(); 
    
    Particle p;  //declare the Particle variable
    for (int i = 0; i < 5; i++) {
        p = Particle(glm::vec3(i*0.2f, 0.5f, 0.0f), glm::vec3(0.0f, i*0.05f, 0.0f));

        particles.push_back(p);
    }
}

void MPMSimulation::step(float dt) {

    // std::cout << "Step start -- dt: " << dt << std::endl;
    // std::cout << "First particle pos: " << particles[0].position.y << ", vel: " << particles[0].velocity.y << std::endl;

    // reset grid
    for (auto& n : grid.nodes) {
        n.velocity = glm::vec3(0.0f);
        n.force = glm::vec3(0.0f);
        n.mass = 0.0f;
    }
    
    transferParticlesToGrid();
    updateGrid(dt);
    transferGridToParticles(dt);
    updateParticles(dt);
}

// STEP 1 IN MPM GUIDE
void MPMSimulation::transferParticlesToGrid() { 

    for (const auto& p : particles) {
        // std::cout << "Particle Pos: (" << p.position.x << ", " << p.position.y << ", " << p.position.z << ")" << std::endl;

        
        glm::vec3 cellIdx = (p.position / grid.spacing);
        
        // bottom left node relative to particle
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));
        
        // Loop through the 3x3x3 neighborhood using -1 to 1 range directly
        // 1 node apart.
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {

                    glm::ivec3 nodeIdx = baseNode + glm::ivec3(i, j, k);
                    
                    if (nodeIdx.x < 0 || nodeIdx.x >= grid.size || 
                        nodeIdx.y < 0 || nodeIdx.y >= grid.size || 
                        nodeIdx.z < 0 || nodeIdx.z >= grid.size) {


                        continue; // if not in the grid bounds, skip it 
                    }


                    int linearIdx = nodeIdx.x + (nodeIdx.y * grid.size) + (nodeIdx.z * grid.size * grid.size); // finding the index in the <node>
                
                    // node(0, 0, 0) is at (0, 0, 0)w
                    glm::vec3 nodePos = (glm::vec3(nodeIdx) * grid.spacing);
                    
                    float weight = computeWeight(p.position, nodePos); // where we will implement B-spline

    
                    grid.nodes[linearIdx].mass += p.mass * weight; // setting the mass and velocity in the array 
                    grid.nodes[linearIdx].velocity += p.mass * p.velocity * weight;
                    
                    // WE STILL NEED TO DO EITHER APIC OR MLS FOR PRESERVED MOMENTUM (INSTEAD OF LINEAR VELOCITY LIKE WE ARE RIGHT NOW)
                }
            }
        }
    }
    
    //have to normalize velocity to convert accumulated momentum back to velocity
    for (auto& node : grid.nodes) {
        if (node.mass > 0.0f) {
            node.velocity /= node.mass;
        }
    }

    float totalMass = 0.0f;
    for (const auto& n : grid.nodes) {
        totalMass += n.mass;
    }
    std::cout << "Total grid mass after P2G: " << totalMass << std::endl;

    // for (int i = 0; i < grid.nodes.size(); ++i) {
    //     if (grid.nodes[i].mass > 0.0f) {
    //         std::cout << "Node " << i << " vel.y: " << grid.nodes[i].velocity.y << std::endl;
    //     }
    // }
}

// STEP 2 IN THE MPM GUIDE
void MPMSimulation::updateGrid(float dt ) { 
    //EVENTUALLY THIS WILL HAVE INTERNAL FORCES FROM STRESS, BOUNDARY CONDITIONS, DAMPING 
    for (auto& node : grid.nodes) {
        if (node.mass > 0.0f) {
            // Apply gravity
            node.force += glm::vec3(0.0f, -9.8f, 0.0f) * node.mass;
            node.velocity += (node.force / node.mass) * dt;  // <-- Add this line

        }
    
  }
}

// STEP 3 IN MPM GUIDE (basically same as P2G but in opposite direction?)
void MPMSimulation::transferGridToParticles(float dt) { 
    for (auto& p : particles) {
        
        glm::vec3 cellIdx = (p.position / grid.spacing);
        
        // bottom left node relative to particle
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));


        glm::vec3 newVelocity(0.0f);
        
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    glm::ivec3 nodeIdx = baseNode + glm::ivec3(i, j, k);
                    if (nodeIdx.x < 0 || nodeIdx.x >= grid.size || nodeIdx.y < 0 || nodeIdx.y >= grid.size ||nodeIdx.z < 0 || nodeIdx.z >= grid.size) {
                        continue;
                    }

                    int linearIdx = nodeIdx.x + nodeIdx.y * grid.size + nodeIdx.z * grid.size * grid.size;
                    glm::vec3 nodePos = glm::vec3(nodeIdx) * grid.spacing;
                    
                    float weight = computeWeight(p.position, nodePos);
                    newVelocity += grid.nodes[linearIdx].velocity * weight;
                }
            }
        }
        
        p.velocity = newVelocity;
    }
}


void MPMSimulation::updateParticles(float dt) {
    float gridBoundary = grid.size * grid.spacing;
    
    for (auto& p : particles) {
        // Update velocity based on forces
        p.position += p.velocity * dt;
        
        // Floor collision
        if (p.position.y < 0.0f) {
            p.position.y = 0.0f;
            p.velocity.y *= -0.5f; // Bounce with energy loss
        }
    }
}

float MPMSimulation::computeWeight(const glm::vec3& particlePos, const glm::vec3& nodePos) { 
    //RIGHT NOW WE ARE DOING LINEAR, WE NEED TO IMPLEMENT B-SPLINE HERE
    glm::vec3 diff = (particlePos - nodePos) / grid.spacing; // Normalize distance by grid spacing
    glm::vec3 w;
    
    for (int i = 0; i < 3; ++i) {
        
        float d = fabs(diff[i]);
        
        if (d < 0.5f) {
            w[i] = 0.75f - d * d;
        }
        else if (d < 1.5f) {
            float t = 1.5f - d;
            w[i] = 0.5f * t * t;
        } 
        else {
            w[i] = 0.0f;
        }
    }

    return w.x * w.y * w.z;  // Tensor product of 1D kernels
}

// glm::vec3 MPMSimulation::computeWeightGradient(const glm::vec3& particlePos, const glm::vec3& nodePos) {
//     //NEED TO USE THIS FUNCTION LATER ON IN APIC/MLS AND B-SPLINE CALCULATIONS
//     return glm::vec3(0.0f);
// } 
