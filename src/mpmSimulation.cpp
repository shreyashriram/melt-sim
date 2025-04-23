#include "mpmSimulation.h"
#include <iostream>
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "shader_utils.h"

// Add a clamp implementation if not using C++17
template <typename T>
T clamp(const T& value, const T& lower, const T& upper) {
    return std::max(lower, std::min(value, upper));
}

MPMSimulation::MPMSimulation() : 
    dx(1.0f / gridSize),
    invDx(1.0f / dx),
    particleVAO(0),
    particleVBO(0) {
    
    // Compute Lame parameters
    mu = youngsModulus / (2.0f * (1.0f + poissonRatio));
    lambda = youngsModulus * poissonRatio / ((1.0f + poissonRatio) * (1.0f - 2.0f * poissonRatio));
    
    // Initialize grid
    grid.resize((gridSize + 1) * (gridSize + 1) * (gridSize + 1));
    
    // Create shader for particles
    const char* vertexShaderSource = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aColor;
        layout (location = 2) in float aSize;
        
        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        
        out vec3 particleColor;
        
        void main() {
            gl_Position = projection * view * vec4(aPos, 1.0);
            gl_PointSize = aSize;
            particleColor = aColor;
        }
    )";
    
    const char* fragmentShaderSource = R"(
        #version 330 core
        in vec3 particleColor;
        out vec4 FragColor;
        
        void main() {
            // Calculate distance from center of point
            vec2 coord = gl_PointCoord - vec2(0.5);
            if(length(coord) > 0.5)
                discard;
                
            FragColor = vec4(particleColor, 1.0);
        }
    )";
    
    // Create and compile shaders
    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    
    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    
    // Create shader program
    particleShaderProgram = glCreateProgram();
    glAttachShader(particleShaderProgram, vertexShader);
    glAttachShader(particleShaderProgram, fragmentShader);
    glLinkProgram(particleShaderProgram);
    
    // Delete shaders after linking
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    
    // Create vertex array for particles
    glGenVertexArrays(1, &particleVAO);
    glGenBuffers(1, &particleVBO);
}

MPMSimulation::~MPMSimulation() {
    glDeleteVertexArrays(1, &particleVAO);
    glDeleteBuffers(1, &particleVBO);
    glDeleteProgram(particleShaderProgram);
}

void MPMSimulation::initialize(const std::vector<Particle>& initialParticles) {
    // Convert our particles to material points
    particles.clear();
    particles.reserve(initialParticles.size());
    
    for (const auto& p : initialParticles) {
        MaterialPoint mp;
        mp.position = p.position;
        mp.velocity = p.velocity;
        mp.temperature = p.temperature / 100.0f; // Scale to 0-1 range
        mp.color = glm::vec3(0.7f, 0.2f, 0.2f); // Red for solids
        
        // Handle particle mass if needed
        // The MPM simulation can use the particle mass from the input particles
        // if particleMass needs to be variable instead of constant
        
        // Handle fixed particles if necessary
        if (p.isFixed) {
            // You may want to mark these particles differently or apply special constraints
        }
        
        particles.push_back(mp);
    }
    
    std::cout << "MPM simulation initialized with " << particles.size() << " particles" << std::endl;
}

void MPMSimulation::update(float deltaTime) {
    // Split the time step if needed for stability
    const float maxDt = 1e-4f;
    int steps = std::max(1, static_cast<int>(deltaTime / maxDt));
    float dt = deltaTime / steps;
    
    for (int step = 0; step < steps; step++) {
        resetGrid();
        particleToGrid(dt);  // Pass dt to particleToGrid
        updateGrid(dt);
        gridToParticle(dt);
    }
}

void MPMSimulation::resetGrid() {
    // Reset grid to initial state
    for (auto& node : grid) {
        node.velocity = glm::vec3(0.0f);
        node.mass = 0.0f;
    }
}

void MPMSimulation::particleToGrid(float dt) {  // Add dt parameter
    // Transfer mass and momentum from particles to grid
    for (auto& p : particles) {
        // Convert position to grid space
        glm::vec3 gridPos = p.position * invDx;
        
        // Base grid node (bottom-left-back corner)
        glm::ivec3 baseNode = glm::ivec3(gridPos - 0.5f);
        
        // Fractional offset from base node
        glm::vec3 fracPos = gridPos - glm::vec3(baseNode) - 0.5f;
        
        // Quadratic kernels for weights
        glm::vec3 weights[3];
        weights[0] = 0.5f * glm::pow(1.5f - fracPos, glm::vec3(2.0f));
        weights[1] = 0.75f - glm::pow(fracPos - glm::vec3(1.0f), glm::vec3(2.0f));
        weights[2] = 0.5f * glm::pow(fracPos - glm::vec3(0.5f), glm::vec3(2.0f));
        
        // Update material state (solid->fluid transition)
        updateMaterialState(p, 1e-4f);
        
        // Compute stress based on material model
        float J = glm::determinant(p.deformationGradient);
        
        // Decompose deformation gradient into rotation and strain
        glm::mat3 R, S;
        // TODO: Implement polar decomposition
        // For now we'll just use identity for R
        R = glm::mat3(1.0f);
        
        // Compute stress tensor
        glm::mat3 stress;
        
        if (p.fluidity < 0.01f) {
            // Solid - use corotated elasticity model
            stress = -dt * particleVolume * (
                2.0f * mu * (p.deformationGradient - R) * glm::transpose(p.deformationGradient) + 
                lambda * (J - 1.0f) * J * glm::mat3(1.0f)
            );
        } else if (p.fluidity > 0.99f) {
            // Fluid - use weak compressibility
            float pressure = lambda * (J - 1.0f);
            stress = -dt * particleVolume * (pressure * J * glm::mat3(1.0f));
        } else {
            // Blend between solid and fluid
            glm::mat3 solidStress = 2.0f * mu * (p.deformationGradient - R) * glm::transpose(p.deformationGradient) + 
                                   lambda * (J - 1.0f) * J * glm::mat3(1.0f);
            float pressure = lambda * (J - 1.0f);
            glm::mat3 fluidStress = pressure * J * glm::mat3(1.0f);
            stress = -dt * particleVolume * ((1.0f - p.fluidity) * solidStress + p.fluidity * fluidStress);
        }
        
        // Affine momentum matrix
        glm::mat3 affine = stress + particleMass * p.affineVelocity;
        
        // Transfer to grid nodes
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    int nodeIdx = (baseNode.x + i) * (gridSize + 1) * (gridSize + 1) + 
                                 (baseNode.y + j) * (gridSize + 1) + 
                                 (baseNode.z + k);
                    
                    // Skip out-of-bounds nodes
                    if (baseNode.x + i < 0 || baseNode.x + i > gridSize ||
                        baseNode.y + j < 0 || baseNode.y + j > gridSize ||
                        baseNode.z + k < 0 || baseNode.z + k > gridSize) {
                        continue;
                    }
                    
                    // Compute weight
                    float weight = weights[i].x * weights[j].y * weights[k].z;
                    
                    // Distance from node to particle
                    glm::vec3 dist = (glm::vec3(i, j, k) - fracPos) * dx;
                    
                    // Momentum and mass transfer
                    grid[nodeIdx].velocity += weight * (p.velocity * particleMass + glm::vec3(affine * dist));
                    grid[nodeIdx].mass += weight * particleMass;
                }
            }
        }
    }
}

void MPMSimulation::updateGrid(float dt) {
    // Update grid velocities and apply forces
    for (int i = 0; i <= gridSize; i++) {
        for (int j = 0; j <= gridSize; j++) {
            for (int k = 0; k <= gridSize; k++) {
                int idx = i * (gridSize + 1) * (gridSize + 1) + j * (gridSize + 1) + k;
                
                auto& node = grid[idx];
                
                if (node.mass > 0) {
                    // Convert momentum to velocity
                    node.velocity /= node.mass;
                    
                    // Apply gravity
                    node.velocity += dt * glm::vec3(0.0f, -9.8f, 0.0f);
                    
                    // Handle boundary conditions
                    float boundary = 0.05f;
                    float x = (float)i / gridSize;
                    float y = (float)j / gridSize;
                    float z = (float)k / gridSize;
                    
                    // Sticky boundaries (walls)
                    if (x < boundary || x > 1.0f - boundary || 
                        y > 1.0f - boundary || 
                        z < boundary || z > 1.0f - boundary) {
                        node.velocity = glm::vec3(0.0f);
                    }
                    
                    // Separate boundary (floor with friction)
                    if (y < boundary) {
                        node.velocity.y = std::max(0.0f, node.velocity.y);
                        // Apply friction
                        node.velocity.x *= 0.95f;
                        node.velocity.z *= 0.95f;
                    }
                }
            }
        }
    }
}

void MPMSimulation::gridToParticle(float dt) {
    // Transfer velocities back to particles and update
    for (auto& p : particles) {
        // Reset particle velocity and affine velocity
        p.velocity = glm::vec3(0.0f);
        p.affineVelocity = glm::mat3(0.0f);
        
        // Convert position to grid space
        glm::vec3 gridPos = p.position * invDx;
        
        // Base grid node
        glm::ivec3 baseNode = glm::ivec3(gridPos - 0.5f);
        
        // Fractional offset from base node
        glm::vec3 fracPos = gridPos - glm::vec3(baseNode) - 0.5f;
        
        // Quadratic kernels for weights
        glm::vec3 weights[3];
        weights[0] = 0.5f * glm::pow(1.5f - fracPos, glm::vec3(2.0f));
        weights[1] = 0.75f - glm::pow(fracPos - glm::vec3(1.0f), glm::vec3(2.0f));
        weights[2] = 0.5f * glm::pow(fracPos - glm::vec3(0.5f), glm::vec3(2.0f));
        
        // Gather from grid nodes
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    int nodeIdx = (baseNode.x + i) * (gridSize + 1) * (gridSize + 1) + 
                                 (baseNode.y + j) * (gridSize + 1) + 
                                 (baseNode.z + k);
                    
                    // Skip out-of-bounds nodes
                    if (baseNode.x + i < 0 || baseNode.x + i > gridSize ||
                        baseNode.y + j < 0 || baseNode.y + j > gridSize ||
                        baseNode.z + k < 0 || baseNode.z + k > gridSize) {
                        continue;
                    }
                    
                    // Compute weight
                    float weight = weights[i].x * weights[j].y * weights[k].z;
                    
                    // Grid node to local space
                    glm::vec3 dist = glm::vec3(i, j, k) - fracPos;
                    
                    // Gather velocity
                    glm::vec3 nodeVel = grid[nodeIdx].velocity;
                    p.velocity += weight * nodeVel;
                    
                    // Affine velocity update (APIC)
                    p.affineVelocity += 4.0f * invDx * weight * glm::outerProduct(nodeVel, dist);
                }
            }
        }
        
        // Update position
        p.position += dt * p.velocity;
        
        // Update deformation gradient
        glm::mat3 F = p.deformationGradient;
        glm::mat3 updatedF = (glm::mat3(1.0f) + dt * p.affineVelocity) * F;
        
        // Different treatment for solid and fluid
        if (p.fluidity < 0.5f) {
            // For solid, apply plasticity
            // TODO: Implement SVD and clamping of singular values
            
            // Just do simple clamping for now
            float J = glm::determinant(updatedF);
            float clampedJ = clamp(J, 0.6f, 20.0f);  // Use our template function instead of std::clamp
            updatedF *= std::pow(clampedJ / J, 1.0f/3.0f);
            
            p.plasticVolume = clampedJ;
        } else {
            // For fluid, only keep volume change
            float J = glm::determinant(updatedF);
            updatedF = glm::mat3(std::pow(J, 1.0f/3.0f));
            p.plasticVolume = 1.0f;
        }
        
        p.deformationGradient = updatedF;
        
        // Update color based on material state
        if (p.fluidity < 0.1f) {
            // Solid
            p.color = glm::vec3(0.7f, 0.2f, 0.2f);
        } else if (p.fluidity < 0.5f) {
            // Transitioning
            float t = p.fluidity / 0.5f;
            p.color = (1.0f - t) * glm::vec3(0.7f, 0.2f, 0.2f) + t * glm::vec3(0.2f, 0.4f, 0.7f);
        } else {
            // Fluid
            p.color = glm::vec3(0.2f, 0.4f, 0.7f);
        }
    }
}

void MPMSimulation::updateMaterialState(MaterialPoint& p, float dt) {
    // Compute stress for failure criterion
    float J = glm::determinant(p.deformationGradient);
    
    // Simple von Mises stress approximation
    float vonMisesStress = 0.0f;
    
    // Increase temperature based on deformation
    p.temperature += vonMisesStress * 0.001f;
    
    // Apply heat at bottom
    if (p.position.y < 0.3f) {
        p.temperature += dt * 0.1f;
    }
    
    // Determine if material should start melting
    if (vonMisesStress > yieldStress || p.temperature > 0.4f) {
        p.fluidity += decayRate * dt * 10.0f;
        p.fluidity = clamp(p.fluidity, 0.0f, 1.0f);
    }
    if (p.position.y < 0.3f) {
        p.temperature += dt * 0.5f;  // 5x stronger heating (from 0.1f)
    }
}

void MPMSimulation::render(const glm::mat4& view, const glm::mat4& projection) {
    // Prepare particle data for rendering
    std::vector<float> particleData;
    particleData.reserve(particles.size() * 7); // pos(3) + color(3) + size(1)
    
    for (const auto& p : particles) {
        // Position
        particleData.push_back(p.position.x);
        particleData.push_back(p.position.y);
        particleData.push_back(p.position.z);
        
        // Color
        particleData.push_back(p.color.r);
        particleData.push_back(p.color.g);
        particleData.push_back(p.color.b);
        
        // Size (larger for fluid, smaller for solid)
        float size = 5.0f + p.fluidity * 5.0f;
        particleData.push_back(size);
    }
    
    // Update VBO
    glBindVertexArray(particleVAO);
    glBindBuffer(GL_ARRAY_BUFFER, particleVBO);
    glBufferData(GL_ARRAY_BUFFER, particleData.size() * sizeof(float), particleData.data(), GL_DYNAMIC_DRAW);
    
    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    // Color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    
    // Size attribute
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
    
    // Render particles
    glUseProgram(particleShaderProgram);
    glUniformMatrix4fv(glGetUniformLocation(particleShaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(particleShaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(glGetUniformLocation(particleShaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(glm::mat4(1.0f)));
    
    // Enable point sprites and blending
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glDrawArrays(GL_POINTS, 0, particles.size());
    
    // Cleanup
    glDisable(GL_BLEND);
    glBindVertexArray(0);
}