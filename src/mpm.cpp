#include "mpm.h"
#include "particle.h"
#include "grid.h"
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>



bool alreadyPrinted = false; //for debugging 

MPMSimulation::MPMSimulation() 
    : youngsModulus(1.4e5f), poissonsRatio(0.2f), grid(5, 0.25f),yieldThreshold(0.01f), meltRate(0.05f), globalMeltProgress(0.0f) {
    
    //calculate Lamé parameters from Young's modulus and Poisson's ratio
    shearModulus = youngsModulus / (2.0f * (1.0f + poissonsRatio));
    bulkModulus = youngsModulus * poissonsRatio / ((1.0f + poissonsRatio) * (1.0f - 2.0f * poissonsRatio));
    
    if (poissonsRatio > 0.495f) { //if poisson's ratio is too close to 0.5, set bulk modulus to a large value
        bulkModulus = 1000.0f * shearModulus; //approximate for numerical stability
    }
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
    //std::cout << "Total grid mass after P2G: " << totalMass << std::endl;

    // for (int i = 0; i < grid.nodes.size(); ++i) {
    //     if (grid.nodes[i].mass > 0.0f) {
    //         std::cout << "Node " << i << " vel.y: " << grid.nodes[i].velocity.y << std::endl;
    //     }
    // }
}

// STEP 2 IN THE MPM GUIDE
void MPMSimulation::updateGrid(float dt ) { 
    //TODO: computeStress and update deformation gradient in this function 
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
    //TODO: call apply plasticity in this function 
    float gridBoundary = grid.size * grid.spacing;
    
    for (auto& p : particles) {
        // Update velocity based on forces
        p.position += p.velocity * dt;
        
        // Floor collision
        if (p.position.y < 0.0f) {
            p.position.y = 0.0f;
            p.velocity.y *= -0.5f;  // Bounce with energy loss
        }
    }
}

float MPMSimulation::computeWeight(const glm::vec3& particlePos, const glm::vec3& nodePos) { 
    glm::vec3 diff = (particlePos - nodePos) / grid.spacing;
    glm::vec3 w;
    
    for (int i = 0; i < 3; ++i) { //B-spline weight function
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

    return w.x * w.y * w.z;
}


glm::vec3 MPMSimulation::computeWeightGradient(const glm::vec3& particlePos, const glm::vec3& nodePos) {
    glm::vec3 rel = (particlePos - nodePos) / grid.spacing;
    glm::vec3 grad(0.0f);

    std::cout << "=== Weight Gradient Debug ===" << std::endl;
    std::cout << "Particle Position: (" << particlePos.x << ", " << particlePos.y << ", " << particlePos.z << ")" << std::endl;
    std::cout << "Node Position: (" << nodePos.x << ", " << nodePos.y << ", " << nodePos.z << ")" << std::endl;
    std::cout << "Difference (rel): (" << rel.x << ", " << rel.y << ", " << rel.z << ")" << std::endl;

    for (int dim = 0; dim < 3; ++dim) {
        float x = rel[dim];
        float absX = fabs(x);
        float sign = (x >= 0.0f) ? 1.0f : -1.0f;

        // Compute gradient along this dimension
        float gradW = 0.0f;
        if (absX < 0.5f) {
            gradW = -2.0f * x;
        } else if (absX < 1.5f) {
            gradW = -sign * (1.5f - absX);
        }

        int a = (dim + 1) % 3;
        int b = (dim + 2) % 3;

        auto weight = [](float x) -> float {
            float absX = fabs(x);
            if (absX < 0.5f) {
                return 0.75f - absX * absX;
            } else if (absX < 1.5f) {
                float d = 1.5f - absX;
                return 0.5f * d * d;
            } else {
                return 0.0f;
            }
        };

        float wa = weight(rel[a]);
        float wb = weight(rel[b]);

        grad[dim] = gradW * wa * wb;

        std::cout << "dim " << dim << ": x=" << x << ", gradW=" << gradW << ", weight[" << a << "]=" << wa << ", weight[" << b << "]=" << wb << std::endl;
        std::cout << "Partial Gradient[" << dim << "] = " << grad[dim] << std::endl;
    }

    grad /= grid.spacing;
    std::cout << "Final Gradient: (" << grad.x << ", " << grad.y << ", " << grad.z << ")" << std::endl;

    return grad;
}




void MPMSimulation::polarDecomposition(const glm::mat3& F, glm::mat3& R, glm::mat3& S){
    static bool firstDecomp = true;
    if (firstDecomp) {
        std::cout << "\n=== Polar Decomposition Debug ===" << std::endl;
        std::cout << "Input Matrix F:" << std::endl;
        for(int i = 0; i < 3; i++) {
            std::cout << "[ ";
            for(int j = 0; j < 3; j++) {
                std::cout << F[i][j] << " ";
            }
            std::cout << "]" << std::endl;
        }
    }
    
    Eigen::Matrix3f F_eigen; //converting the glm::mat3 matrix to Eigen matrix
    for(int i = 0; i < 3; i++) { 
        for(int j = 0; j < 3; j++){
            F_eigen(i, j) = F[i][j];
        }
    }

    //compute the SVD
    Eigen::JacobiSVD<Eigen::Matrix3f> svd(F_eigen, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3f U = svd.matrixU();
    Eigen::Matrix3f V = svd.matrixV();

    Eigen::Matrix3f R_eigen = U * V.transpose(); //get the rotation matrix

    //compute the S matrix (S = V * s * V^T)
    Eigen::Vector3f singularValues = svd.singularValues();
    Eigen::Matrix3f S_eigen = V * singularValues.asDiagonal() * V.transpose(); //get the diagonal matrix of singular values

    //convert back to glm::mat3
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++){
            R[i][j] = R_eigen(i, j);
            S[i][j] = S_eigen(i, j);
        }
    }

    if (firstDecomp) {
        std::cout << "\nSingular Values: " << singularValues.transpose() << std::endl;
        
        std::cout << "\nRotation Matrix R:" << std::endl;
        for(int i = 0; i < 3; i++) {
            std::cout << "[ ";
            for(int j = 0; j < 3; j++) {
                std::cout << R[i][j] << " ";
            }
            std::cout << "]" << std::endl;
        }
        
        std::cout << "\nStretch Matrix S:" << std::endl;
        for(int i = 0; i < 3; i++) {
            std::cout << "[ ";
            for(int j = 0; j < 3; j++) {
                std::cout << S[i][j] << " ";
            }
            std::cout << "]" << std::endl;
        }
        
        // Verify properties
        glm::mat3 reconstructed = R * S;
        float maxDiff = 0.0f;
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                maxDiff = std::max(maxDiff, std::abs(reconstructed[i][j] - F[i][j]));
            }
        }
        std::cout << "\nMax difference between F and R*S: " << maxDiff << std::endl;
        
        firstDecomp = false;
    }

}

glm::mat3 MPMSimulation::computeVelocityGradient(const Particle& p){
    return glm::mat3(0.0f);
}

void MPMSimulation::updatePlasticity(Particle& p){}


glm::mat3 MPMSimulation::computeStress(const Particle& p){
    return glm::mat3(0.0f);
}

void MPMSimulation::runTests() {
    std::cout << "\n=== Starting MPM Tests ===" << std::endl;
    
    // Test 1: Weight Gradient
    std::cout << "\nTest 1: Weight Gradient" << std::endl;
    grid.spacing = 0.25f;
    glm::vec3 particlePos(0.003f, 0.003f, 0.003f); // this is within the kernel support range
    glm::vec3 nodePos(0.0f, 0.0f, 0.0f); // node is at the origin    
    glm::vec3 gradient = computeWeightGradient(particlePos, nodePos);
    
    std::cout << "Expected: Gradient should point from node to particle" << std::endl;
    std::cout << "Gradient: (" << gradient.x << ", " << gradient.y << ", " << gradient.z << ")" << std::endl;
    
    // Test 2: Polar Decomposition
    std::cout << "\nTest 2: Polar Decomposition" << std::endl;
    // Create a known rotation matrix (45 degrees around Z-axis)
    float angle = glm::radians(45.0f);
    float c = cos(angle);
    float s = sin(angle);
    glm::mat3 rotation(
        c, -s, 0.0f,
        s,  c, 0.0f,
        0.0f, 0.0f, 1.0f
    );
    
    // Create a known stretch matrix
    glm::mat3 stretch(
        2.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 1.0f
    );
    
    // Compose deformation gradient F = R * S
    glm::mat3 F = rotation * stretch;
    
    // Decompose
    glm::mat3 R, S;
    polarDecomposition(F, R, S);
    
    // Verify results
    std::cout << "\nVerification:" << std::endl;
    std::cout << "Input matrix should be rotation * stretch" << std::endl;
    std::cout << "R should be close to rotation matrix" << std::endl;
    std::cout << "S should be close to stretch matrix" << std::endl;
    
    // Compute errors
    float rotationError = 0.0f;
    float stretchError = 0.0f;
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            rotationError = std::max(rotationError, std::abs(R[i][j] - rotation[i][j]));
            stretchError = std::max(stretchError, std::abs(S[i][j] - stretch[i][j]));
        }
    }
    
    std::cout << "Rotation error: " << rotationError << std::endl;
    std::cout << "Stretch error: " << stretchError << std::endl;
    
    std::cout << "\n=== Tests Complete ===" << std::endl;
}
/*FUNCTIONS TO ADD:
COMPUTE STRESS
UPDATE DEFORMATION GRADIENT 
APPLY PLASTICITY
COMPUTE INTERNAL FORCES (helper to simplify the update grid function)
POLAR DECOMPOSITION (needed for the math in the update grid function)*/
