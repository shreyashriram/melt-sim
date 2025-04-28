#define GLM_ENABLE_EXPERIMENTAL  // Add this line before any GLM includes
#include "mpm.h"
#include "particle.h"
#include "grid.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>


bool alreadyPrinted = false; //for debugging 

MPMSimulation::MPMSimulation() 
    : youngsModulus(1.0e5f),  // Much lower stiffness for liquid-like behavior
      poissonsRatio(0.49f),    // Very close to 0.5 for incompressible fluid
      grid(5, 0.5f), 
      yieldThreshold(0.01f),   // Very low yield threshold for fluid-like flow
      meltRate(2.0f),          // Increased melt rate for more noticeable changes
      globalMeltProgress(0.0f) {
    
    grid.setupBuffers();
    
    //calculate Lamé parameters from Young's modulus and Poisson's ratio
    shearModulus = youngsModulus / (2.0f * (1.0f + poissonsRatio));
    bulkModulus = youngsModulus * poissonsRatio / ((1.0f + poissonsRatio) * (1.0f - 2.0f * poissonsRatio));
    
    if (poissonsRatio > 0.495f) {
        bulkModulus = 1000.0f * shearModulus;
    }
}


// TODO: make sure its in grid, currently manually translated
void MPMSimulation::addMeshParticles(std::vector<Vector3> sampledPoints) {
    Particle p;
    for (auto& pt : sampledPoints) {
        // Give initial downward velocity and slight horizontal motion
        float mesh_translate = 1.25f;

        p = Particle(glm::vec3(pt.x()+0.75, pt.y()+mesh_translate, pt.z()+0.75), 
                    glm::vec3(0.0f, -2.0f, 0.0f));  // Moderate initial velocity
        p.meltStatus = 0.0f;  // Start as solid
        particles.push_back(p);
    }
}

// ? Debugging: Manually create particles
void MPMSimulation::initializeParticles() {
    particles.clear(); 
    
    Particle p;  //declare the Particle variable
    for (int i = 0; i < 5; i++) {
        p = Particle(glm::vec3(i*0.2f, 0.5f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f));

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
        // n.fluidFraction = 0.0f;  // Reset fluid fraction
    }
    
    transferParticlesToGrid();
    updateGrid(dt);
    transferGridToParticles(dt);
    updateParticles(dt);
}

// STEP 1 IN MPM GUIDE
void MPMSimulation::transferParticlesToGrid() { 
    // Add a local counter for each grid cell
    std::vector<int> nodeParticleCounts(grid.nodes.size(), 0);

    for (const auto& p : particles) {
        
        
        glm::vec3 cellIdx = (p.position / grid.spacing);
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));

        
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    glm::ivec3 nodeIdx = baseNode + glm::ivec3(i, j, k);
                    
                    if (nodeIdx.x < 0 || nodeIdx.x >= grid.size || 
                        nodeIdx.y < 0 || nodeIdx.y >= grid.size || 
                        nodeIdx.z < 0 || nodeIdx.z >= grid.size) {
                        continue;
                    }

                    int linearIdx = nodeIdx.x + (nodeIdx.y * grid.size) + (nodeIdx.z * grid.size * grid.size);
                    nodeParticleCounts[linearIdx]++;
                
                    glm::vec3 nodePos = (glm::vec3(nodeIdx) * grid.spacing);
                    float weight = computeWeight(p.position, nodePos);
    
                    // ! DAMPING 
                    // Minimal damping based on particle count
                    float damping = 1.0f;
                    if (nodeParticleCounts[linearIdx] > 4) {
                        damping = 0.99f;
                    }
                    if (nodeParticleCounts[linearIdx] > 8) {
                        damping = 0.98f;
                    }

                    grid.nodes[linearIdx].mass += p.mass * weight;
                    grid.nodes[linearIdx].velocity += p.mass * p.velocity * weight * damping;                    
                }
            }
        }
    }
    
    // Normalize velocity by mass to get average 
    for (size_t i = 0; i < grid.nodes.size(); ++i) {
        if (grid.nodes[i].mass > 0.0f) {
            grid.nodes[i].velocity /= grid.nodes[i].mass;
            
            // Normalize fluid fraction by mass
            // grid.nodes[i].fluidFraction /= grid.nodes[i].mass;
            
            // ! DAMPING 
            // Minimal damping for overcrowded nodes
            if (nodeParticleCounts[i] > 4) {
                grid.nodes[i].velocity *= 0.99f;
            }
            if (nodeParticleCounts[i] > 8) {
                grid.nodes[i].velocity *= 0.98f;
            }
        }
    }
}

// STEP 2 IN THE MPM GUIDE
void MPMSimulation::updateGrid(float dt) {
    // First compute forces from stress for each particle
    for (const auto& p : particles) {
        glm::mat3 stress = computeStress(p);
        glm::vec3 cellIdx = p.position / grid.spacing;
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));
        
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    glm::ivec3 nodeIdx = baseNode + glm::ivec3(i, j, k);
                    
                    if (nodeIdx.x < 0 || nodeIdx.x >= grid.size || 
                        nodeIdx.y < 0 || nodeIdx.y >= grid.size || 
                        nodeIdx.z < 0 || nodeIdx.z >= grid.size) {
                        continue;
                    }
                    
                    glm::vec3 nodePos = glm::vec3(nodeIdx) * grid.spacing;
                    glm::vec3 weightGrad = computeWeightGradient(p.position, nodePos);
                    
                    // Force = -volume * stress * weightGrad
                    glm::vec3 force = -p.volume * (stress * weightGrad);
                    
                    // ! Scale down forces for stability
                    force *= 0.5f;
                    
                    int linearIdx = nodeIdx.x + nodeIdx.y * grid.size + nodeIdx.z * grid.size * grid.size;
                    grid.nodes[linearIdx].force += force;
                }
            }
        }
    }
    
    // Then update grid velocities using forces with fluid-aware diffusion
    for (size_t idx = 0; idx < grid.nodes.size(); idx++) {
        auto& node = grid.nodes[idx];
        if (node.mass > 0.0f) {
            // // Add gravity
            node.force += glm::vec3(0.0f, -9.8f, 0.0f) * node.mass;
            
            // // Apply velocity diffusion for fluid-like behavior
            // // This simulates viscosity in fluid regions
            
            glm::vec3 avgNeighborVel(0.0f);
            float totalWeight = 0.0f;
            
            // Extract 3D indices from linear index
            int z = idx / (grid.size * grid.size);
            int remainder = idx % (grid.size * grid.size);
            int y = remainder / grid.size;
            int x = remainder % grid.size;
            glm::ivec3 nodeIdx(x, y, z);
            
            // Look at neighboring grid nodes
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    for (int k = -1; k <= 1; ++k) {
                        if (i == 0 && j == 0 && k == 0) continue; // Skip self
                        
                        glm::ivec3 neighborIdx = nodeIdx + glm::ivec3(i, j, k);
                        if (neighborIdx.x < 0 || neighborIdx.x >= grid.size || 
                            neighborIdx.y < 0 || neighborIdx.y >= grid.size || 
                            neighborIdx.z < 0 || neighborIdx.z >= grid.size) {
                            continue;
                        }
                        
                        int linearIdx = neighborIdx.x + neighborIdx.y * grid.size + 
                                        neighborIdx.z * grid.size * grid.size;
                                        
                        if (grid.nodes[linearIdx].mass > 0.0f) {
                            avgNeighborVel += grid.nodes[linearIdx].velocity;
                            totalWeight += 1.0f;
                        }
                    }
                }
            }
            
            if (totalWeight > 0.0f) {
                avgNeighborVel /= totalWeight;
                // Blend velocity with neighbors based on fluid fraction
                // This creates viscosity-like behavior
                float diffusionStrength = 0.5f;
                node.velocity = glm::mix(node.velocity, avgNeighborVel, diffusionStrength);
            }
    
    
            // Apply minimal damping and update velocity
            // node.velocity *= 0.999f;  // Minimal damping
            node.velocity += (node.force / node.mass) * dt;
        }
    }
}

// STEP 3 IN MPM GUIDE (basically same as P2G but in opposite direction?)
void MPMSimulation::transferGridToParticles(float dt) { 
    static bool firstUpdate = true;
    
    for (auto& p : particles) {
        glm::vec3 cellIdx = (p.position / grid.spacing);
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));
        
        glm::vec3 newVelocity(0.0f);
        
        // Compute PIC velocity transfer
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    glm::ivec3 nodeIdx = baseNode + glm::ivec3(i, j, k);
                    if (nodeIdx.x < 0 || nodeIdx.x >= grid.size || 
                        nodeIdx.y < 0 || nodeIdx.y >= grid.size ||
                        nodeIdx.z < 0 || nodeIdx.z >= grid.size) {
                        continue;
                    }

                    int linearIdx = nodeIdx.x + nodeIdx.y * grid.size + nodeIdx.z * grid.size * grid.size;
                    glm::vec3 nodePos = glm::vec3(nodeIdx) * grid.spacing;
                    
                    float weight = computeWeight(p.position, nodePos);
                    glm::vec3 nodeVelocity = grid.nodes[linearIdx].velocity;
                    
                    // PIC update
                    newVelocity += nodeVelocity * weight;
                }
            }
        }
        
        p.velocity = newVelocity;
        p.velocityGradient = computeVelocityGradient(p);
        
        // Update deformation gradient
        glm::mat3 F_new = (glm::mat3(1.0f) + dt * p.velocityGradient) * p.F;
        p.F = F_new;
        p.J = glm::determinant(p.F);
    }
}

void MPMSimulation::updateParticles(float dt) {
    float gridBoundary = grid.size * grid.spacing;
    
    for (auto& p : particles) {
        // Position-based melting
        // Define a transition height - above this height particles are solid,
        // below this height they are liquid
        float transitionHeight = 0.05f; // Lowered from 1.0f to 0.2f
        float transitionWidth = 0.2f;  // Controls how gradual the transition is
        
        // Calculate melt status based on height (y-position)
        // 1.0 = fully liquid, 0.0 = fully solid
        if (p.position.y < transitionHeight - transitionWidth) {
            // Below transition zone - fully liquid
            p.meltStatus = 1.0f;
        } 
        else if (p.position.y > transitionHeight + transitionWidth) {
            // Above transition zone - fully solid
            p.meltStatus = 0.0f;
        }
        else {
            // In transition zone - gradual change from liquid to solid
            float t = (p.position.y - (transitionHeight - transitionWidth)) / (2.0f * transitionWidth);
            p.meltStatus = 1.0f - t; // Smoothly transition from 1.0 to 0.0
        }
        
        // Update position based on velocity
        p.position += p.velocity * dt;
        
        // Apply plasticity to limit deformation
        updatePlasticity(p, dt);
        
        // Floor collision with melt-dependent behavior
        if (p.position.y < 0.0f) {
            p.position.y = 0.0f;
            
            // Calculate impact velocity
            float impactVelocity = std::abs(p.velocity.y);
            
            // Adjust bounce based on melt status
            float solidBounce = 0.8f;
            float liquidBounce = 0.3f;
            float bounce = glm::mix(solidBounce, liquidBounce, p.meltStatus);
            
            // Apply bounce to vertical velocity
            p.velocity.y = impactVelocity * bounce;
            
            // Adjust horizontal damping based on melt status
            float solidDamping = 0.99f;
            float liquidDamping = 0.95f;
            float damping = glm::mix(solidDamping, liquidDamping, p.meltStatus);
            
            p.velocity.x *= damping;
            p.velocity.z *= damping;
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

    /*std::cout << "=== Weight Gradient Debug ===" << std::endl;
    std::cout << "Particle Position: (" << particlePos.x << ", " << particlePos.y << ", " << particlePos.z << ")" << std::endl;
    std::cout << "Node Position: (" << nodePos.x << ", " << nodePos.y << ", " << nodePos.z << ")" << std::endl;
    std::cout << "Difference (rel): (" << rel.x << ", " << rel.y << ", " << rel.z << ")" << std::endl;
*/
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

        //std::cout << "dim " << dim << ": x=" << x << ", gradW=" << gradW << ", weight[" << a << "]=" << wa << ", weight[" << b << "]=" << wb << std::endl;
        //std::cout << "Partial Gradient[" << dim << "] = " << grad[dim] << std::endl;
    }

    grad /= grid.spacing;
    //std::cout << "Final Gradient: (" << grad.x << ", " << grad.y << ", " << grad.z << ")" << std::endl;

    return grad;
}

void MPMSimulation::polarDecomposition(const glm::mat3& F, glm::mat3& R, glm::mat3& S){
    /*static bool firstDecomp = true;
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
    }*/
    
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

   /*if (firstDecomp) {
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
    }*/

}

glm::mat3 MPMSimulation::computeVelocityGradient(const Particle& p) {
    glm::vec3 cellIdx = (p.position / grid.spacing);
    glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));
    glm::mat3 velocityGradient(0.0f);
    
    // Loop through neighboring nodes
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            for (int k = -1; k <= 1; ++k) {
                glm::ivec3 nodeIdx = baseNode + glm::ivec3(i, j, k);
                
                // Skip if outside grid
                if (nodeIdx.x < 0 || nodeIdx.x >= grid.size || 
                    nodeIdx.y < 0 || nodeIdx.y >= grid.size ||
                    nodeIdx.z < 0 || nodeIdx.z >= grid.size) {
                    continue;
                }

                int linearIdx = nodeIdx.x + nodeIdx.y * grid.size + nodeIdx.z * grid.size * grid.size;
                glm::vec3 nodePos = glm::vec3(nodeIdx) * grid.spacing;
                
                glm::vec3 weightGrad = computeWeightGradient(p.position, nodePos);
                glm::vec3 nodeVelocity = grid.nodes[linearIdx].velocity;
                
                // Compute velocity gradient
                for (int d = 0; d < 3; d++) {
                    for (int g = 0; g < 3; g++) {
                        velocityGradient[g][d] += nodeVelocity[d] * weightGrad[g];
                    }
                }
            }
        }
    }
    
    return velocityGradient;
}

void MPMSimulation::updatePlasticity(Particle& p, float dt) {
    // Skip if deformation is invalid
    if (p.J <= 0.0f) {
        p.F = glm::mat3(1.0f);
        p.J = 1.0f;
        return;
    }
    
    // For highly melted material, reset the deviatoric part of F
    // This is key to fluid-like behavior - fluids don't remember shear deformation
    if (p.meltStatus > 0.8f) {
        // Preserve only the volumetric component (J)
        float J = p.J;
        
        // Reset F to a pure volume-changing matrix (scaled identity)
        // This cubic root ensures volume is preserved exactly
        float scaleFactor = std::pow(J, 1.0f/3.0f);
        p.F = glm::mat3(1.0f) * scaleFactor;
        return;
    }
    
    // For partially melted or solid material, use existing plasticity model
    // but with yield threshold that scales with melt status
    glm::mat3 R, S;
    polarDecomposition(p.F, R, S);
    
    // Calculate von Mises equivalent strain
    float traceS = S[0][0] + S[1][1] + S[2][2];
    glm::mat3 identity(1.0f);
    glm::mat3 devS = S - (traceS / 3.0f) * identity;
    
    // Compute Frobenius norm of deviatoric stretch
    float normDevS = 0.0f;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            normDevS += devS[i][j] * devS[i][j];
        }
    }
    normDevS = sqrt(normDevS);
    
    // Yield threshold decreases with melt status
    float effectiveYieldThreshold = yieldThreshold * (1.0f - 0.9f * p.meltStatus);
    
    // Check if yielding occurs
    if (normDevS > effectiveYieldThreshold) {
        // Calculate how much to scale back the elastic deformation
        float scale = effectiveYieldThreshold / normDevS;
        
        // Apply hardening that decreases with melt status
        float hardeningFactor = 0.5f * (1.0f - p.meltStatus);
        scale = scale * hardeningFactor + (1.0f - hardeningFactor);
        
        // Create the modified stretch matrix
        glm::mat3 newS = (traceS / 3.0f) * identity + scale * devS;
        
        // Allow more deformation for melted material
        float minS = 0.2f + 0.3f * p.meltStatus;  // Lower bound decreases with melt
        float maxS = 1.5f + 1.0f * p.meltStatus;  // Upper bound increases with melt
        
        for (int i = 0; i < 3; i++) {
            if (newS[i][i] > maxS) newS[i][i] = maxS;
            if (newS[i][i] < minS) newS[i][i] = minS;
        }
        
        // Reconstruct F with the yielded S
        p.F = R * newS;
        p.J = glm::determinant(p.F);
    }
}

glm::mat3 MPMSimulation::computeStress(const Particle& p) {
    // // Check for numerical instability in deformation gradient
    // float maxF = 0.0f;
    // for (int i = 0; i < 3; i++) {
    //     for (int j = 0; j < 3; j++) {
    //         maxF = std::max(maxF, std::abs(p.F[i][j]));
    //     }
    // }
    
    // // If F has grown too large, apply stabilization
    // if (maxF > 1000.0f) {
    //     // Create a temporary copy with stabilized F
    //     glm::mat3 stabilizedF;
    //     float scaleFactor = 1.0f / maxF;
        
    //     for (int i = 0; i < 3; i++) {
    //         for (int j = 0; j < 3; j++) {
    //             stabilizedF[i][j] = p.F[i][j] * scaleFactor;
    //         }
    //     }
        
    //     // Use the stabilized F for stress computation
    //     glm::mat3 R, S;
    //     polarDecomposition(stabilizedF, R, S);
        
    //     // Scale shear modulus based on melt status
    //     // For fluid-like behavior, reduce shear resistance as melt increases
    //     float effectiveShear = shearModulus * (1.0f - 0.99f * p.meltStatus);
        
    //     // Keep bulk modulus high for volume preservation
    //     float effectiveBulk = bulkModulus;
        
    //     // Compute strain: E = F - R (linear strain for co-rotational model)
    //     glm::mat3 strain = stabilizedF - R;
        
    //     // Compute trace of strain
    //     float trace = strain[0][0] + strain[1][1] + strain[2][2];
        
    //     // Compute deviatoric strain: E' = E - (1/3)*trace(E)*I
    //     glm::mat3 identity(1.0f);
    //     glm::mat3 strainDeviatoric = strain - (trace/3.0f) * identity;
        
    //     // Compute stabilized stress: σ = 2μ*E' + κ*trace(E)*I
    //     glm::mat3 elasticStress = 2.0f * effectiveShear * strainDeviatoric + 
    //                              effectiveBulk * trace * identity;
        
    //     // Scale the stress back up to compensate for the F scaling
    //     return elasticStress;
    // }
    
    // Normal stress computation for stable cases
    glm::mat3 R, S;
    polarDecomposition(p.F, R, S);
    
    // Scale shear modulus based on melt status
    // When fully melted (meltStatus = 1.0), shear resistance approaches zero
    float effectiveShear = shearModulus * (0.1f * p.meltStatus);
    
    // Keep bulk modulus high to maintain volume preservation
    float effectiveBulk = bulkModulus;
    
    // Compute strain: E = F - R (linear strain for co-rotational model)
    glm::mat3 strain = p.F - R;
    
    // Compute trace of strain
    float trace = strain[0][0] + strain[1][1] + strain[2][2];
    
    // Compute deviatoric strain: E' = E - (1/3)*trace(E)*I
    glm::mat3 identity(1.0f);
    // glm::mat3 strainDeviatoric = strain - (trace/3.0f) * identity;
    
    // Compute stress: σ = 2μ*E' + κ*trace(E)*I
    glm::mat3 elasticStress = 2.0f * effectiveShear * strain + 
                             effectiveBulk * trace * identity;
    
    return elasticStress;
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

void MPMSimulation::testComputeStress() {
    std::cout << "=== Testing computeStress function ===\n";
    
    // Create a test particle with known F
    Particle testParticle;
    
    // Test 1: Identity deformation (no stress)
    testParticle.F = glm::mat3(1.0f); // Identity matrix
    testParticle.meltStatus = 0.0f;
    
    std::cout << "Test 1: Identity deformation\n";
    std::cout << "Input F = Identity matrix\n";
    
    glm::mat3 stress1 = computeStress(testParticle);
    
    std::cout << "Resulting stress:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << stress1[i][j] << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "Expected: All zeros (or very small values)\n\n";
    
    // Test 2: Simple stretch
    testParticle.F = glm::mat3(1.0f);
    testParticle.F[0][0] = 1.1f; // 10% stretch in x direction
    
    std::cout << "Test 2: 10% stretch in x direction\n";
    std::cout << "Input F:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << testParticle.F[i][j] << " ";
        }
        std::cout << "]\n";
    }
    
    glm::mat3 stress2 = computeStress(testParticle);
    
    std::cout << "Resulting stress:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << stress2[i][j] << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "Expected: Positive stress in x direction\n\n";
    
    // Test 3: Simple shear
    testParticle.F = glm::mat3(1.0f);
    testParticle.F[0][1] = 0.1f; // Shear in xy plane
    
    std::cout << "Test 3: Shear in xy plane\n";
    std::cout << "Input F:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << testParticle.F[i][j] << " ";
        }
        std::cout << "]\n";
    }
    
    glm::mat3 stress3 = computeStress(testParticle);
    
    std::cout << "Resulting stress:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << stress3[i][j] << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "Expected: Shear stress in xy plane\n\n";
    
    // Test 4: Large deformation
    testParticle.F = glm::mat3(1.0f);
    testParticle.F[0][0] = 2.0f; // Large stretch in x direction
    
    std::cout << "Test 4: Large stretch (100% stretch in x direction)\n";
    std::cout << "Input F:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << testParticle.F[i][j] << " ";
        }
        std::cout << "]\n";
    }
    
    glm::mat3 stress4 = computeStress(testParticle);
    
    std::cout << "Resulting stress:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << stress4[i][j] << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "Expected: Large positive stress in x direction\n\n";
    
    // Test 5: Rotation (should give zero stress in co-rotational model)
    float angle = glm::radians(45.0f);
    float c = cos(angle);
    float s = sin(angle);
    testParticle.F = glm::mat3(
        c, -s, 0.0f,
        s,  c, 0.0f,
        0.0f, 0.0f, 1.0f
    );
    
    std::cout << "Test 5: Pure rotation (45 degrees)\n";
    std::cout << "Input F:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << testParticle.F[i][j] << " ";
        }
        std::cout << "]\n";
    }
    
    glm::mat3 stress5 = computeStress(testParticle);
    
    std::cout << "Resulting stress:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << stress5[i][j] << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "Expected: All zeros (no stress from pure rotation)\n\n";
    
    // Test 6: Extremely large deformation
    testParticle.F = glm::mat3(1000.0f);
    
    std::cout << "Test 6: Extremely large deformation\n";
    std::cout << "Input F: All elements = 1000.0\n";
    
    glm::mat3 stress6 = computeStress(testParticle);
    
    std::cout << "Resulting stress:\n";
    for (int i = 0; i < 3; i++) {
        std::cout << "[ ";
        for (int j = 0; j < 3; j++) {
            std::cout << stress6[i][j] << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "Expected: Large but finite values (should not be NaN or inf)\n";
    
    std::cout << "=== End of computeStress testing ===\n";
}

/*FUNCTIONS TO ADD:
COMPUTE STRESS
UPDATE DEFORMATION GRADIENT 
APPLY PLASTICITY
COMPUTE INTERNAL FORCES (helper to simplify the update grid function)
POLAR DECOMPOSITION (needed for the math in the update grid function)*/
