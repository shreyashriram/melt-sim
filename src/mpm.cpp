

#define GLM_ENABLE_EXPERIMENTAL  
#include "mpm.h"
#include "particle.h"
#include "grid.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <random>  

float randomFloat(float min, float max) {
    static std::mt19937 rng(std::random_device{}()); // good random generator
    std::uniform_real_distribution<float> dist(min, max);
    return dist(rng);
}

void debugGroundImpact(const Particle& p) {
    std::string typeStr;
    switch (p.materialType) {
        case MaterialType::Solid:    typeStr = "Solid";    break;
        case MaterialType::Liquid:   typeStr = "Liquid";   break;
        case MaterialType::Melting:  typeStr = "Melting";  break;
        default:                     typeStr = "Unknown";  break;
    }

    std::cout << "Hit Ground | material: " << typeStr
              << " | velocity: (" << p.velocity.x << ", " << p.velocity.y << ", " << p.velocity.z << ")"
              << " | J: " << p.J << std::endl;
}

MPMSimulation::MPMSimulation() 
    : youngsModulus(1.0e4f),        // Higher stiffness for better stability
      poissonsRatio(0.4f),          // Less extreme for better stability
      grid(10, 0.25), 
      yieldThreshold(0.15f),        // Higher threshold for more stability
      meltRate(0.2f),               // More moderate melt rate
      globalMeltProgress(0.0f) {    // Start fully solid
    
    grid.setupBuffers();
    
    //calculate Lamé parameters from Young's modulus and Poisson's ratio
    shearModulus = youngsModulus / (2.0f * (1.0f + poissonsRatio));
    bulkModulus = youngsModulus * poissonsRatio / ((1.0f + poissonsRatio) * (1.0f - 2.0f * poissonsRatio));
}

void MPMSimulation::addMeshParticles(std::vector<Vector3> sampledPoints, MaterialType type) {
    Particle p;
    for (auto& pt : sampledPoints) {
        // Position the mesh higher for a longer fall
        
        glm::vec3 translate = glm::vec3(1.5f, 1.0f, 1.0f);
        glm::vec3 scale = glm::vec3(10.0f);


        glm::mat3 rot90Y = glm::mat3(glm::rotate(glm::mat4(1.0f), glm::radians(100.0f), glm::vec3(0, 1, 0)));

        glm::vec3 new_pos = rot90Y * (glm::vec3(pt.x(), pt.y(), pt.z()) * scale) + translate;


        p = Particle(new_pos, glm::vec3(0.0f, -1.0f, 0.0f));  // Initial downward velocity
        
        p.materialType = type;
        if (type == MaterialType::Solid) {
            p.meltStatus = 0.0f;
        } else if (type == MaterialType::Liquid) {
            p.meltStatus = 1.0f;
        } else if (type == MaterialType::Melting) {
            p.meltStatus = 0.0f; // will melt over time
        }
        
        particles.push_back(p);
    }
}

void MPMSimulation::spawnCube(MaterialType type, glm::vec3 center, float spacing, int countPerAxis) {
    particles.clear();

    float halfExtent = 0.5f * spacing * (countPerAxis - 1);

    for (int i = 0; i < countPerAxis; ++i) {
        for (int j = 0; j < countPerAxis; ++j) {
            for (int k = 0; k < countPerAxis; ++k) {
                glm::vec3 offset = glm::vec3(
                    (i * spacing) - halfExtent,
                    (j * spacing) - halfExtent,
                    (k * spacing) - halfExtent
                );
                Particle p(center + offset, glm::vec3(0.0f, -1.0f, 0.0f));
                p.materialType = type;
                if (type == MaterialType::Solid) {
                    p.meltStatus = 0.0f;
                } else if (type == MaterialType::Liquid) {
                    p.meltStatus = 1.0f;
                } else if (type == MaterialType::Melting) {
                    p.meltStatus = 0.0f; // will melt over time
                }
                particles.push_back(p);
            }
        }
    }
}

void MPMSimulation::step(float dt) {
    // Limit dt for stability
    float stableDt = glm::min(dt, 0.05f);
    
    // Progress the global melting over time (more controlled)
    globalMeltProgress += stableDt * 0.05f; // Even slower, more stable melting
    globalMeltProgress = glm::clamp(globalMeltProgress, 0.0f, 1.0f);

    // Reset grid
    for (auto& n : grid.nodes) {
        n.velocity = glm::vec3(0.0f);
        n.force = glm::vec3(0.0f);
        n.mass = 0.0f;
        // n.fluidFraction = 0.0f;  // Reset fluid fraction
    }
    
    transferParticlesToGrid();
    updateGrid(stableDt);
    transferGridToParticles(stableDt);
    updateParticles(stableDt);

    // resolveParticleCollisions(0.03f, 0.9f);
}

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
    
                    // APIC: Add the affine velocity term
                    glm::vec3 distToNode = nodePos - p.position;
                    glm::vec3 particleVelocity = p.velocity + p.C * distToNode;

                    // Dynamic damping based on melt status
                    float damping = 1.0f - (p.meltStatus * 0.02f);
                    
                    // Additional damping for crowded cells
                    if (nodeParticleCounts[linearIdx] > 4) {
                        damping *= 0.98f;
                    }
                    if (nodeParticleCounts[linearIdx] > 8) {
                        damping *= 0.98f;
                    }

                    grid.nodes[linearIdx].mass += p.mass * weight;
                    grid.nodes[linearIdx].velocity += p.mass * particleVelocity * weight * damping;                    
                }
            }
        }
    }
    
    // Normalize velocity by mass to get average 
    for (size_t i = 0; i < grid.nodes.size(); ++i) {
        if (grid.nodes[i].mass > 0.0f) {
            grid.nodes[i].velocity /= grid.nodes[i].mass;
            
            // Additional velocity smoothing based on crowding
            // Acts as a simple incompressibility constraint
            float damping = 1.0f;
            if (nodeParticleCounts[i] > 4) {
                damping = 0.98f;
            }
            if (nodeParticleCounts[i] > 8) {
                damping = 0.95f;
            }
            
            grid.nodes[i].velocity *= damping;
        }
    }
}

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
                    // Scale forces for stability and smooth transition
                    float stressScale = 1.0f - (p.meltStatus * 0.3f); // Less reduction for stability
                    

                    // & Force Damping!!
                    // Add force damping for stability
                    glm::vec3 force = -p.volume * (stress * weightGrad) * stressScale * 0.2f;
                    
                    int linearIdx = nodeIdx.x + nodeIdx.y * grid.size + nodeIdx.z * grid.size * grid.size;
                    grid.nodes[linearIdx].force += force;
                }
            }
        }
    }
    
    // Update grid velocities using forces with fluid-aware diffusion
    for (size_t idx = 0; idx < grid.nodes.size(); idx++) {
        auto& node = grid.nodes[idx];
        if (node.mass > 0.0f) {
            // Add reduced gravity for better stability
            node.force += glm::vec3(0.0f, -5.8f, 0.0f) * node.mass;
            
            // Apply velocity diffusion for fluid-like behavior
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
                // Moderate diffusion strength for stability
                float diffusionStrength = 0.96f;
                node.velocity = glm::mix(node.velocity, avgNeighborVel, diffusionStrength);
            }
            
            // Apply mild velocity damping for stability
            node.velocity *= 0.995f;
            
            // Update velocity with forces (with damping)
            node.velocity += (node.force / node.mass) * dt * 0.8f;
            
            // Clamp extreme velocities for stability
            float maxVel = 5.0f;
            if (glm::length(node.velocity) > maxVel) {
                node.velocity = glm::normalize(node.velocity) * maxVel;
            }
        }
    }
}

void MPMSimulation::transferGridToParticles(float dt) { 
    static bool firstUpdate = true;
    
    for (auto& p : particles) {
        glm::vec3 cellIdx = (p.position / grid.spacing);
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));
        
        glm::vec3 newVelocity(0.0f);
        glm::mat3 B(0.0f); // Temporary matrix for computing C
        
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
                    
                    // Compute affine matrix B for APIC
                    glm::vec3 distToNode = nodePos - p.position;
                    B += glm::mat3(
                        nodeVelocity.x * distToNode.x, nodeVelocity.x * distToNode.y, nodeVelocity.x * distToNode.z,
                        nodeVelocity.y * distToNode.x, nodeVelocity.y * distToNode.y, nodeVelocity.y * distToNode.z,
                        nodeVelocity.z * distToNode.x, nodeVelocity.z * distToNode.y, nodeVelocity.z * distToNode.z
                    ) * weight;
                }
            }
        }
        
        p.velocity = newVelocity;
        
        // Update C (affine velocity field)
        float scale = 4.0f / (grid.spacing * grid.spacing);
        p.C = B * scale;
        
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
        if (p.materialType == MaterialType::Solid) {
            updateSolidParticle(p, dt);
        } else if (p.materialType == MaterialType::Liquid) {
            updateLiquidParticle(p, dt);
        } else if (p.materialType == MaterialType::Melting) {
            updateMeltingParticle(p, dt);
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
    }

    grad /= grid.spacing;
    return grad;
}

void MPMSimulation::polarDecomposition(const glm::mat3& F, glm::mat3& R, glm::mat3& S){
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
    if (p.J <= 0.0f || p.J > 10.0f) {
        p.F = glm::mat3(1.0f);
        p.J = 1.0f;
        return;
    }
    
    // More conservative threshold for stability
    if (p.meltStatus > 0.5f) {
        // Preserve volumetric component (J) with safety clamp
        float J = glm::clamp(p.J, 0.5f, 2.0f);
        
        // Add very gentle randomness for fluid-like behavior
        float jitterAmount = 0.01f * p.meltStatus;
        float jitter = 1.0f + (glm::sin(p.position.x * 31.0f + p.position.y * 23.0f + p.position.z * 17.0f) * jitterAmount);
        J *= jitter;
        
        // Reset F to a pure volume-changing matrix (scaled identity)
        float scaleFactor = std::pow(J, 1.0f/3.0f);
        
        // More gradual blend towards fluid for stability
        glm::mat3 fluidF = glm::mat3(1.0f) * scaleFactor;
        float blendFactor = (p.meltStatus - 0.5f) * 2.0f; // 0 at meltStatus=0.5, 1 at meltStatus=1.0
        
        // Blend between current F and fluid F based on melt status
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                p.F[i][j] = glm::mix(p.F[i][j], fluidF[i][j], blendFactor * 0.3f);
            }
        }
        
        // If fully melted, just use the fluid F directly
        if (p.meltStatus > 0.9f) {
            p.F = fluidF;
        }
        
        // Recalculate J after the blend
        p.J = glm::determinant(p.F);
        return;
    }
    
    // For solid material, use a more stable plasticity model
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
    
    // More conservative yielding threshold for stability
    float effectiveYieldThreshold = yieldThreshold * (1.0f - 0.7f * p.meltStatus);
    
    // Check if yielding occurs - use more stable approach
    if (normDevS > effectiveYieldThreshold) {
        // Calculate how much to scale back the elastic deformation
        float scale = effectiveYieldThreshold / normDevS;
        
        // Apply more conservative hardening
        float hardeningFactor = 0.3f * (1.0f - p.meltStatus);
        scale = scale * hardeningFactor + (1.0f - hardeningFactor);
        
        // Create the modified stretch matrix
        glm::mat3 newS = (traceS / 3.0f) * identity + scale * devS;
        
        // Use more conservative deformation limits
        float minS = 0.3f + 0.2f * p.meltStatus;  // Higher min for stability
        float maxS = 1.1f + 0.6f * p.meltStatus;  // Lower max for stability
        
        for (int i = 0; i < 3; i++) {
            if (newS[i][i] > maxS) newS[i][i] = maxS;
            if (newS[i][i] < minS) newS[i][i] = minS;
        }
        
        // Reconstruct F with the yielded S
        p.F = R * newS;
        
        // Clamp determinant for stability
        p.J = glm::determinant(p.F);
        p.J = glm::clamp(p.J, 0.5f, 2.0f);
        
        // Rescale F to match clamped J if needed
        if (p.J != glm::determinant(p.F)) {
            float currentJ = glm::determinant(p.F);
            float scale = std::pow(p.J / currentJ, 1.0f/3.0f);
            p.F *= scale;
        }
    }
}

glm::mat3 MPMSimulation::computeStress(const Particle& p) {
    // Get rotation and stretch components
    glm::mat3 R, S;
    polarDecomposition(p.F, R, S);
    
    // More conservative stress calculation for stability
    
    // Scale shear modulus based on melt status, with less extreme reduction
    float effectiveShear = shearModulus * (1.0f - 0.9f * p.meltStatus);
    
    // Keep high bulk modulus but with slight reduction for melted state
    float effectiveBulk = bulkModulus * (1.0f - 0.1f * p.meltStatus);
    
    // Compute strain: E = F - R (linear strain for co-rotational model)
    glm::mat3 strain = p.F - R;
    
    // Compute trace of strain
    float trace = strain[0][0] + strain[1][1] + strain[2][2];
    
    // Safety clamp on trace for stability
    trace = glm::clamp(trace, -0.5f, 0.5f);
    
    // Compute deviatoric strain: E' = E - (1/3)*trace(E)*I
    glm::mat3 identity(1.0f);
    glm::mat3 strainDeviatoric = strain - (trace/3.0f) * identity;
    
    // Clamp extreme deviatoric strain components for stability
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            strainDeviatoric[i][j] = glm::clamp(strainDeviatoric[i][j], -0.2f, 0.2f);
        }
    }
    
    // Compute stress: σ = 2μ*E' + κ*trace(E)*I
    glm::mat3 elasticStress = 2.0f * effectiveShear * strainDeviatoric + 
                             effectiveBulk * trace * identity;
    
    // Apply global stress damping factor for stability
    return elasticStress * 0.8f;

    
}

void MPMSimulation::updateSolidParticle(Particle& p, float dt) {
    // TODO: fill this in next

    float maxVel = 12.0f;
    if (glm::length(p.velocity) > maxVel) {
        p.velocity = glm::normalize(p.velocity) * maxVel;
    }
    
    // Update position based on velocity with sub-stepping for stability
    int substeps = 3;
    float subDt = dt / float(substeps);
    for (int step = 0; step < substeps; step++) {
        p.position += p.velocity * subDt;
        
        // Apply boundary conditions in each substep
        // FLOOR COLLISIONS
        if (p.position.y < 0.0f) {
            
            p.position.y = 0.0f;
            
            // ** JITTERING 
            float splashJitter = 0.001f; // small randomness
            p.position.x += randomFloat(-splashJitter, splashJitter);
            p.position.z += randomFloat(-splashJitter, splashJitter);

            
            // ** BOUNCING
            float impactVelocity = std::abs(p.velocity.y);

            // 0.2 solid bounce -> 0.5 liquid bounce based on meltStatus
            float baseBounce = glm::mix(0.2f, 0.5f, p.meltStatus);
            
            // bounce varies from 0.05 (small impact) up to baseBounce (full impact)
            float velocityFactor = glm::clamp(impactVelocity / 15.0f, 0.0f, 1.0f);
            float bounce = glm::mix(0.05f, baseBounce, velocityFactor);
            
            p.velocity.y = impactVelocity * bounce;
        

            // ** DAMPING
            // Horizontal damping depends on melt status: [0.9: solid --> 0.98: liquid]
            float baseDamping = glm::mix(0.5f, 0.98f, p.meltStatus);
            p.velocity.x *= baseDamping;
            p.velocity.z *= baseDamping;


            debugGroundImpact(p);
        } 
    }

    updatePlasticity(p, dt);
}

void MPMSimulation::updateLiquidParticle(Particle& p, float dt) {
    // TODO: fill this in next

    float maxVel = 15.0f;
    if (glm::length(p.velocity) > maxVel) {
        p.velocity = glm::normalize(p.velocity) * maxVel;
    }
    
    // Update position based on velocity with sub-stepping for stability
    int substeps = 3;
    float subDt = dt / float(substeps);
    for (int step = 0; step < substeps; step++) {
        p.position += p.velocity * subDt;
        
        // Apply boundary conditions in each substep
        // FLOOR COLLISIONS
        if (p.position.y < 0.0f) {
            
            p.position.y = 0.0f;
            
            // ** JITTERING 
            float splashJitter = 0.015f; // small randomness
            p.position.x += randomFloat(-splashJitter, splashJitter);
            p.position.z += randomFloat(-splashJitter, splashJitter);

            
            // ** BOUNCING
            float impactVelocity = std::abs(p.velocity.y);

            float baseBounce = 0.05f;
            float velocityFactor = glm::clamp(impactVelocity / 15.0f, 0.0f, 1.0f);
            float bounce = glm::mix(0.02f, baseBounce, velocityFactor);
            
            p.velocity.y = impactVelocity * bounce;

            // ** SPLASHING
            if (impactVelocity > 1.0f) {
                float randX = glm::sin(p.position.x * 43.0f + p.position.z * 17.0f) * 0.5f + 0.5f;
                float randZ = glm::cos(p.position.x * 23.0f + p.position.z * 31.0f) * 0.5f + 0.5f;

                float splashMagnitude = impactVelocity * 0.25f;

                p.velocity.x += (randX * 2.0f - 1.0f) * splashMagnitude;
                p.velocity.z += (randZ * 2.0f - 1.0f) * splashMagnitude;


                // ** DAMPING
                // Horizontal damping depends on melt status: [0.9: solid --> 0.98: liquid]
                float baseDamping = 0.98f;
                p.velocity.x *= baseDamping;
                p.velocity.z *= baseDamping;
            } 

            debugGroundImpact(p);
        }
    }
}

void MPMSimulation::updateMeltingParticle(Particle& p, float dt) {
    
    // ** MELTING RATE
    // Update melt status
    float floorHeight = 0.0f;
    float fallHeight = 1.75f;
    float maxMeltDist = fallHeight - floorHeight;

    // How close to ground
    float heightFactor = 1.0f - glm::clamp((p.position.y - floorHeight) / maxMeltDist, 0.0f, 1.0f);

    // Combine height-based melting and global progress
    float targetMeltStatus = heightFactor * globalMeltProgress * 0.7f;

    // Slowly blend toward target
    float meltRate = 0.1f;
    p.meltStatus = glm::mix(p.meltStatus, targetMeltStatus, dt * meltRate);

    // ** 

    float maxVel = 6.0f;
    if (glm::length(p.velocity) > maxVel) {
        p.velocity = glm::normalize(p.velocity) * maxVel;
    }

    int substeps = 3;
    float subDt = dt / float(substeps);
    for (int step = 0; step < substeps; step++) {
        p.position += p.velocity * subDt;

        if (p.position.y < 0.0f) {
            p.position.y = 0.0f;

            // ** JITTERING 
            float splashJitter = glm::mix(0.008f, 0.030f, p.meltStatus); // blend between solid jitter and liquid jitter
            p.position.x += randomFloat(-splashJitter, splashJitter);
            p.position.z += randomFloat(-splashJitter, splashJitter);

            float impactVelocity = std::abs(p.velocity.y);

            // ** BOUNCING
            float solidBounce = 0.2f;
            float fluidBounce = 0.05f;
            float baseBounce = glm::mix(solidBounce, fluidBounce, p.meltStatus);

            float velocityFactor = glm::clamp(impactVelocity / 15.0f, 0.0f, 1.0f);
            float bounce = glm::mix(0.05f, baseBounce, velocityFactor);

            p.velocity.y = impactVelocity * bounce;

            // ** SPLASHING
            if (impactVelocity > 1.0f) {
                float randX = glm::sin(p.position.x * 43.0f + p.position.z * 17.0f) * 0.5f + 0.5f;
                float randZ = glm::cos(p.position.x * 23.0f + p.position.z * 31.0f) * 0.5f + 0.5f;

                float splashMagnitude = impactVelocity * glm::mix(0.15f, 0.25f, p.meltStatus); // bigger splashes when more melted
                p.velocity.x += (randX * 2.0f - 1.0f) * splashMagnitude;
                p.velocity.z += (randZ * 2.0f - 1.0f) * splashMagnitude;
            }

            // ** DAMPING
            float baseDamping = glm::mix(0.9f, 0.98f, p.meltStatus);
            p.velocity.x *= baseDamping;
            p.velocity.z *= baseDamping;

            debugGroundImpact(p);
        }
    }


    // && UPDATE PLASTICITY AT MELTING STATUS
    if (p.meltStatus < 0.9f) {
        // updatePlasticity(p, dt);
    }

    // CHANGE MATERIAL 
    if (p.meltStatus >= 0.95f) {
        p.materialType = MaterialType::Liquid;
    }
}
