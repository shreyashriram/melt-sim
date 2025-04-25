#include "mpm.h"
#include "particle.h"
#include "grid.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>

//------------------------------------------------------------------------------
// Global debug flag
bool alreadyPrinted = false; // for debugging: toggles one‐time prints

//------------------------------------------------------------------------------
// Constructor: allocate and zero out the grid
MPMSimulation::MPMSimulation()
    : youngsModulus(1.4e5f), poissonsRatio(0.2f), gridSize(64), gridSpacing(0.25f)
{
    // allocate flat array of gridSize³ nodes
    grid.resize(gridSize * gridSize * gridSize);
    resetGrid();
}

//------------------------------------------------------------------------------
// resetGrid: clear mass, velocity, and force on every node
void MPMSimulation::resetGrid() {
    for (auto &n : grid) {
        n.mass     = 0.0f;
        n.velocity = glm::vec3(0.0f);
        n.force    = glm::vec3(0.0f);
    }
}

//------------------------------------------------------------------------------
// addMeshParticles: seed particles from a mesh sample
void MPMSimulation::addMeshParticles(std::vector<Vector3> sampledPoints) {
    for (auto &pt : sampledPoints) {
        // position offset by +1 in y to sit above floor
        Particle p(glm::vec3(pt.x(), pt.y() + 1.0f, pt.z()), glm::vec3(0.0f));
        particles.push_back(p);
    }
}

//------------------------------------------------------------------------------
// initializeParticles: build a small test line of particles
void MPMSimulation::initializeParticles() {
    particles.clear();
    for (int i = 0; i < 5; ++i) {
        Particle p(
            glm::vec3(i * 0.2f, 0.5f, 0.0f),  // start positions
            glm::vec3(0.0f, i * 0.05f, 0.0f)); // initial upward velocities
        particles.push_back(p);
    }
}

//------------------------------------------------------------------------------
// step: full MPM pipeline in one call
void MPMSimulation::step(float dt) {
    resetGrid();               // 1) clear grid from last frame
    transferParticlesToGrid();// 2) P2G: scatter particle mass+momentum
    updateGrid(dt);           // 3) apply forces & integrate on grid
    transferGridToParticles();// 4) G2P: gather back new velocities + affine
    updateParticles(dt, gridSize, gridSpacing); // 5) advect & collide
}

//------------------------------------------------------------------------------
// transferParticlesToGrid: Particle→Grid (P2G) using APIC
void MPMSimulation::transferParticlesToGrid() {
    float halfSize = gridSize * gridSpacing * 0.5f;
    for (auto &p : particles) {
        glm::vec3 offsetPos = p.position + glm::vec3(halfSize);
        glm::vec3 cellIdx   = offsetPos / gridSpacing - 0.5f;
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));
        for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
            glm::ivec3 nodeIdx = baseNode + glm::ivec3(dx, dy, dz);
            if (nodeIdx.x < 0 || nodeIdx.x >= gridSize ||
                nodeIdx.y < 0 || nodeIdx.y >= gridSize ||
                nodeIdx.z < 0 || nodeIdx.z >= gridSize) continue;
            int idx = nodeIdx.x + nodeIdx.y * gridSize + nodeIdx.z * gridSize * gridSize;
            glm::vec3 nodePos = glm::vec3(nodeIdx) * gridSpacing - glm::vec3(halfSize);
            float    w       = computeWeight(p.position, nodePos);
            glm::vec3 d      = nodePos - p.position;
            grid[idx].mass     += p.mass * w;
            grid[idx].velocity += p.mass * (p.velocity + p.C * d) * w;
        }
    }
    for (auto &n : grid)
        if (n.mass > 0.0f)
            n.velocity /= n.mass;
}

//------------------------------------------------------------------------------
// updateGrid: apply forces (gravity) and integrate velocity
void MPMSimulation::updateGrid(float dt) {
    for (auto &n : grid) {
        if (n.mass > 0.0f) {
            n.force += glm::vec3(0.0f, -9.8f, 0.0f) * n.mass;
            n.velocity += (n.force / n.mass) * dt;
        }
    }
}

//------------------------------------------------------------------------------
// transferGridToParticles: Grid→Particle (G2P) using APIC
void MPMSimulation::transferGridToParticles() {
    float halfSize = gridSize * gridSpacing * 0.5f;
    for (auto &p : particles) {
        glm::vec3 newV(0.0f);
        glm::mat3 newC(0.0f);
        glm::vec3 offsetPos = p.position + glm::vec3(halfSize);
        glm::vec3 cellIdx   = offsetPos / gridSpacing - 0.5f;
        glm::ivec3 baseNode = glm::ivec3(floor(cellIdx));
        for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
            glm::ivec3 nodeIdx = baseNode + glm::ivec3(dx, dy, dz);
            if (nodeIdx.x < 0 || nodeIdx.x >= gridSize ||
                nodeIdx.y < 0 || nodeIdx.y >= gridSize ||
                nodeIdx.z < 0 || nodeIdx.z >= gridSize) continue;
            int idx = nodeIdx.x + nodeIdx.y * gridSize + nodeIdx.z * gridSize * gridSize;
            glm::vec3 nodePos = glm::vec3(nodeIdx) * gridSpacing - glm::vec3(halfSize);
            float    w       = computeWeight(p.position, nodePos);
            glm::vec3 vi     = grid[idx].velocity;
            glm::vec3 d      = nodePos - p.position;
            newV += w * vi;
            newC += w * glm::outerProduct(vi, d);
        }
        p.velocity = newV;
        p.C = (4.0f / (gridSpacing * gridSpacing)) * newC;
    }
}

//------------------------------------------------------------------------------
// updateParticles: explicit integration, floor bounce, and simple collisions
void MPMSimulation::updateParticles(float dt, int gridSize, float gridSpacing) {
    // 1) gravity & advect
    for (auto &p : particles) {
        p.velocity += glm::vec3(0.0f, -9.81f, 0.0f) * dt;
        p.position += p.velocity * dt;
    }

    // 2) floor bounce
    for (auto &p : particles) {
        if (p.position.y < 0.0f) {
            p.position.y = 0.0f;
            p.velocity.y *= -0.5f; // restitution
        }
    }

    // 3) inter-particle collisions
    float collisionRadius = gridSpacing * 0.5f;
    float rad2 = collisionRadius * collisionRadius;
    int n = (int)particles.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            auto &p1 = particles[i];
            auto &p2 = particles[j];
            glm::vec3 d = p1.position - p2.position;
            float dist2 = glm::dot(d,d);
            if (dist2 > 0.0f && dist2 < rad2) {
                float dist = std::sqrt(dist2);
                glm::vec3 nrm = d / dist;
                float overlap = collisionRadius - dist;
                // separate
                p1.position += 0.5f * overlap * nrm;
                p2.position -= 0.5f * overlap * nrm;
                // 1D elastic on normal
                float v1n = glm::dot(p1.velocity, nrm);
                float v2n = glm::dot(p2.velocity, nrm);
                float m1 = p1.mass, m2 = p2.mass;
                float vn1 = (v1n*(m1-m2) + 2.0f*m2*v2n)/(m1+m2);
                float vn2 = (v2n*(m2-m1) + 2.0f*m1*v1n)/(m1+m2);
                p1.velocity += (vn1 - v1n) * nrm;
                p2.velocity += (vn2 - v2n) * nrm;
            }
        }
    }
}

//------------------------------------------------------------------------------
// computeWeight: trilinear hat basis
float MPMSimulation::computeWeight(const glm::vec3 &particlePos,
                                   const glm::vec3 &nodePos) {
    glm::vec3 d = (particlePos - nodePos) / gridSpacing;
    float wx = std::max(0.0f, 1.0f - std::abs(d.x));
    float wy = std::max(0.0f, 1.0f - std::abs(d.y));
    float wz = std::max(0.0f, 1.0f - std::abs(d.z));
    return wx * wy * wz;
}

//------------------------------------------------------------------------------
// computeWeightGradient: placeholder for B-spline / MLS
glm::vec3 MPMSimulation::computeWeightGradient(const glm::vec3 &particlePos,
                                               const glm::vec3 &nodePos) {
    return glm::vec3(0.0f);
}
