#include "particle.h"

void Particle::initializeParticles() {
    particles.clear(); 
    
    for (int i = 0; i < 5; i++) {
        Particle p = Particle(glm::vec3(i*0.2f, 0.5f, 0.0f), glm::vec3(0.0f, i*0.05f, 0.0f));

        particles.push_back(p);
    }
}

void Particle::updateParticles(float dt) {
    for (auto& p : particles) {
        // Update velocity
        p.velocity += (p.force / p.mass) * dt;
        
        // Update position
        p.position += p.velocity * dt;
        
        // Modified boundary conditions to match grid size
        float gridBoundary = gridSize * gridSpacing;
        
        // Keep particles within [0, gridBoundary) instead of [-1, 1]
        for (int i = 0; i < 3; ++i) {
            if (p.position[i] < 0.0f) {
                p.position[i] = 0.0f;
                p.velocity[i] = 0.0f;
            }
            if (p.position[i] >= gridBoundary) {
                p.position[i] = gridBoundary - 0.001f; // Slightly inside boundary
                p.velocity[i] = 0.0f;
            }
        }
    }
}
