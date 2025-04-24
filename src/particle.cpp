#include "particle.h"
#include "grid.h"

void Particles::initializeParticles() {
    particles.clear(); 
    
    for (int i = 0; i < 5; i++) {
        Particle p = Particle(glm::vec3(i*0.2f, 0.5f, 0.0f), glm::vec3(0.0f, i*0.05f, 0.0f));

        particles.push_back(p);
    }
}

void Particles::updateParticles(float dt, int gridSize, float gridSpacing) {
    float gridBoundary = gridSize * gridSpacing;
    
    for (auto& p : particles) {
        // Update velocity based on forces
        p.velocity += (p.force / p.mass) * dt;
        
        // Update position
        p.position += p.velocity * dt;
        
        // Simple boundary conditions
        for (int i = 0; i < 3; i++) {
            if (p.position[i] < 0.0f) {
                p.position[i] = 0.0f;
                p.velocity[i] = 0.0f;
            }
            else if (p.position[i] >= gridBoundary) {
                p.position[i] = gridBoundary - 0.001f; // Slightly inside boundary
                p.velocity[i] = 0.0f;
            }
        }
    }
}
