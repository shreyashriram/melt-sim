#include "particle.h"
#include "grid.h"

void Particles::initializeParticles() {
    particles.clear(); 
    
    for (int i = 0; i < 5; i++) {
        Particle p = Particle(glm::vec3(i*0.2f, 0.5f, 0.0f), glm::vec3(0.0f, i*0.05f, 0.0f));

        particles.push_back(p);
    }
}


