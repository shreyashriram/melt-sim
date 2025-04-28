
#pragma once
#include <vector>
#include <glm/glm.hpp>

struct Particle {
    glm::vec3 position;
    glm::vec3 velocity;
    glm::vec3 force;
    float mass;


    glm::mat3 F; //deformation gradient
    glm::mat3 F_p; //plasticity deformation gradient
    glm::mat3 F_e; //elasticity deformation gradient
    float meltStatus; //a number between 0 and 1 indicating the amount of melting
    float volume;
    float J; //determinant of the deformation gradient
    glm::mat3 C; //APIC velocity field matrix
    glm::mat3 velocityGradient; //velocity gradient tensor

    
    Particle(float density = 1000.0f) : position(0.0f), velocity(0.0f), force(0.0f), mass(1.0f), J(1.0f), volume(1.0f), meltStatus(0.0f),F(1.0f), F_p(1.0f), F_e(1.0f), C(0.0f), velocityGradient(0.0f) {
        volume = mass / density;
        J = 1.0f;
    };
          
    //Constructor with position
    Particle(const glm::vec3& pos, float density = 1000.0f) : position(pos), velocity(0.0f), force(0.0f), mass(1.0f), J(1.0f), volume(1.0f), meltStatus(0.0f),F(1.0f), F_p(1.0f), F_e(1.0f), C(0.0f), velocityGradient(0.0f)  {
        volume = mass / density;
        J = 1.0f;
    };
    
    // Constructor with position and velocity   
    Particle(const glm::vec3& pos, const glm::vec3& vel, float density = 1000.0f) : position(pos), velocity(vel), force(0.0f), mass(1.0f), J(1.0f), volume(1.0f), meltStatus(0.0f), F(1.0f), F_p(1.0f), F_e(1.0f), C(0.0f), velocityGradient(0.0f)  {
        volume = mass / density;
        J = 1.0f;
    };
};
