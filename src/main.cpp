// File: src/main.cpp
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "shader_utils.h"
#include "mesh.h"
#include "floor.h"
#include "input.h"

// Particle simulation includes
#include "particleSampler.h"
#include "particleRenderer.h"
#include "mpmSimulator.h" // Add the MPM simulator header

#include <iostream>
#include <chrono> // For timing

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

// Timing variables
float deltaTime = 0.0f;
float lastFrame = 0.0f;

int main() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    #ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif

    GLFWwindow *window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "MeltSim", NULL, NULL);
    if (window == NULL) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    unsigned int shaderProgram = createShaderProgram("../src/assets/shaders/vertex_shader.glsl", "../src/assets/shaders/fragment_shader.glsl");

    Mesh cowMesh("../src/assets/models/cow.obj");
    Mesh cubeMesh("../src/assets/models/cube.obj");
    
    Mesh squareMesh("../src/assets/models/square.obj"); //to test 2D

    // Create particles from the cube mesh
    ParticleSampler particleSampler;
    std::vector<Particle> cubeParticles = particleSampler.sampleMeshVolume("../src/assets/models/cube.obj", 0.05f);
    std::cout << "Created " << cubeParticles.size() << " particles from cube mesh" << std::endl;

    // Initialize the MPM simulator
    MPMSimulator mpmSimulator(0.1f, 0.016f); // Grid size 0.1, time step 0.016 (approx 60fps)
    mpmSimulator.initialize(cubeParticles);
    
    // Set simulation parameters
    mpmSimulator.setGravity(glm::vec3(0.0f, -9.81f, 0.0f));
    mpmSimulator.setViscosity(0.05f); // Lower viscosity for water-like behavior
    mpmSimulator.setDensity(1000.0f);
    mpmSimulator.setRestDensity(1000.0f);
    mpmSimulator.setBoundaries(glm::vec3(-2.0f, 0.05f, -2.0f), glm::vec3(2.0f, 4.0f, 2.0f));

    // Initialize particle renderer
    ParticleRenderer particleRenderer;
    particleRenderer.init("../src/assets/shaders");

    unsigned int floorVAO, floorVBO, floorEBO;
    setupFloor(floorVAO, floorVBO, floorEBO);

    glEnable(GL_DEPTH_TEST);

    glm::mat4 model = glm::mat4(1.0f);

    // View: camera level with cube/cow
    glm::vec3 cameraPos = glm::vec3(0.0f, 1.5f, 5.0f);  // higher Y
    //Rotate camera 45 degrees around Y axis
    cameraPos = glm::vec3(glm::rotate(glm::mat4(1.0f), glm::radians(-45.0f), glm::vec3(0.0f, 1.0f, 0.0f)) * glm::vec4(cameraPos, 1.0f));
    glm::vec3 target = glm::vec3(0.0f, 1.0f, 0.0f);     // still looking at the cow
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

    glm::mat4 view = glm::lookAt(cameraPos, target, up);    
    
    // glm::mat4 view = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, -5.0f));
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);

    glm::vec3 lightPos(1.2f, 1.0f, 2.0f);
    glm::vec3 lightColor(1.0f, 1.0f, 1.0f);
    glm::vec3 objectColor(0.0f, 0.4f, 0.7f);

    // Simulation flag - pause/play with spacebar
    bool simulationRunning = true;
    
    // Simulation rate - can be slower than rendering if needed
    const float simulationTimeStep = 0.016f; // 60 Hz
    float accumulatedTime = 0.0f;

    while (!glfwWindowShouldClose(window)) {
        // Calculate delta time
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        
        // Process input
        processInput(window);
        
        // Toggle simulation with spacebar
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) {
            static float lastToggleTime = 0.0f;
            if (currentFrame - lastToggleTime > 0.3f) { // Debounce
                simulationRunning = !simulationRunning;
                lastToggleTime = currentFrame;
                std::cout << "Simulation " << (simulationRunning ? "running" : "paused") << std::endl;
            }
        }
        
        // Update simulation at fixed time step
        if (simulationRunning) {
            accumulatedTime += deltaTime;
            
            // Run simulation step(s)
            while (accumulatedTime >= simulationTimeStep) {
                mpmSimulator.update(cubeParticles);
                accumulatedTime -= simulationTimeStep;
            }
        }
        
        // Rendering
        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        glUseProgram(shaderProgram);
        // Upload camera and lighting uniforms
        glUniform3fv(glGetUniformLocation(shaderProgram, "viewPos"), 1, glm::value_ptr(cameraPos));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniform3fv(glGetUniformLocation(shaderProgram, "lightPos"), 1, glm::value_ptr(lightPos));
        glUniform3fv(glGetUniformLocation(shaderProgram, "lightColor"), 1, glm::value_ptr(lightColor));
        
        // Draw floor at origin
        glm::mat4 floorModel = glm::mat4(1.0f);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(floorModel));
        glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.8f, 0.8f, 0.8f);
        glBindVertexArray(floorVAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        
        // Draw cube at origin - commented out as we're using particles instead
        // glm::mat4 model = glm::mat4(1.0f); // Identity matrix
        // glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
        // glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.0f, 0.4f, 0.7f);
        // cubeMesh.draw(); //uncomment this line to draw the cube
        
        // Draw cow on top of the cube (assumes cube height is 1.0)
        // glm::mat4 cowModel = glm::translate(model, glm::vec3(0.0f, 0.8f, 0.2f));
        // cowModel = glm::rotate(cowModel, glm::radians(-45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        // glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(cowModel));
        // glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.8f, 0.2f, 0.2f);
        // cowMesh.draw(); //uncomment this line to draw the cow
        
        // Use render function to draw particles - with updated positions from MPM simulation
        particleRenderer.render(cubeParticles, view, projection);

        // Display simulation status
        // For a real UI, you'd want to use ImGui or similar
        std::string windowTitle = "MeltSim - MPM Fluid [";
        windowTitle += simulationRunning ? "Running" : "Paused";
        windowTitle += "] - FPS: " + std::to_string(1.0f / deltaTime);
        glfwSetWindowTitle(window, windowTitle.c_str());
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    // Clean up
    glDeleteVertexArrays(1, &floorVAO);
    glDeleteBuffers(1, &floorVBO);
    glDeleteBuffers(1, &floorEBO);
    glDeleteProgram(shaderProgram);
    
    glfwTerminate();
    return 0;
}