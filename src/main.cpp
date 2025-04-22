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

#include "mpmSimulation.h"   

#include <iostream>

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

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
    MPMSimulation mpmSim;

    // Create particles from the cow mesh
    ParticleSampler particleSampler;
    std::vector<Particle> cubeParticles = particleSampler.sampleMeshVolume("../src/assets/models/cube.obj", 0.05f);
    std::cout << "Created " << cubeParticles.size() << " particles from cube mesh" << std::endl;

    // Initialize MPM simulation
    mpmSim.initialize(cubeParticles);


    // Initialize particle renderer
    ParticleRenderer particleRenderer;
    particleRenderer.init("../src/assets/shaders");

    unsigned int floorVAO, floorVBO, floorEBO;
    setupFloor(floorVAO, floorVBO, floorEBO);

    glEnable(GL_DEPTH_TEST);

    glm::mat4 model = glm::mat4(1.0f);

    // View: camera level with cube/cow
    glm::vec3 cameraPos = glm::vec3(0.0f, 1.5f, 5.0f);  // higher Y
    glm::vec3 target = glm::vec3(0.0f, 1.0f, 0.0f);     // still looking at the cow
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

    glm::mat4 view = glm::lookAt(cameraPos, target, up);    
    
    // glm::mat4 view = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, -5.0f));
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);

    glm::vec3 lightPos(1.2f, 1.0f, 2.0f);
    glm::vec3 lightColor(1.0f, 1.0f, 1.0f);
    glm::vec3 objectColor(0.0f, 0.4f, 0.7f);

    while (!glfwWindowShouldClose(window)) {
        // Process input
        processInput(window);
        
        // Update MPM simulation (add this line)
        mpmSim.update(0.016f); // assuming ~60 FPS, or use actual deltaTime
        
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
        
        // // Draw cube at origin - commented out as we're using particles instead
        glm::mat4 model = glm::mat4(1.0f); // Identity matrix
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
        glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.0f, 0.4f, 0.7f);
        // cubeMesh.draw(); //uncomment this line to draw the cube
        
        // // Draw cow on top of the cube (assumes cube height is 1.0)
        glm::mat4 cowModel = glm::translate(model, glm::vec3(0.0f, 0.8f, 0.2f));
        cowModel = glm::rotate(cowModel, glm::radians(-45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(cowModel));
        glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.8f, 0.2f, 0.2f);
        // cowMesh.draw(); //uncomment this line to draw the cow
        
        // Draw particles - Updated to use MPM simulation
        mpmSim.render(view, projection); 
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
}
