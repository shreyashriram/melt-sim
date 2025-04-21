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
#include "MPMSimulator.h"
#include "particle_renderer.h"

#include <iostream>

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    // Set OpenGL version to 3.3 core
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // Create GLFW window
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "MeltSim", nullptr, nullptr);
    if (!window) {
        std::cerr << "GLFW Error: Failed to create window." << std::endl;
        std::cerr << "Possible causes:\n";
        std::cerr << "- Unsupported OpenGL version (3.3 core requested)\n";
        std::cerr << "- Graphics driver issue\n";
        std::cerr << "- Running in headless environment without display\n";
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Load OpenGL function pointers using GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);  // To allow gl_PointSize in shader

    // Load shaders
    unsigned int floorShader = createShaderProgram("../src/assets/shaders/vertex_shader.glsl", "../src/assets/shaders/fragment_shader.glsl");
    unsigned int particleShader = createShaderProgram("../src/assets/shaders/particle.vert", "../src/assets/shaders/particle.frag");

    // Set up objects
    unsigned int floorVAO, floorVBO, floorEBO;
    setupFloor(floorVAO, floorVBO, floorEBO);
    Mesh cubeMesh("../src/assets/models/cube.obj");

    // Camera setup
    glm::vec3 cameraPos = glm::vec3(0.0f, 1.5f, 5.0f);
    glm::vec3 target = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::mat4 view = glm::lookAt(cameraPos, target, up);
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / SCR_HEIGHT, 0.1f, 100.0f);

    // MPM simulation
    MPMSimulator simulator;
    ParticleRenderer particleRenderer;

    // Render loop
    while (!glfwWindowShouldClose(window)) {
        processInput(window);
        simulator.step();
        auto& particles = simulator.getParticles();

        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Draw floor
        glUseProgram(floorShader);
        glUniformMatrix4fv(glGetUniformLocation(floorShader, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(floorShader, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniformMatrix4fv(glGetUniformLocation(floorShader, "model"), 1, GL_FALSE, glm::value_ptr(glm::mat4(1.0f)));
        glUniform3f(glGetUniformLocation(floorShader, "objectColor"), 0.8f, 0.8f, 0.8f);
        glUniform3fv(glGetUniformLocation(floorShader, "viewPos"), 1, glm::value_ptr(cameraPos));
        glUniform3f(glGetUniformLocation(floorShader, "lightPos"), 1.2f, 1.0f, 2.0f);
        glUniform3f(glGetUniformLocation(floorShader, "lightColor"), 1.0f, 1.0f, 1.0f);
        glBindVertexArray(floorVAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // Draw particles
        glUseProgram(particleShader);
        glUniformMatrix4fv(glGetUniformLocation(particleShader, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(particleShader, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        particleRenderer.update(particles);
        particleRenderer.draw();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}