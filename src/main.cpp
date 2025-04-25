// File: src/main.cpp
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "shader_utils.h"
#include "mesh.h"
#include "floor.h"
#include "plane.h"
#include "input.h"

#include <iostream>

#include "leaven/surfaceSampler.h"  
#include "leaven/typedef.h"
using namespace leaven;
#include <Eigen/Dense>
using scalar = float;
using Vector3 = Eigen::Matrix<scalar, 3, 1>;

#include "particle.h"
#include "particleRenderer.h"
#include "grid.h"
#include "mpm.h"

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

void drawAxes(GLuint shaderProgram, float axisLength = 10.0f) {
    static GLuint vao = 0, vbo = 0;
    static bool initialized = false;

    if (!initialized) {
        glm::vec3 axes[] = {
            // X axis (red)
            glm::vec3(-(axisLength), 0.0f, 0.0f), glm::vec3(axisLength, 0.0f, 0.0f),
            // Y axis (green)
            glm::vec3(0.0f, -(axisLength), 0.0f), glm::vec3(0.0f, axisLength, 0.0f),
            // Z axis (blue)
            glm::vec3(0.0f, 0.0f, -(axisLength)), glm::vec3(0.0f, 0.0f, axisLength)
        };

        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(axes), axes, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
        glBindVertexArray(0);

        initialized = true;
    }

    glUseProgram(shaderProgram);

    // Model matrix = identity (origin-centered)
    glm::mat4 model = glm::mat4(1.0f);
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));

    glBindVertexArray(vao);

    // X axis (red)
    glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 1.0f, 0.0f, 0.0f);
    glDrawArrays(GL_LINES, 0, 2);

    // Y axis (green)
    glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.0f, 1.0f, 0.0f);
    glDrawArrays(GL_LINES, 2, 2);

    // Z axis (blue)
    glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.0f, 0.0f, 1.0f);
    glDrawArrays(GL_LINES, 4, 2);

    glBindVertexArray(0);
}


void drawGridLines(GLuint shaderProgram, int gridSize = 5, float spacing = 0.25f) {
    static GLuint vao = 0, vbo = 0;
    static bool initialized = false;

    if (!initialized) {
        std::vector<glm::vec3> lines;

        float start = 0.0f; // bottom-left-front at (0,0,0)
        float end = gridSize * spacing;

        // Generate lines
        for (int i = 0; i < gridSize + 1; ++i) {
            float offset = i * spacing + start;

            // Y-Z planes (lines along X)
            for (int j = 0; j < gridSize + 1; ++j) {
                float y = j * spacing + start;
                lines.emplace_back(start, y, start - i * spacing); // notice start - offset (for Z)
                lines.emplace_back(end,   y, start - i * spacing);
            }

            // X-Z planes (lines along Y)
            for (int j = 0; j < gridSize + 1; ++j) {
                float x = j * spacing + start;
                lines.emplace_back(x, start, start - i * spacing);
                lines.emplace_back(x, end,   start - i * spacing);
            }

            // X-Y planes (lines along Z)
            for (int j = 0; j < gridSize + 1; ++j) {
                float x = j * spacing + start;
                lines.emplace_back(x, offset, start);
                lines.emplace_back(x, offset, start - end);
            }
        }

        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, lines.size() * sizeof(glm::vec3), lines.data(), GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
        glBindVertexArray(0);

        initialized = true;
    }

    glUseProgram(shaderProgram);

    glm::mat4 model = glm::mat4(1.0f);
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
    glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.6f, 0.6f, 0.6f); // gray lines

    glBindVertexArray(vao);
    glDrawArrays(GL_LINES, 0, (gridSize + 1) * (gridSize + 1) * 3 * 2);
    glBindVertexArray(0);
}



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


    Plane myPlane(glm::vec3(0, 0, 0), 0, 10.0f);
    
    // ! Mesh Subsampling 
    Mesh myMesh("../src/assets/models/cow.obj");
    std::vector<Vector3> sampledPoints = myMesh.sampleSurfacePoints(0.05f,60,100.0f, 1);

    
    MPMSimulation mpmSim;
    mpmSim.addMeshParticles(sampledPoints);
    
    // ! Particle Setup
    ParticleRenderer particleRenderer;
    particleRenderer.init(mpmSim.particles);

    
    float deltaTime = 0.0045f;

    // * Rendering Matrices
    glm::mat4 model = glm::mat4(1.0f);

    // View: camera level with cube/cow
    glm::vec3 cameraPos = glm::vec3(-0.5f, 2.0f, 2.5f);  // higher Y
    glm::vec3 target = glm::vec3(0.5f, 0.5f, -0.5f);     // still looking at the cow
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

    glm::mat4 view = glm::lookAt(cameraPos, target, up);    

    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);

    glm::vec3 lightPos(1.2f, 1.0f, 2.0f);
    glm::vec3 lightColor(1.0f, 1.0f, 1.0f);
    glm::vec3 objectColor(0.0f, 0.4f, 0.7f);
    

    while (!glfwWindowShouldClose(window)) {
        glDisable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        processInput(window);
    
        glClearColor(0.6f, 0.8f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
        glUseProgram(shaderProgram);
    
        // Upload camera and lighting uniforms
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
        glUniform3fv(glGetUniformLocation(shaderProgram, "viewPos"), 1, glm::value_ptr(cameraPos));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniform3fv(glGetUniformLocation(shaderProgram, "lightPos"), 1, glm::value_ptr(lightPos));
        glUniform3fv(glGetUniformLocation(shaderProgram, "lightColor"), 1, glm::value_ptr(lightColor));
    
        // Update simulation
        // mpmSim.step(deltaTime);
        particleRenderer.update(mpmSim.particles);
        drawGridLines(shaderProgram);
        drawAxes(shaderProgram, 10.0f);


        // ! Draw Plane
        // glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.8f, 0.8f, 0.8f);
        // myPlane.draw();

        // ! Draw Particles
        glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 1.0f, 0.0f, 0.0f); // red
        particleRenderer.draw();

        // ! Draw Mesh
        // glm::mat4 meshModel = glm::translate(model, glm::vec3(0.0f, 1.0f, 0.0f));
        // // // meshModel = glm::rotate(meshModel, glm::radians(-45.0f), glm::vec3(0.0f, 1.0f, 0.0f));        
        // glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(meshModel));
        // glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.2f, 0.5f, 1.0f);
        // myMesh.draw();
    
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
}
