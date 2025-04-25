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
#include <fstream>
#include <sstream>
#include <string>

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
#include "gridRenderer.h"

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
    
    // Linking shaders with enhanced error checking
    std::cout << "Compiling and linking grid shader program..." << std::endl;
    unsigned int gridShaderProgram = createShaderProgramEnhanced(
        "../src/assets/shaders/grid_vertex_shader.glsl", 
        "../src/assets/shaders/grid_fragment_shader.glsl"
    );
    // if (gridShaderProgram == 0) {
    //     std::cout << "Failed to compile/link grid shader program!" << std::endl;
    //     glfwTerminate();
    //     return -1;
    // }
    
    std::cout << "Compiling and linking main shader program..." << std::endl;
    unsigned int shaderProgram = createShaderProgramEnhanced(
        "../src/assets/shaders/vertex_shader.glsl", 
        "../src/assets/shaders/fragment_shader.glsl"
    );
    // if (shaderProgram == 0) {
    //     std::cout << "Failed to compile/link main shader program!" << std::endl;
    //     glfwTerminate();
    //     return -1;
    // }
    
    // // Validate shader programs
    // std::cout << "Validating grid shader program..." << std::endl;
    // if (!validateShaderProgram(gridShaderProgram)) {
    //     std::cout << "Grid shader program validation failed!" << std::endl;
    // }
    
    // std::cout << "Validating main shader program..." << std::endl;
    // if (!validateShaderProgram(shaderProgram)) {
    //     std::cout << "Main shader program validation failed!" << std::endl;
    // }
    Plane myPlane(glm::vec3(0, 0, 0), 0, 10.0f);
    
    // ! Mesh Subsampling 
    Mesh myMesh("../src/assets/models/cube.obj");
    // std::vector<Vector3> sampledPoints = myMesh.sampleSurfacePoints(0.05f,60,100.0f, 1);
    std::vector<Vector3> sampledPoints = myMesh.sampleVolumePoints(1000);
    std::cout << "Sampled " << sampledPoints.size() << " points from the mesh." << std::endl;
    
    MPMSimulation mpmSim;
    mpmSim.addMeshParticles(sampledPoints);
    // ! Particle Setup

    ParticleRenderer particleRenderer;
    particleRenderer.init(mpmSim.particles);

    // ! Grid Setup
    GridRenderer gridRenderer(mpmSim.getGridSize(), mpmSim.getGridSpacing());
    gridRenderer.init();
    
    float deltaTime = 0.01f;

    // * Rendering Matrices
    glm::mat4 model = glm::mat4(1.0f);

    // View: camera level with cube/cow
    glm::vec3 cameraPos = glm::vec3(0.0f, 0.75f, 2.5f);  // higher Y
    glm::vec3 target = glm::vec3(0.0f, 0.5f, 0.0f);     // still looking at the cow
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

    glm::mat4 view = glm::lookAt(cameraPos, target, up);    

    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);

    glm::vec3 lightPos(1.2f, 1.0f, 2.0f);
    glm::vec3 lightColor(1.0f, 1.0f, 1.0f);
    glm::vec3 objectColor(0.0f, 0.4f, 0.7f);
    
    // Print all the uniform locations for debugging
    std::cout << "Checking uniform locations in main shader program..." << std::endl;
    glUseProgram(shaderProgram);
    GLint modelLoc = getAndCheckUniform(shaderProgram, "model");
    GLint viewPosLoc = getAndCheckUniform(shaderProgram, "viewPos");
    GLint viewLoc = getAndCheckUniform(shaderProgram, "view");
    GLint projectionLoc = getAndCheckUniform(shaderProgram, "projection");
    GLint lightPosLoc = getAndCheckUniform(shaderProgram, "lightPos");
    GLint lightColorLoc = getAndCheckUniform(shaderProgram, "lightColor");
    GLint objectColorLoc = getAndCheckUniform(shaderProgram, "objectColor");
    
    std::cout << "Checking uniform locations in grid shader program..." << std::endl;
    glUseProgram(gridShaderProgram);
    GLint gridModelLoc = getAndCheckUniform(gridShaderProgram, "model");
    GLint gridViewLoc = getAndCheckUniform(gridShaderProgram, "view");
    GLint gridProjectionLoc = getAndCheckUniform(gridShaderProgram, "projection");
    
    while (!glfwWindowShouldClose(window)) {
        glDisable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        processInput(window);
    
        glClearColor(0.6f, 0.8f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
        // Use and validate main shader program
        glUseProgram(shaderProgram);
        
        // Upload camera and lighting uniforms with error checking
        if (modelLoc != -1)
            glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
        if (viewPosLoc != -1)
            glUniform3fv(viewPosLoc, 1, glm::value_ptr(cameraPos));
        if (viewLoc != -1)
            glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
        if (projectionLoc != -1)
            glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
        if (lightPosLoc != -1)
            glUniform3fv(lightPosLoc, 1, glm::value_ptr(lightPos));
        if (lightColorLoc != -1)
            glUniform3fv(lightColorLoc, 1, glm::value_ptr(lightColor));
    
        // Single-step pipeline
        mpmSim.step(deltaTime);

        // Use and validate grid shader program
        glUseProgram(gridShaderProgram);
        
        // Upload matrices for grid shader
        if (gridModelLoc != -1)
            glUniformMatrix4fv(gridModelLoc, 1, GL_FALSE, glm::value_ptr(model));
        if (gridViewLoc != -1)
            glUniformMatrix4fv(gridViewLoc, 1, GL_FALSE, glm::value_ptr(view));
        if (gridProjectionLoc != -1)
            glUniformMatrix4fv(gridProjectionLoc, 1, GL_FALSE, glm::value_ptr(projection));

        // Draw grid-node debug
        gridRenderer.update(mpmSim.getGridNodes()); 
        gridRenderer.draw(model, view, projection, gridShaderProgram);  

        // Draw particles with main shader program
        glUseProgram(shaderProgram);
        particleRenderer.update(mpmSim.particles);        

        // Draw Plane with main shader
        if (objectColorLoc != -1)
            glUniform3f(objectColorLoc, 0.8f, 0.8f, 0.8f);
        myPlane.draw();

        // Draw Particles with main shader
        if (objectColorLoc != -1)
            glUniform3f(objectColorLoc, 1.0f, 0.0f, 0.0f); // red
        particleRenderer.draw();

        // Draw Mesh with main shader (commented out in original code)
        // glm::mat4 meshModel = glm::translate(model, glm::vec3(0.0f, 1.0f, 0.0f));
        // // meshModel = glm::rotate(meshModel, glm::radians(-45.0f), glm::vec3(0.0f, 1.0f, 0.0f));        
        // if (modelLoc != -1)
        //     glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(meshModel));
        // if (objectColorLoc != -1)
        //     glUniform3f(objectColorLoc, 0.2f, 0.5f, 1.0f);
        // myMesh.draw();
    
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
    // Clean up
    glDeleteProgram(shaderProgram);
    glDeleteProgram(gridShaderProgram);
    glfwTerminate();
    return 0;
}