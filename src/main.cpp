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
#include "particleSplatter.h"
#include "gridRenderer.h"

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

void drawAxes(GLuint shaderProgram, float axisLength = 10.0f)
{
    static GLuint vao = 0, vbo = 0;
    static bool initialized = false;

    if (!initialized)
    {
        glm::vec3 axes[] = {
            // X axis (red)
            glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(axisLength, 0.0f, 0.0f),
            // Y axis (green)
            glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, axisLength, 0.0f),
            // Z axis (blue)
            glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, axisLength)};

        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(axes), axes, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void *)0);
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

int main()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow *window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "MeltSim", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    unsigned int shaderProgram = createShaderProgram("../src/assets/shaders/vertex_shader.glsl", "../src/assets/shaders/fragment_shader.glsl");

    unsigned int gridShaderProgram = createShaderProgramEnhanced(
        "../src/assets/shaders/grid_vertex_shader.glsl",
        "../src/assets/shaders/grid_fragment_shader.glsl");
    
    Plane myPlane(glm::vec3(0, 0, 0), 0, 10.0f);

    // ! Mesh Subsampling
    Mesh myMesh("../src/assets/models/cube.obj");
    std::vector<Vector3> sampledPoints = myMesh.sampleSurfacePoints(0.05f,60,100.0f, 1);

    // std::vector<Vector3> sampledPointsVolume = myMesh.sampleVolumePoints(1000);
    // //append both vectors
    // sampledPoints.insert(sampledPoints.end(), sampledPointsVolume.begin(), sampledPointsVolume.end());
    // std::cout << "Sampled " << sampledPoints.size() << " points from the mesh." << std::endl;

    MPMSimulation mpmSim;
    mpmSim.addMeshParticles(sampledPoints);

     // ! Grid Setup
     GridRenderer gridRenderer(mpmSim.grid.size, mpmSim.grid.spacing);
     gridRenderer.init();

    // ! Particle Setup
    ParticleRenderer particleRenderer;
    particleRenderer.init(mpmSim.particles);

    // ! Particle Splatter Setup
    // Initialize the particle splatter with desired settings
    /*
     * particleRadius: Controls the size of each particle splat
     * - Range: 0.01 - 0.15
     * - Default: 0.05
     * - Effects:
     *   - Smaller values (0.01-0.04): More detailed but potentially grainy fluid
     *   - Medium values (0.05-0.08): Good balance of detail and smoothness
     *   - Larger values (0.09-0.15): Smoother but less detailed fluid appearance
     */

    /*
     * smoothingKernel: Controls the edge softness of particles
     * - Range: 0.5 - 1.0
     * - Default: 0.8
     * - Effects:
     *   - Lower values (0.5-0.6): Sharper particle edges, more distinct particles
     *   - Medium values (0.7-0.8): Natural soft edges
     *   - Higher values (0.9-1.0): Very soft blending between particles
     */
    ParticleSplatter particleSplatter;
    particleSplatter.init(0.15f, 1.00f); // Particle radius and smoothing

    // Set metaball parameters
    /*
     * metaballThreshold: Threshold for when particles start to blend together
     * - Range: 0.5 - 2.0
     * - Default: 1.0
     * - Effects:
     *   - Lower values (0.5-0.8): More unified/blobby fluid appearance
     *   - Medium values (0.9-1.3): Balanced cohesion
     *   - Higher values (1.4-2.0): More distinct particles, less unified
     */

    /*
     * metaballStrength: Intensity of the metaball effect between particles
     * - Range: 0.1 - 1.0
     * - Default: 0.5
     * - Effects:
     *   - Lower values (0.1-0.3): Subtle blending, particles mostly distinct
     *   - Medium values (0.4-0.6): Natural fluid cohesion
     *   - Higher values (0.7-1.0): Strong cohesion, more unified fluid surface
     */

    // particleSplatter.setMetaballThreshold(0.5f); // Threshold for metaball effect
    // particleSplatter.setMetaballStrength(1.0f);  // Strength of the metaball effect

    // // Set water droplet parameters
    // particleSplatter.enableWaterDroplets(false); // Enable water droplet textures
    // particleSplatter.setDropletScale(20.0f);    // Adjust scale of droplet pattern
    // particleSplatter.setDropletIntensity(2.0f);

    // Time step for simulation
    float deltaTime = 0.008f;

    // * Rendering Matrices
    glm::mat4 model = glm::mat4(1.0f);

// View: camera level with cube/cow
glm::vec3 cameraPos = glm::vec3(0.5f, 2.0f, 5.0f);  // higher Y
glm::vec3 target = glm::vec3(0.75f, 0.75f, 0.75f);     // still looking at the cow
glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

    glm::mat4 view = glm::lookAt(cameraPos, target, up);

    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);

    glm::vec3 lightPos(1.2f, 1.0f, 2.0f);
    glm::vec3 lightColor(1.0f, 1.0f, 1.0f);
    glm::vec3 objectColor(0.0f, 0.4f, 0.7f);

    while (!glfwWindowShouldClose(window))
    {
        glDisable(GL_CULL_FACE);
        glEnable(GL_DEPTH_TEST);
        processInput(window);

        // Background color
        // glClearColor(.9f, 0.9f, 0.8f, 1.0f); //off-white
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f); //White

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
        mpmSim.step(deltaTime);
        particleRenderer.update(mpmSim.particles);

        // // ! Util Drawing (Legacy static implementation)
        // drawAxes(shaderProgram, 10.0f);
        // glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.0f, 0.0f, 0.0f);
        // mpmSim.grid.draw();

        // // Draw grid-node debug
        // gridRenderer.update(mpmSim.grid.nodes);
        // gridRenderer.draw(model, view, projection, gridShaderProgram);
        
        // Draw particles with main shader program
        glUseProgram(shaderProgram);

        // ! Draw Plane
        glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.8f, 0.8f, 0.8f);
        myPlane.draw();

        // // Draw Particles with main shader
        // if (objectColorLoc != -1)
        //     glUniform3f(objectColorLoc, 1.0f, 0.0f, 0.0f); // red
        // particleRenderer.draw();

        // In your rendering loop:
        // Update particles with the simulation data
        particleSplatter.update(mpmSim.particles);
        // Draw the particles with the enhanced fluid effects
        particleSplatter.draw(model, view, projection);

        // // ! Draw Mesh
        // glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.2f, 0.5f, 1.0f);
        // myMesh.draw();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}
