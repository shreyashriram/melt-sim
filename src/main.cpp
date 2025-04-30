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
#include <fstream>
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

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>


const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;


std::vector<Vector3> loadPoints(const std::string& filename) {
    std::vector<Vector3> points;
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Failed to open file for reading: " << filename << std::endl;
        return points;
    }
    float x, y, z;
    while (in >> x >> y >> z) {
        points.emplace_back(x, y, z);
    }
    return points;
}

void savePoints(const std::vector<Vector3>& points, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }
    for (const auto& p : points) {
        out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    out.close();
}

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
    // * Window Setup
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
    unsigned int gridShaderProgram = createShaderProgramEnhanced("../src/assets/shaders/grid_vertex_shader.glsl","../src/assets/shaders/grid_fragment_shader.glsl");

    // ! Mesh Subsampling
    Mesh myMesh("../src/assets/models/stanford-bunny.obj");

    std::vector<Vector3> sampledPoints;
    std::vector<Vector3> sampledPointsVolume;

    bool loadSamples = true; // you can toggle this manually or check file existence

    if (loadSamples) {
        std::cout << "loading point" << std::endl;
        sampledPoints = loadPoints("bunny_surface_samples.txt");
        sampledPointsVolume = loadPoints("bunny_surface_samples.txt");

        std::cout << "worked" << std::endl;
    } else {
        std::cout << "sampling point" << std::endl;
        sampledPoints = myMesh.sampleSurfacePoints(0.005f, 60, 100.0f, 1);
        sampledPointsVolume = myMesh.sampleVolumePoints(500);
        savePoints(sampledPoints, "bunny_surface_samples.txt");
        savePoints(sampledPointsVolume, "bunny_volume_samples.txt");

        std::cout << "worked" << std::endl;
    }


    MPMSimulation mpmSim;
    mpmSim.addMeshParticles(sampledPoints, MaterialType::Solid);
    mpmSim.addMeshParticles(sampledPointsVolume, MaterialType::Solid);

    // mpmSim.spawnCube(MaterialType::Solid, glm::vec3(1.50f, 2.5f, 1.50f), 0.1f, 10);
    // mpmSim.spawnCube(MaterialType::Solid, glm::vec3(1.50f, 3.0f, 1.50f), 0.1f, 15);

     // ! Grid Setup
    //  GridRenderer gridRenderer(mpmSim.grid.size, mpmSim.grid.spacing);
    //  gridRenderer.init();

    // ! Particle Setup
    // ParticleRenderer particleRenderer;
    // particleRenderer.init(mpmSim.particles);

    // ! Particle Splatter Setup
    ParticleSplatter particleSplatter;
    particleSplatter.init(0.09f, 1.00f); // Particle radius and smoothing

    // ! Time Step
    float deltaTime = 0.04f;
    // float deltaTime = 0.00f;

    // * Rendering Matrices
    glm::mat4 model = glm::mat4(1.0f);

    // View: camera level with cube/cow
    glm::vec3 cameraPos = glm::vec3(7.5f, 2.5f, 1.5f);   // higher Y
    glm::vec3 target = glm::vec3(1.0f, 1.25f, 1.5f);     // still looking at the cow
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

    glm::mat4 view = glm::lookAt(cameraPos, target, up);

    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);

    glm::vec3 lightPos(3.0f, -6.0f, 5.0f);
    glm::vec3 lightColor(1.0f, 1.0f, 1.0f);
    glm::vec3 objectColor(0.0f, 0.4f, 0.7f);

    while (!glfwWindowShouldClose(window))
    {
        glDisable(GL_CULL_FACE);
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE); 
        

        processInput(window);
    
        glClearColor(.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glUseProgram(shaderProgram);

        // Upload camera and lighting uniforms
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
        glUniform3fv(glGetUniformLocation(shaderProgram, "viewPos"), 1, glm::value_ptr(cameraPos));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniform3fv(glGetUniformLocation(shaderProgram, "lightPos"), 1, glm::value_ptr(lightPos));
        glUniform3fv(glGetUniformLocation(shaderProgram, "lightColor"), 1, glm::value_ptr(lightColor));


        mpmSim.step(deltaTime);
        // particleRenderer.update(mpmSim.particles);

        // // ! Util Drawing (Legacy static implementation)
        // drawAxes(shaderProgram, 10.0f);
        glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.0f, 0.0f, 0.0f);
        mpmSim.grid.draw();

        // // Draw grid-node debug
        // gridRenderer.update(mpmSim.grid.nodes);
        // gridRenderer.draw(model, view, projection, gridShaderProgram);
        
        // Draw particles with main shader program
        glUseProgram(shaderProgram);

        // ! Draw Plane
        // glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.8f, 0.8f, 0.8f);
        // myPlane.draw();

        // // Draw Particles with main shader
        // if (objectColorLoc != -1)
        //     glUniform3f(objectColorLoc, 1.0f, 0.0f, 0.0f); // red
        // particleRenderer.draw();

        // In your rendering loop:
        // Update particles with the simulation data
        particleSplatter.update(mpmSim.particles);
        // Draw the particles with the enhanced fluid effects
        particleSplatter.draw(model, view, projection);

        // gridRenderer.update(mpmSim.grid.nodes);
        // gridRenderer.draw(model, view, projection, gridShaderProgram);

        // // ! Draw Mesh
        // glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.2f, 0.5f, 1.0f);
        // myMesh.draw();



        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}