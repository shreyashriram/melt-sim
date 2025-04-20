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

#include <iostream>

#include "leaven/surfaceSampler.h"   // ✅ This is key
#include "leaven/typedef.h"
using namespace leaven;
#include <Eigen/Dense>
using scalar = float;
using Vector3 = Eigen::Matrix<scalar, 3, 1>;


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

    unsigned int floorVAO, floorVBO, floorEBO;
    setupFloor(floorVAO, floorVBO, floorEBO);

    glEnable(GL_DEPTH_TEST);

    // ! Mesh Subsampling 
    Mesh cowMesh("../src/assets/models/cow.obj");
    // Mesh cowMesh("../src/assets/models/cube.obj");

    size_t numVertices = cowMesh.vertices.size() / 6;
    Eigen::Matrix<scalar, 3, Eigen::Dynamic> eigenVertices(3, numVertices);

    for (size_t i = 0; i < numVertices; ++i) {
        float x = cowMesh.vertices[i * 6 + 0];
        float y = cowMesh.vertices[i * 6 + 1];
        float z = cowMesh.vertices[i * 6 + 2];
        eigenVertices.col(i) = Vector3(x, y, z);
    }
    
    // ? Debugging 
    // std::cout << "Loaded " << numVertices << " vertices.\n";
    // for (size_t i = 0; i < std::min(numVertices, size_t(5)); ++i) {
    //     std::cout << "Vertex " << i << ": " << eigenVertices.col(i).transpose() << "\n";
    // }


    size_t numTriangles = cowMesh.indices.size() / 3;
    Eigen::Matrix<unsigned int, 3, Eigen::Dynamic> eigenIndices(3, numTriangles);

    for (size_t i = 0; i < numTriangles; ++i) {
        eigenIndices.col(i) = Eigen::Matrix<unsigned int, 3, 1>(
         cowMesh.indices[i * 3 + 0],
         cowMesh.indices[i * 3 + 1],
         cowMesh.indices[i * 3 + 2]
        );
    }
    
    // ? Debugging    
    // std::cout << "Loaded " << numTriangles << " triangles.\n";
    // for (size_t i = 0; i < std::min(numTriangles, size_t(5)); ++i) {
    //     std::cout << "Triangle " << i << ": "
    //             << eigenIndices(0, i) << ", "
    //             << eigenIndices(1, i) << ", "
    //             << eigenIndices(2, i) << "\n";
    // }

    scalar minRadius = 0.01f;             // spacing between samples (smaller = more points)
    unsigned int numTrials = 60;         // number of sampling attempts
    scalar initialDensity = 5.0f;        // how many initial points to generate
    unsigned int distanceNorm = 2;       // 2 = Euclidean distance (common for 3D)

    SurfaceSampler sampler;

    std::vector<Vector3> sampledPoints = sampler.sampleMesh(
        eigenVertices,
        eigenIndices,
        minRadius,
        numTrials,
        initialDensity,
        distanceNorm
    );

    // ? Debugging
    // std::cout << "Sampled " << sampledPoints.size() << " points from the cube mesh.\n";
    // for (int i = 0; i < std::min(10, (int)sampledPoints.size()); ++i) {
    //     std::cout << "Point " << i << ": " << sampledPoints[i].transpose() << "\n";
    // }


    std::vector<float> pointData;
    for (const auto& p : sampledPoints) {
        pointData.push_back(p.x());
        pointData.push_back(p.y());
        pointData.push_back(p.z());
    }

    unsigned int pointsVAO, pointsVBO;
    glGenVertexArrays(1, &pointsVAO);
    glGenBuffers(1, &pointsVBO);

    glBindVertexArray(pointsVAO);
    glBindBuffer(GL_ARRAY_BUFFER, pointsVBO);
    glBufferData(GL_ARRAY_BUFFER, pointData.size() * sizeof(float), pointData.data(), GL_STATIC_DRAW);

    // layout (location = 0) = vec3 position
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(0);



    glm::mat4 model = glm::mat4(1.0f);

    // View: camera level with cube/cow
    glm::vec3 cameraPos = glm::vec3(0.0f, 0.5f, 3.0f);  // higher Y
    glm::vec3 target = glm::vec3(0.0f, 0.0f, 0.0f);     // still looking at the cow
    glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);

    glm::mat4 view = glm::lookAt(cameraPos, target, up);    


    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);

    glm::vec3 lightPos(1.2f, 1.0f, 2.0f);
    glm::vec3 lightColor(1.0f, 1.0f, 1.0f);
    glm::vec3 objectColor(0.0f, 0.4f, 0.7f);

    while (!glfwWindowShouldClose(window)) {
        processInput(window);
    
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
    
        // Draw cube at origin
        // glm::mat4 model = glm::mat4(1.0f); // Identity matrix
        // glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
        // glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.0f, 0.4f, 0.7f);
        // cowMesh.draw();
    
        glBindVertexArray(pointsVAO);
        glPointSize(7.0f);  // increase for visibility
        glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 1.0f, 0.2f, 0.2f); // red glow
        glDrawArrays(GL_POINTS, 0, sampledPoints.size());


        // Draw cow on top of the cube (assumes cube height is 1.0)
        // glm::mat4 cowModel = glm::translate(model, glm::vec3(0.0f, 0.8f, 0.2f));
        // cowModel = glm::rotate(cowModel, glm::radians(-45.0f), glm::vec3(0.0f, 1.0f, 0.0f));        
        // glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(cowModel));
        // glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.8f, 0.2f, 0.2f);
        // cowMesh.draw();
    
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    
}
