#ifndef SHADER_UTILS_H
#define SHADER_UTILS_H

#include <string>
#include <glad/glad.h>

// Function to load shader source from file
std::string loadShaderSource(const std::string& filepath);

// Basic shader program creation (original function)
unsigned int createShaderProgram(const std::string& vertexPath, const std::string& fragmentPath);

// Enhanced shader program creation with detailed error checking
unsigned int createShaderProgramEnhanced(const std::string& vertexPath, const std::string& fragmentPath);

// Helper function to validate shader program
bool validateShaderProgram(unsigned int shaderProgram);

// Helper function to check and print uniform location
GLint getAndCheckUniform(unsigned int shaderProgram, const char* uniformName);

#endif