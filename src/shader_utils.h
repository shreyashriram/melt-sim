#ifndef SHADER_UTILS_H
#define SHADER_UTILS_H

#include <string>

unsigned int createShaderProgram(const std::string& vertexPath, const std::string& fragmentPath);
std::string loadShaderSource(const std::string& filepath);

#endif