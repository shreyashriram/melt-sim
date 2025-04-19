#ifndef MESH_LOADER_H
#define MESH_LOADER_H

#include <vector>
#include <string>

void loadMeshData(const std::string& path, std::vector<float>& vertices, std::vector<unsigned int>& indices);

#endif