#include "mesh.h"
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <iostream>

#include "leaven/surfaceSampler.h"
#include "leaven/typedef.h"


Mesh::Mesh(const std::string& path) {
    
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenNormals);

    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "ASSIMP ERROR: " << importer.GetErrorString() << std::endl;
        return;
    }

    aiMesh* mesh = scene->mMeshes[0];


    // TODO: split up vertices and normals
    for (unsigned int i = 0; i < mesh->mNumVertices; ++i) {
        vertices.push_back(mesh->mVertices[i].x);
        vertices.push_back(mesh->mVertices[i].y);
        vertices.push_back(mesh->mVertices[i].z);
        vertices.push_back(mesh->mNormals[i].x);
        vertices.push_back(mesh->mNormals[i].y);
        vertices.push_back(mesh->mNormals[i].z);
    }

    for (unsigned int i = 0; i < mesh->mNumFaces; ++i) {
        for (unsigned int j = 0; j < mesh->mFaces[i].mNumIndices; ++j)
            indices.push_back(mesh->mFaces[i].mIndices[j]);
    }
    
    
    setupBuffers();
}

Mesh::~Mesh() {
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
}

void Mesh::setupBuffers() {
   
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    // TODO: split up vertices and normals
    // vertices
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

    // indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);


}

void Mesh::draw() const {
    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
}

void Mesh::bind() const {
    glBindVertexArray(VAO);
}

void Mesh::unBind() const {
    glBindVertexArray(0);
}

std::vector<Vector3> Mesh::sampleSurfacePoints(
    scalar minRadius,
    unsigned int numTrials,
    scalar initialDensity,
    unsigned int distanceNorm
) const {
    using namespace leaven;

    size_t numVertices = vertices.size() / 6;
    Eigen::Matrix<scalar, 3, Eigen::Dynamic> eigenVertices(3, numVertices);

    for (size_t i = 0; i < numVertices; ++i) {
        eigenVertices.col(i) = Vector3(
            vertices[i * 6 + 0],
            vertices[i * 6 + 1],
            vertices[i * 6 + 2]
        );
    }

    size_t numTriangles = indices.size() / 3;
    Eigen::Matrix<unsigned int, 3, Eigen::Dynamic> eigenIndices(3, numTriangles);
    for (size_t i = 0; i < numTriangles; ++i) {
        eigenIndices.col(i) = Eigen::Matrix<unsigned int, 3, 1>(
            indices[i * 3 + 0],
            indices[i * 3 + 1],
            indices[i * 3 + 2]
        );
    }

    SurfaceSampler sampler;
    return sampler.sampleMesh(eigenVertices, eigenIndices, minRadius, numTrials, initialDensity, distanceNorm);
}