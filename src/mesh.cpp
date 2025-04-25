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

// Volume sampling using ray casting
// This function samples points inside the volume of the mesh using a ray-casting method.
// It generates random points within the bounding box of the mesh and checks if they are inside the mesh
// by casting a ray in a fixed direction (e.g., +x direction) and counting the number of intersections with the mesh triangles.
// If the number of intersections is odd, the point is inside the mesh; if even, it's outside.
// The function returns a vector of sampled points inside the mesh volume.
std::vector<Vector3> Mesh::sampleVolumePoints(unsigned int numPoints) const {
    std::vector<Vector3> volumePoints;

    // Step 1: Compute bounding box of the mesh
    Vector3 minBB = Vector3::Constant(std::numeric_limits<float>::max());
    Vector3 maxBB = Vector3::Constant(std::numeric_limits<float>::lowest());

    for (size_t i = 0; i < vertices.size(); i += 6) {
        Vector3 v(vertices[i], vertices[i + 1], vertices[i + 2]);
        minBB = minBB.cwiseMin(v);
        maxBB = maxBB.cwiseMax(v);
    }

    // Step 2: Flatten triangle vertices
    std::vector<Vector3> triangleVertices;
    for (size_t i = 0; i < indices.size(); ++i) {
        size_t idx = indices[i] * 6;
        triangleVertices.emplace_back(vertices[idx], vertices[idx + 1], vertices[idx + 2]);
    }

    // Step 3: Sample random points and use ray casting
    std::default_random_engine rng(std::random_device{}());
    std::uniform_real_distribution<float> distX(minBB.x(), maxBB.x());
    std::uniform_real_distribution<float> distY(minBB.y(), maxBB.y());
    std::uniform_real_distribution<float> distZ(minBB.z(), maxBB.z());

    while (volumePoints.size() < numPoints) {
        Vector3 p(distX(rng), distY(rng), distZ(rng));
        Vector3 dir(1.0f, 0.0f, 0.0f); // Ray in +x direction

        int intersections = 0;
        for (size_t i = 0; i < triangleVertices.size(); i += 3) {
            if (rayIntersectsTriangle(p, dir,
                triangleVertices[i],
                triangleVertices[i + 1],
                triangleVertices[i + 2])) {
                ++intersections;
            }
        }

        if (intersections % 2 == 1) {
            volumePoints.push_back(p); // Inside mesh
        }
    }

    return volumePoints;
}

//Uses the Möller–Trumbore algorithm for ray-triangle intersection
// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_algorithm
// This function checks if a ray intersects a triangle in 3D space
// orig: ray origin
// dir: ray direction
// v0, v1, v2: vertices of the triangle
// Returns true if the ray intersects the triangle, false otherwise
bool Mesh::rayIntersectsTriangle(
    const Vector3& orig,
    const Vector3& dir,
    const Vector3& v0,
    const Vector3& v1,
    const Vector3& v2
) const {
    const float EPSILON = 1e-6f;
    Vector3 edge1 = v1 - v0;
    Vector3 edge2 = v2 - v0;
    Vector3 h = dir.cross(edge2);
    float a = edge1.dot(h);
    if (std::fabs(a) < EPSILON)
        return false;

    float f = 1.0f / a;
    Vector3 s = orig - v0;
    float u = f * s.dot(h);
    if (u < 0.0f || u > 1.0f)
        return false;

    Vector3 q = s.cross(edge1);
    float v = f * dir.dot(q);
    if (v < 0.0f || u + v > 1.0f)
        return false;

    float t = f * edge2.dot(q);
    return t > EPSILON;
}