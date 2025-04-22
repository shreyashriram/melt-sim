// particleSampler.cpp
#include "particleSampler.h"
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <glad/glad.h>
#include <iostream>
#include <random>
#include <cmath>

ParticleSampler::ParticleSampler() {}

ParticleSampler::~ParticleSampler() {}

std::vector<Particle> ParticleSampler::sampleMeshVolume(const std::string& meshPath, float particleRadius) {
    std::vector<Particle> particles;
    
    std::vector<glm::vec3> vertices = extractVertices(meshPath);
    std::vector<unsigned int> indices = extractIndices(meshPath);
    
    if (vertices.empty() || indices.empty()) {
        std::cerr << "Failed to load mesh for particle sampling" << std::endl;
        return particles;
    }
    
    // Calculate bounding box
    glm::vec3 bbMin = calculateBoundingBoxMin(vertices);
    glm::vec3 bbMax = calculateBoundingBoxMax(vertices);
    
    // Add a small margin to the bounding box
    bbMin -= glm::vec3(particleRadius);
    bbMax += glm::vec3(particleRadius);
    
    // Calculate grid dimensions based on particle radius
    glm::vec3 dimensions = bbMax - bbMin;
    int numX = std::ceil(dimensions.x / (2.0f * particleRadius));
    int numY = std::ceil(dimensions.y / (2.0f * particleRadius));
    int numZ = std::ceil(dimensions.z / (2.0f * particleRadius));
    
    // Grid spacing
    float spacing = 2.0f * particleRadius;
    
    // Add some jitter to prevent perfect grid alignment
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> jitterDist(-0.1f * particleRadius, 0.1f * particleRadius);
    
    // Create particles inside the mesh
    for (int x = 0; x < numX; x++) {
        for (int y = 0; y < numY; y++) {
            for (int z = 0; z < numZ; z++) {
                glm::vec3 position = bbMin + glm::vec3(
                    x * spacing + jitterDist(gen), 
                    y * spacing + jitterDist(gen), 
                    z * spacing + jitterDist(gen)
                );
                
                // Check if the position is inside the mesh
                if (isPointInMesh(position, vertices, indices)) {
                    particles.emplace_back(position);
                }
            }
        }
    }
    
    std::cout << "Created " << particles.size() << " volume particles from mesh" << std::endl;
    return particles;
}

std::vector<Particle> ParticleSampler::sampleMeshSurface(const std::string& meshPath, float particleSpacing) {
    std::vector<Particle> particles;
    
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(meshPath, 
        aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenNormals);
    
    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "ASSIMP ERROR: " << importer.GetErrorString() << std::endl;
        return particles;
    }
    
    aiMesh* mesh = scene->mMeshes[0];
    
    // Sample points on each triangle face
    for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
        const aiFace& face = mesh->mFaces[i];
        if (face.mNumIndices != 3) continue; // Skip non-triangular faces
        
        // Get triangle vertices
        aiVector3D v0 = mesh->mVertices[face.mIndices[0]];
        aiVector3D v1 = mesh->mVertices[face.mIndices[1]];
        aiVector3D v2 = mesh->mVertices[face.mIndices[2]];
        
        // Get triangle normals
        aiVector3D n0 = mesh->mNormals[face.mIndices[0]];
        aiVector3D n1 = mesh->mNormals[face.mIndices[1]];
        aiVector3D n2 = mesh->mNormals[face.mIndices[2]];
        
        // Calculate triangle area
        aiVector3D edge1 = v1 - v0;
        aiVector3D edge2 = v2 - v0;
        aiVector3D cross = edge1 ^ edge2;
        float area = 0.5f * sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
        
        // Calculate number of particles to place on this triangle
        int numParticles = std::max(1, static_cast<int>(area / (particleSpacing * particleSpacing)));
        
        // Generate random barycentric coordinates for particles
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);
        
        for (int j = 0; j < numParticles; j++) {
            // Generate barycentric coordinates
            float u = dist(gen);
            float v = dist(gen);
            if (u + v > 1.0f) {
                u = 1.0f - u;
                v = 1.0f - v;
            }
            float w = 1.0f - u - v;
            
            // Calculate position using barycentric coordinates
            aiVector3D position = v0 * u + v1 * v + v2 * w;
            aiVector3D normal = n0 * u + n1 * v + n2 * w;
            normal.Normalize();
            
            // Create particle
            Particle particle(glm::vec3(position.x, position.y, position.z));
            particles.push_back(particle);
        }
    }
    
    std::cout << "Created " << particles.size() << " surface particles from mesh" << std::endl;
    return particles;
}

bool ParticleSampler::isPointInMesh(const glm::vec3& point, 
                                    const std::vector<glm::vec3>& vertices, 
                                    const std::vector<unsigned int>& indices) {
    // Ray casting algorithm for point-in-mesh test
    // Cast a ray from the point in any direction and count intersections
    // Odd number of intersections = inside, even = outside
    
    glm::vec3 rayDir(1.0f, 0.0f, 0.0f); // Cast ray in positive X direction
    int intersections = 0;
    
    // Check each triangle for intersection
    for (size_t i = 0; i < indices.size(); i += 3) {
        glm::vec3 v0 = vertices[indices[i]];
        glm::vec3 v1 = vertices[indices[i+1]];
        glm::vec3 v2 = vertices[indices[i+2]];
        
        // Calculate triangle normal
        glm::vec3 edge1 = v1 - v0;
        glm::vec3 edge2 = v2 - v0;
        glm::vec3 normal = glm::normalize(glm::cross(edge1, edge2));
        
        // Ray-plane intersection
        float denom = glm::dot(normal, rayDir);
        if (std::abs(denom) > 0.0001f) { // Non-parallel ray and plane
            float t = glm::dot(v0 - point, normal) / denom;
            if (t > 0.0f) { // Intersection is in positive ray direction
                // Calculate intersection point
                glm::vec3 intersection = point + t * rayDir;
                
                // Check if intersection point is inside triangle
                // Using barycentric coordinates
                glm::vec3 edge3 = intersection - v0;
                float dot00 = glm::dot(edge1, edge1);
                float dot01 = glm::dot(edge1, edge2);
                float dot02 = glm::dot(edge1, edge3);
                float dot11 = glm::dot(edge2, edge2);
                float dot12 = glm::dot(edge2, edge3);
                
                float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
                float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
                float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
                
                if (u >= 0.0f && v >= 0.0f && u + v <= 1.0f) {
                    intersections++;
                }
            }
        }
    }
    
    // Odd number of intersections means the point is inside the mesh
    return (intersections % 2) == 1;
}

std::vector<glm::vec3> ParticleSampler::extractVertices(const std::string& meshPath) {
    std::vector<glm::vec3> vertices;
    
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(meshPath, 
        aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenNormals);
    
    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "ASSIMP ERROR: " << importer.GetErrorString() << std::endl;
        return vertices;
    }
    
    aiMesh* mesh = scene->mMeshes[0];
    
    vertices.reserve(mesh->mNumVertices);
    for (unsigned int i = 0; i < mesh->mNumVertices; i++) {
        aiVector3D vertex = mesh->mVertices[i];
        vertices.emplace_back(vertex.x, vertex.y, vertex.z);
    }
    
    return vertices;
}

std::vector<unsigned int> ParticleSampler::extractIndices(const std::string& meshPath) {
    std::vector<unsigned int> indices;
    
    Assimp::Importer importer;
    const aiScene* scene = importer.ReadFile(meshPath, 
        aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_GenNormals);
    
    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "ASSIMP ERROR: " << importer.GetErrorString() << std::endl;
        return indices;
    }
    
    aiMesh* mesh = scene->mMeshes[0];
    
    indices.reserve(mesh->mNumFaces * 3);
    for (unsigned int i = 0; i < mesh->mNumFaces; i++) {
        const aiFace& face = mesh->mFaces[i];
        for (unsigned int j = 0; j < face.mNumIndices; j++) {
            indices.push_back(face.mIndices[j]);
        }
    }
    
    return indices;
}

glm::vec3 ParticleSampler::calculateBoundingBoxMin(const std::vector<glm::vec3>& vertices) {
    glm::vec3 min(std::numeric_limits<float>::max());
    for (const auto& vertex : vertices) {
        min.x = std::min(min.x, vertex.x);
        min.y = std::min(min.y, vertex.y);
        min.z = std::min(min.z, vertex.z);
    }
    return min;
}

glm::vec3 ParticleSampler::calculateBoundingBoxMax(const std::vector<glm::vec3>& vertices) {
    glm::vec3 max(std::numeric_limits<float>::lowest());
    for (const auto& vertex : vertices) {
        max.x = std::max(max.x, vertex.x);
        max.y = std::max(max.y, vertex.y);
        max.z = std::max(max.z, vertex.z);
    }
    return max;
}

unsigned int ParticleSampler::createParticleVAO(const std::vector<Particle>& particles) {
    unsigned int VAO, VBO;
    
    std::vector<float> particleData;
    particleData.reserve(particles.size() * 3); // Only position data
    
    for (const auto& particle : particles) {
        particleData.push_back(particle.position.x);
        particleData.push_back(particle.position.y);
        particleData.push_back(particle.position.z);
    }
    
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    
    glBindVertexArray(VAO);
    
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, particleData.size() * sizeof(float), particleData.data(), GL_STATIC_DRAW);
    
    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    glBindVertexArray(0);
    
    return VAO;
}