#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <glad/glad.h>

#include <Eigen/Dense>
#include <limits>
#include <random>

using scalar = float;
using Vector3 = Eigen::Matrix<scalar, 3, 1>;

class Mesh {
public:
    Mesh(const std::string& path);
    ~Mesh();

    void draw() const;
    void bind() const;
    void unBind() const;

    std::vector<float> vertices;
    std::vector<unsigned int> indices;

    // Samples evenly spaced points on the surface using Poisson Disk sampling
    std::vector<Vector3> sampleSurfacePoints(
        scalar minRadius = 0.05f,
        unsigned int numTrials = 60,
        scalar initialDensity = 5.0f,
        unsigned int distanceNorm = 2
    ) const;

    // Samples points inside the volume using point-in-mesh testing
    std::vector<Vector3> sampleVolumePoints(
        unsigned int numPoints = 1000
    ) const;

private:
    unsigned int VAO, VBO, EBO;
    void setupBuffers();

    // Internal helper for ray-triangle intersection (used for volume sampling)
    bool rayIntersectsTriangle(
        const Vector3& orig,
        const Vector3& dir,
        const Vector3& v0,
        const Vector3& v1,
        const Vector3& v2
    ) const;
};

#endif // MESH_H