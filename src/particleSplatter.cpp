#include "particleSplatter.h"
#include <iostream>
#include <glm/gtc/type_ptr.hpp>
#include "shader_utils.h" // Import your shader utilities
#include <cmath>
#include <chrono>

ParticleSplatter::ParticleSplatter() 
    : VAO(0), VBO(0), instanceVBO(0), 
      shaderProgram(0), numParticles(0),
      radius(0.1f), smoothing(0.8f),
      surfaceTensionStrength(1.2f), blendDistance(0.15f),
      waterColor(0.2f, 0.5f, 0.8f),
      reflectionStrength(0.8f), specularPower(30.0f), specularStrength(2.0f),
      fresnelBias(0.25f), fresnelScale(1.0f), fresnelPower(5.0f),
      rippleStrength(0.15f), rippleSpeed(1.0f),
      reflectionsEnabled(true), refractionsEnabled(true),
      normalMapTexture(0), environmentMapTexture(0), texturesLoaded(false),
      currentTime(0.0f) {}

ParticleSplatter::~ParticleSplatter() {
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &instanceVBO);
    if (shaderProgram != 0) {
        glDeleteProgram(shaderProgram);
    }
    if (texturesLoaded) {
        glDeleteTextures(1, &normalMapTexture);
        glDeleteTextures(1, &environmentMapTexture);
    }
}

void ParticleSplatter::init(float particleRadius, float smoothingKernel) {
    radius = particleRadius;
    smoothing = smoothingKernel;
    
    // Load external shader files
    std::cout << "Compiling and linking water shader program..." << std::endl;
    shaderProgram = createShaderProgramEnhanced(
        "../src/assets/shaders/splatter_vertex_shader.glsl", 
        "../src/assets/shaders/splatter_fragment_shader.glsl"
    );
    
    if (shaderProgram == 0) {
        std::cerr << "ERROR: Failed to create water shader program!" << std::endl;
        return;
    }
    
    // Create a simple quad for each particle
    static const float quadVertices[] = {
        -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, // Added normals and texture coords
         0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f,
         0.5f,  0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f,
         0.5f,  0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f,
        -0.5f,  0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f,
        -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f
    };
    
    // Setup VBO for quad vertices
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &instanceVBO);
    
    glBindVertexArray(VAO);
    
    // Bind quad vertices with positions, normals, and texture coordinates
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
    
    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    // Normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    
    // Texture coordinate attribute
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
    
    // Prepare instanced buffer for particle data - includes position and velocity
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, 1000 * sizeof(glm::vec4) * 2, nullptr, GL_DYNAMIC_DRAW);
    
    // Attribute 3: Per-instance particle position data (position + size)
    glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4) * 2, (void*)0);
    glEnableVertexAttribArray(3);
    glVertexAttribDivisor(3, 1); // This makes it instanced
    
    // Attribute 4: Per-instance particle velocity data
    glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4) * 2, (void*)(sizeof(glm::vec4)));
    glEnableVertexAttribArray(4);
    glVertexAttribDivisor(4, 1); // This makes it instanced
    
    glBindVertexArray(0);
    
    // Create required textures
    normalMapTexture = createNormalMapTexture();
    environmentMapTexture = createEnvironmentMapTexture();
    texturesLoaded = (normalMapTexture != 0 && environmentMapTexture != 0);
    
    std::cout << "Water renderer initialized with radius: " << radius << std::endl;
    std::cout << "Water appearance: reflection = " << reflectionStrength 
              << ", specular power = " << specularPower 
              << ", fresnel = " << fresnelScale << std::endl;
    std::cout << "Surface tension = " << surfaceTensionStrength
              << ", blend distance = " << blendDistance << std::endl;
}

unsigned int ParticleSplatter::createNormalMapTexture() {
    // Size of the texture
    const int texWidth = 512;
    const int texHeight = 512;
    
    // Allocate memory for texture data
    unsigned char* texData = new unsigned char[texWidth * texHeight * 4]; // RGBA
    
    // Generate a procedural normal map for water surface
    for (int y = 0; y < texHeight; y++) {
        for (int x = 0; x < texWidth; x++) {
            float u = (float)x / texWidth;
            float v = (float)y / texHeight;
            
            // Generate a perlin-like noise pattern
            float scale1 = 5.0f, scale2 = 15.0f, scale3 = 25.0f;
            float noise1 = sin(u * scale1) * cos(v * scale1) * 0.5f;
            float noise2 = sin(u * scale2 + 0.3f) * cos(v * scale2 + 0.7f) * 0.3f;
            float noise3 = sin(u * scale3 + 1.2f) * cos(v * scale3 + 2.3f) * 0.2f;
            
            // Combine noises for varied surface
            float combinedNoise = noise1 + noise2 + noise3;
            
            // Calculate normal from the heightmap
            // This is a simplified approach - a real normal map would use gradient calculations
            float heightUp = combinedNoise - 0.01f;
            float heightDown = combinedNoise + 0.01f;
            float heightLeft = combinedNoise - 0.01f;
            float heightRight = combinedNoise + 0.01f;
            
            // Generate normal from height differences (simplified gradient)
            float nx = (heightRight - heightLeft) * 0.5f;
            float ny = (heightUp - heightDown) * 0.5f;
            float nz = 1.0f; // Strong Z component for subtle normals
            
            // Normalize
            float length = sqrt(nx*nx + ny*ny + nz*nz);
            nx /= length;
            ny /= length;
            nz /= length;
            
            // Convert from [-1, 1] to [0, 1] range for storage
            nx = nx * 0.5f + 0.5f;
            ny = ny * 0.5f + 0.5f;
            nz = nz * 0.5f + 0.5f;
            
            // Set RGBA values (RGB = normal, A = 1)
            int index = (y * texWidth + x) * 4;
            texData[index + 0] = (unsigned char)(nx * 255.0f); // R = normal.x
            texData[index + 1] = (unsigned char)(ny * 255.0f); // G = normal.y
            texData[index + 2] = (unsigned char)(nz * 255.0f); // B = normal.z
            texData[index + 3] = 255; // A = 1
        }
    }
    
    // Generate OpenGL texture
    unsigned int textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    
    // Set texture parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
    // Upload texture data
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texWidth, texHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData);
    glGenerateMipmap(GL_TEXTURE_2D);
    
    // Free memory
    delete[] texData;
    
    std::cout << "Water normal map created with ID: " << textureID << std::endl;
    return textureID;
}

unsigned int ParticleSplatter::createEnvironmentMapTexture() {
    // Size of the texture
    const int texWidth = 512;
    const int texHeight = 512;
    
    // Allocate memory for texture data
    unsigned char* texData = new unsigned char[texWidth * texHeight * 4]; // RGBA
    
    // Generate a simple gradient skybox for reflection
    for (int y = 0; y < texHeight; y++) {
        for (int x = 0; x < texWidth; x++) {
            float u = (float)x / texWidth;
            float v = (float)y / texHeight;
            
            // Create a simple gradient from blue to white (sky-like)
            float blueGradient = v; // 0 at bottom, 1 at top
            
            // Add some cloud-like patterns
            float cloudNoise = (sin(u * 10.0f) * 0.5f + 0.5f) * (cos(v * 8.0f) * 0.5f + 0.5f);
            cloudNoise *= (sin(u * 5.0f + v * 7.0f) * 0.5f + 0.5f);
            
            // Mix blue sky with white clouds
            float cloudMix = cloudNoise * 0.7f;
            float skyBlue = (1.0f - cloudMix) * blueGradient;
            
            // Set RGBA values
            int index = (y * texWidth + x) * 4;
            texData[index + 0] = (unsigned char)((0.5f + cloudMix * 0.5f) * 255.0f); // R
            texData[index + 1] = (unsigned char)((0.7f + cloudMix * 0.3f) * 255.0f); // G
            texData[index + 2] = (unsigned char)((0.9f + skyBlue * 0.1f) * 255.0f);  // B
            texData[index + 3] = 255; // A
        }
    }
    
    // Generate OpenGL texture
    unsigned int textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    
    // Set texture parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
    // Upload texture data
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texWidth, texHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData);
    glGenerateMipmap(GL_TEXTURE_2D);
    
    // Free memory
    delete[] texData;
    
    std::cout << "Environment map created with ID: " << textureID << std::endl;
    return textureID;
}

void ParticleSplatter::updateTime() {
    // Get current time for animations
    static auto startTime = std::chrono::high_resolution_clock::now();
    auto currentTimePoint = std::chrono::high_resolution_clock::now();
    float timeDiff = std::chrono::duration<float, std::chrono::seconds::period>(currentTimePoint - startTime).count();
    
    currentTime = timeDiff;
}

void ParticleSplatter::update(const std::vector<Particle>& particles) {
    numParticles = particles.size();
    if (numParticles == 0) return;
    
    // Update time for animations
    updateTime();
    
    // We store both position and velocity data for each particle
    particleData.resize(numParticles * 2);
    
    // Update particle data: position, size, velocity, and padding
    for (size_t i = 0; i < numParticles; i++) {
        // Add small random offsets to Z position to avoid Z-fighting
        float randomOffset = ((float)rand() / (float)RAND_MAX) * 0.0001f;
        
        // Store position with radius
        particleData[i*2] = glm::vec4(
            particles[i].position.x, 
            particles[i].position.y, 
            particles[i].position.z + randomOffset, 
            radius
        );
        
        // Store velocity with a padding value
        particleData[i*2 + 1] = glm::vec4(
            particles[i].velocity.x,
            particles[i].velocity.y,
            particles[i].velocity.z,
            0.0f // Padding
        );
    }
    
    // Update instance buffer
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, particleData.size() * sizeof(glm::vec4), particleData.data(), GL_DYNAMIC_DRAW);
}

void ParticleSplatter::draw(const glm::mat4& model, const glm::mat4& view, const glm::mat4& projection) {
    if (numParticles == 0 || shaderProgram == 0) return;
    
    // Save previous OpenGL state
    GLboolean blendEnabled, depthMaskEnabled;
    GLint blendSrcRGB, blendDestRGB;
    glGetBooleanv(GL_BLEND, &blendEnabled);
    glGetBooleanv(GL_DEPTH_WRITEMASK, &depthMaskEnabled);
    glGetIntegerv(GL_BLEND_SRC_RGB, &blendSrcRGB);
    glGetIntegerv(GL_BLEND_DST_RGB, &blendDestRGB);
    
    // Enable transparency and adjust blending for water-like appearance
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE); // Disable depth writes for transparency
    
    // Sort particles back-to-front for better transparency
    if (numParticles > 1) {
        // Extract camera position from view matrix
        glm::mat4 invView = glm::inverse(view);
        glm::vec3 cameraPos = glm::vec3(invView[3]);
        
        // Sort particles by distance to camera (back to front)
        std::vector<std::pair<float, size_t>> distances(numParticles);
        for (size_t i = 0; i < numParticles; i++) {
            glm::vec3 particlePos = glm::vec3(particleData[i*2]);
            distances[i] = std::make_pair(-glm::distance(particlePos, cameraPos), i);
        }
        std::sort(distances.begin(), distances.end());
        
        // Reorder particle data
        std::vector<glm::vec4> sortedData(particleData.size());
        for (size_t i = 0; i < numParticles; i++) {
            sortedData[i*2] = particleData[distances[i].second*2];
            sortedData[i*2 + 1] = particleData[distances[i].second*2 + 1];
        }
        
        // Update the buffer with sorted data
        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferData(GL_ARRAY_BUFFER, sortedData.size() * sizeof(glm::vec4), sortedData.data(), GL_DYNAMIC_DRAW);
    }
    
    // Use the water shader program
    glUseProgram(shaderProgram);
    
    // Set common uniforms
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    
    // Set time for animations
    GLint timeLoc = glGetUniformLocation(shaderProgram, "time");
    if (timeLoc != -1) {
        glUniform1f(timeLoc, currentTime);
    }
    
    // Set camera position for reflections
    GLint viewPosLoc = glGetUniformLocation(shaderProgram, "viewPos");
    if (viewPosLoc != -1) {
        glm::mat4 invView = glm::inverse(view);
        glm::vec3 cameraPos = glm::vec3(invView[3]);
        glUniform3fv(viewPosLoc, 1, glm::value_ptr(cameraPos));
    }
    
    // Set lighting parameters
    GLint lightPosLoc = glGetUniformLocation(shaderProgram, "lightPos");
    if (lightPosLoc != -1) {
        glUniform3f(lightPosLoc, 5.0f, 5.0f, 5.0f); // Light position
    }
    
    GLint lightColorLoc = glGetUniformLocation(shaderProgram, "lightColor");
    if (lightColorLoc != -1) {
        glUniform3f(lightColorLoc, 1.0f, 1.0f, 1.0f); // White light
    }
    
    // Set water color
    GLint waterColorLoc = glGetUniformLocation(shaderProgram, "waterColor");
    if (waterColorLoc != -1) {
        glUniform3fv(waterColorLoc, 1, glm::value_ptr(waterColor));
    }
    
    // Set reflection and refraction parameters
    GLint reflectionStrengthLoc = glGetUniformLocation(shaderProgram, "reflectionStrength");
    if (reflectionStrengthLoc != -1) {
        glUniform1f(reflectionStrengthLoc, reflectionStrength);
    }
    
    GLint specularPowerLoc = glGetUniformLocation(shaderProgram, "specularPower");
    if (specularPowerLoc != -1) {
        glUniform1f(specularPowerLoc, specularPower);
    }
    
    GLint specularStrengthLoc = glGetUniformLocation(shaderProgram, "specularStrength");
    if (specularStrengthLoc != -1) {
        glUniform1f(specularStrengthLoc, specularStrength);
    }
    
    GLint fresnelBiasLoc = glGetUniformLocation(shaderProgram, "fresnelBias");
    if (fresnelBiasLoc != -1) {
        glUniform1f(fresnelBiasLoc, fresnelBias);
    }
    
    GLint fresnelScaleLoc = glGetUniformLocation(shaderProgram, "fresnelScale");
    if (fresnelScaleLoc != -1) {
        glUniform1f(fresnelScaleLoc, fresnelScale);
    }
    
    GLint fresnelPowerLoc = glGetUniformLocation(shaderProgram, "fresnelPower");
    if (fresnelPowerLoc != -1) {
        glUniform1f(fresnelPowerLoc, fresnelPower);
    }
    
    // Set surface ripple parameters
    GLint rippleStrengthLoc = glGetUniformLocation(shaderProgram, "rippleStrength");
    if (rippleStrengthLoc != -1) {
        glUniform1f(rippleStrengthLoc, rippleStrength);
    }
    
    GLint rippleSpeedLoc = glGetUniformLocation(shaderProgram, "rippleSpeed");
    if (rippleSpeedLoc != -1) {
        glUniform1f(rippleSpeedLoc, rippleSpeed);
    }
    
    // Set surface tension and blending parameters
    GLint surfaceTensionStrengthLoc = glGetUniformLocation(shaderProgram, "surfaceTensionStrength");
    if (surfaceTensionStrengthLoc != -1) {
        glUniform1f(surfaceTensionStrengthLoc, surfaceTensionStrength);
    }
    
    GLint blendDistanceLoc = glGetUniformLocation(shaderProgram, "blendDistance");
    if (blendDistanceLoc != -1) {
        glUniform1f(blendDistanceLoc, blendDistance);
    }
    
    // Set feature toggles
    GLint reflectionsEnabledLoc = glGetUniformLocation(shaderProgram, "reflectionsEnabled");
    if (reflectionsEnabledLoc != -1) {
        glUniform1i(reflectionsEnabledLoc, reflectionsEnabled ? 1 : 0);
    }
    
    GLint refractionsEnabledLoc = glGetUniformLocation(shaderProgram, "refractionsEnabled");
    if (refractionsEnabledLoc != -1) {
        glUniform1i(refractionsEnabledLoc, refractionsEnabled ? 1 : 0);
    }
    
    // Pass all particle positions for surface tension and blending
    GLint numParticlesLoc = glGetUniformLocation(shaderProgram, "numParticles");
    if (numParticlesLoc != -1) {
        // Limit number of particles for calculation for performance
        int maxParticles = std::min((int)numParticles, 100);
        glUniform1i(numParticlesLoc, maxParticles);
        
        GLint particlePositionsLoc = glGetUniformLocation(shaderProgram, "particlePositions");
        if (particlePositionsLoc != -1) {
            // Pass particle positions for surface tension calculation
            std::vector<glm::vec4> positionsOnly(maxParticles);
            for (int i = 0; i < maxParticles; i++) {
                positionsOnly[i] = particleData[i*2]; // Every other entry is a position
            }
            glUniform4fv(particlePositionsLoc, maxParticles, glm::value_ptr(positionsOnly[0]));
        }
    }
    
    // Bind textures
    if (texturesLoaded) {
        // Bind normal map texture
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, normalMapTexture);
        GLint normalMapLoc = glGetUniformLocation(shaderProgram, "normalMap");
        if (normalMapLoc != -1) {
            glUniform1i(normalMapLoc, 0); // Texture unit 0
        }
        
        // Bind environment map texture
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, environmentMapTexture);
        GLint environmentMapLoc = glGetUniformLocation(shaderProgram, "environmentMap");
        if (environmentMapLoc != -1) {
            glUniform1i(environmentMapLoc, 1); // Texture unit 1
        }
    }
    
    // Draw particles as instanced quads
    glBindVertexArray(VAO);
    glDrawArraysInstanced(GL_TRIANGLES, 0, 6, numParticles);
    glBindVertexArray(0);
    
    // Restore previous OpenGL state
    if (!blendEnabled) glDisable(GL_BLEND);
    if (depthMaskEnabled) glDepthMask(GL_TRUE);
    glBlendFunc(blendSrcRGB, blendDestRGB);
}