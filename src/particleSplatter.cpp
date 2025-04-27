#include "particleSplatter.h"
#include <iostream>
#include <glm/gtc/type_ptr.hpp>
#include "shader_utils.h" // Import your shader utilities
#include <cmath>

ParticleSplatter::ParticleSplatter() 
    : VAO(0), VBO(0), instanceVBO(0), 
      shaderProgram(0), numParticles(0),
      radius(0.05f), smoothing(0.8f),
      metaballThreshold(1.0f), metaballStrength(0.5f),
      dropletScale(15.0f), dropletIntensity(0.5f),
      waterDropletsEnabled(true), 
      waterTexture(0), waterTextureLoaded(false) {}

ParticleSplatter::~ParticleSplatter() {
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &instanceVBO);
    if (shaderProgram != 0) {
        glDeleteProgram(shaderProgram);
    }
    if (waterTextureLoaded) {
        glDeleteTextures(1, &waterTexture);
    }
}

void ParticleSplatter::init(float particleRadius, float smoothingKernel) {
    radius = particleRadius;
    smoothing = smoothingKernel;
    
    // Load external shader files
    std::cout << "Compiling and linking fluid splatter shader program..." << std::endl;
    shaderProgram = createShaderProgramEnhanced(
        "../src/assets/shaders/splatter_vertex_shader.glsl", 
        "../src/assets/shaders/splatter_fragment_shader.glsl"
    );
    
    if (shaderProgram == 0) {
        std::cerr << "ERROR: Failed to create fluid splatter shader program!" << std::endl;
        return;
    }
    
    // Create a simple quad for each particle
    static const float quadVertices[] = {
        -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f, // Added normals pointing in Z direction
         0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
         0.5f,  0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
         0.5f,  0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
        -0.5f,  0.5f, 0.0f, 0.0f, 0.0f, 1.0f,
        -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 1.0f
    };
    
    // Setup VBO for quad vertices
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &instanceVBO);
    
    glBindVertexArray(VAO);
    
    // Bind quad vertices with positions and normals
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
    
    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    
    // Normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    
    // Prepare instanced buffer for particle data - now includes velocity
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, 1000 * sizeof(glm::vec4) * 2, nullptr, GL_DYNAMIC_DRAW); // Reserve space for position and velocity
    
    // Attribute 2: Per-instance particle position data (position + size)
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4) * 2, (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1); // This makes it instanced
    
    // Attribute 3: Per-instance particle velocity data
    glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4) * 2, (void*)(sizeof(glm::vec4)));
    glEnableVertexAttribArray(3);
    glVertexAttribDivisor(3, 1); // This makes it instanced
    
    glBindVertexArray(0);
    
    // Create water droplet texture if needed
    if (waterDropletsEnabled) {
        waterTexture = createWaterDropletTexture();
        waterTextureLoaded = (waterTexture != 0);
    }
    
    std::cout << "Fluid splatter initialized with radius: " << radius << std::endl;
    std::cout << "Metaball effect: threshold = " << metaballThreshold << ", strength = " << metaballStrength << std::endl;
    std::cout << "Water droplets: " << (waterDropletsEnabled ? "enabled" : "disabled") 
              << ", scale = " << dropletScale 
              << ", intensity = " << dropletIntensity << std::endl;
}

unsigned int ParticleSplatter::createWaterDropletTexture() {
    // Size of the texture
    const int texWidth = 512;
    const int texHeight = 512;
    
    // Allocate memory for texture data
    unsigned char* texData = new unsigned char[texWidth * texHeight * 4]; // RGBA
    
    // Generate a procedural water droplet texture
    for (int y = 0; y < texHeight; y++) {
        for (int x = 0; x < texWidth; x++) {
            float u = (float)x / texWidth;
            float v = (float)y / texHeight;
            
            // Center coordinates
            float cx = u - 0.5f;
            float cy = v - 0.5f;
            
            // Distance from center
            float dist = sqrt(cx*cx + cy*cy);
            
            // Create a basic droplet pattern
            float droplet = std::max(0.0f, 1.0f - dist * 2.0f);
            droplet = pow(droplet, 2.0f); // Shape the falloff
            
            // Add some noise for texture
            float noise = (float)rand() / RAND_MAX * 0.1f;
            
            // Add a highlight effect
            float highlight = std::max(0.0f, 1.0f - dist * 3.0f);
            highlight = pow(highlight, 5.0f); // Sharp highlight
            
            // Final combined effect
            float alpha = std::min(1.0f, droplet + highlight);
            
            // Set RGBA values
            int index = (y * texWidth + x) * 4;
            texData[index + 0] = (unsigned char)(255 * (0.8f + 0.2f * droplet)); // R
            texData[index + 1] = (unsigned char)(255 * (0.9f + 0.1f * droplet)); // G
            texData[index + 2] = (unsigned char)(255 * (1.0f)); // B
            texData[index + 3] = (unsigned char)(255 * alpha); // A
        }
    }
    
    // Generate OpenGL texture
    unsigned int textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    
    // Set texture parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    
    // Upload texture data
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texWidth, texHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData);
    glGenerateMipmap(GL_TEXTURE_2D);
    
    // Free memory
    delete[] texData;
    
    std::cout << "Water droplet texture created with ID: " << textureID << std::endl;
    return textureID;
}

void ParticleSplatter::update(const std::vector<Particle>& particles) {
    numParticles = particles.size();
    if (numParticles == 0) return;
    
    // We now store both position and velocity data for each particle
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
    
    // Enable transparency and disable depth writes
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE); // Disable depth writes for transparency
    
    // Sort particles back-to-front for better transparency (simple approach)
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
    
    // Use the splatter shader program
    glUseProgram(shaderProgram);
    
    // Set common uniforms
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    
    // Set lighting uniforms if using standard shaders
    GLint lightPosLoc = glGetUniformLocation(shaderProgram, "lightPos");
    if (lightPosLoc != -1) {
        glUniform3f(lightPosLoc, 1.2f, 1.0f, 2.0f); // Match your existing light
    }
    
    GLint viewPosLoc = glGetUniformLocation(shaderProgram, "viewPos");
    if (viewPosLoc != -1) {
        glm::mat4 invView = glm::inverse(view);
        glm::vec3 cameraPos = glm::vec3(invView[3]);
        glUniform3fv(viewPosLoc, 1, glm::value_ptr(cameraPos));
    }
    
    GLint lightColorLoc = glGetUniformLocation(shaderProgram, "lightColor");
    if (lightColorLoc != -1) {
        glUniform3f(lightColorLoc, 1.0f, 1.0f, 1.0f); // White light
    }
    
    // Set fluid color
    GLint objectColorLoc = glGetUniformLocation(shaderProgram, "objectColor");
    if (objectColorLoc != -1) {
        glUniform3f(objectColorLoc, 0.2f, 0.6f, 0.9f); // Fluid blue color
    }
    
    // Set fluid-specific uniforms
    GLint smoothingLoc = glGetUniformLocation(shaderProgram, "smoothing");
    if (smoothingLoc != -1) {
        glUniform1f(smoothingLoc, smoothing);
    }
    
    // Set droplet deformation factor
    GLint dropletDeformFactorLoc = glGetUniformLocation(shaderProgram, "dropletDeformFactor");
    if (dropletDeformFactorLoc != -1) {
        glUniform1f(dropletDeformFactorLoc, 0.7f); // Controls how much to stretch particles based on velocity
    }
    
    // Set metaball-specific uniforms
    GLint metaballThresholdLoc = glGetUniformLocation(shaderProgram, "metaballThreshold");
    if (metaballThresholdLoc != -1) {
        glUniform1f(metaballThresholdLoc, metaballThreshold);
    }
    
    GLint metaballStrengthLoc = glGetUniformLocation(shaderProgram, "metaballStrength");
    if (metaballStrengthLoc != -1) {
        glUniform1f(metaballStrengthLoc, metaballStrength);
    }
    
    // Set water droplet-specific uniforms
    GLint enableWaterDropletsLoc = glGetUniformLocation(shaderProgram, "enableWaterDroplets");
    if (enableWaterDropletsLoc != -1) {
        glUniform1i(enableWaterDropletsLoc, waterDropletsEnabled ? 1 : 0);
    }
    
    GLint dropletScaleLoc = glGetUniformLocation(shaderProgram, "dropletScale");
    if (dropletScaleLoc != -1) {
        glUniform1f(dropletScaleLoc, dropletScale);
    }
    
    GLint dropletIntensityLoc = glGetUniformLocation(shaderProgram, "dropletIntensity");
    if (dropletIntensityLoc != -1) {
        glUniform1f(dropletIntensityLoc, dropletIntensity);
    }
    
    // Pass all particle positions for metaball calculation
    GLint numParticlesLoc = glGetUniformLocation(shaderProgram, "numParticles");
    if (numParticlesLoc != -1) {
        // Limit number of particles for metaball calculation for performance
        int maxMetaballParticles = std::min((int)numParticles, 100);
        glUniform1i(numParticlesLoc, maxMetaballParticles);
        
        GLint particlePositionsLoc = glGetUniformLocation(shaderProgram, "particlePositions");
        if (particlePositionsLoc != -1) {
            // Pass every other entry (just positions, not velocities)
            std::vector<glm::vec4> positionsOnly(maxMetaballParticles);
            for (int i = 0; i < maxMetaballParticles; i++) {
                positionsOnly[i] = particleData[i*2]; // Every other entry is a position
            }
            // Pass the first 100 (or fewer) particle positions for metaball calculation
            glUniform4fv(particlePositionsLoc, maxMetaballParticles, glm::value_ptr(positionsOnly[0]));
        }
    }
    
    // Bind water droplet texture if enabled
    if (waterDropletsEnabled && waterTextureLoaded) {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, waterTexture);
        
        GLint waterTextureLoc = glGetUniformLocation(shaderProgram, "waterTexture");
        if (waterTextureLoc != -1) {
            glUniform1i(waterTextureLoc, 0); // Texture unit 0
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