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
    
    // Try to load external shader files first (if they exist)
    std::cout << "Compiling and linking fluid splatter shader program..." << std::endl;
    shaderProgram = createShaderProgramEnhanced(
        "../src/assets/shaders/splatter_vertex_shader.glsl", 
        "../src/assets/shaders/splatter_fragment_shader.glsl"
    );
    
    if (shaderProgram == 0) {
        // If external shader loading failed, create embedded shader program
        std::cout << "Using embedded splatter shaders..." << std::endl;
        shaderProgram = createEmbeddedShaderProgram();
        if (shaderProgram == 0) {
            std::cerr << "ERROR: Failed to create fluid splatter shader program!" << std::endl;
            return;
        }
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
    
    // Prepare instanced buffer for particle data
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, 1000 * sizeof(glm::vec4), nullptr, GL_DYNAMIC_DRAW); // Reserve space
    
    // Attribute 2: Per-instance particle data (position + size)
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), (void*)0);
    glEnableVertexAttribArray(2);
    glVertexAttribDivisor(2, 1); // This makes it instanced
    
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
    
    particleData.resize(numParticles);
    
    // Update particle data: position (xyz) and size (w)
    for (size_t i = 0; i < numParticles; i++) {
        // Add small random offsets to Z position to avoid Z-fighting
        float randomOffset = ((float)rand() / (float)RAND_MAX) * 0.0001f;
        
        // Store position with radius
        particleData[i] = glm::vec4(
            particles[i].position.x, 
            particles[i].position.y, 
            particles[i].position.z + randomOffset, 
            radius
        );
    }
    
    // Update instance buffer
    glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
    glBufferData(GL_ARRAY_BUFFER, numParticles * sizeof(glm::vec4), particleData.data(), GL_DYNAMIC_DRAW);
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
            glm::vec3 particlePos = glm::vec3(particleData[i]);
            distances[i] = std::make_pair(-glm::distance(particlePos, cameraPos), i);
        }
        std::sort(distances.begin(), distances.end());
        
        // Reorder particle data
        std::vector<glm::vec4> sortedData(numParticles);
        for (size_t i = 0; i < numParticles; i++) {
            sortedData[i] = particleData[distances[i].second];
        }
        
        // Update the buffer with sorted data
        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferData(GL_ARRAY_BUFFER, numParticles * sizeof(glm::vec4), sortedData.data(), GL_DYNAMIC_DRAW);
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
            // Pass the first 100 (or fewer) particle positions for metaball calculation
            glUniform4fv(particlePositionsLoc, maxMetaballParticles, glm::value_ptr(particleData[0]));
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

// Create updated shaders using embedded code (fallback if external shaders not found)
unsigned int ParticleSplatter::createEmbeddedShaderProgram() {
    // Vertex shader with metaball and water droplet support
    const char* vertexShaderSource = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aNormal;
        layout (location = 2) in vec4 aParticleData; // xyz = position, w = size

        out vec3 FragPos;
        out vec3 Normal;
        out vec2 TexCoord;
        out float DistFromCenter;
        out vec4 ParticleParams; // Pass particle data to fragment shader

        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;

        void main() {
            // Extract particle data
            vec3 particlePos = aParticleData.xyz;
            float particleSize = aParticleData.w;
            
            // Calculate billboarding vectors
            vec3 camRight = normalize(vec3(view[0][0], view[1][0], view[2][0]));
            vec3 camUp = normalize(vec3(view[0][1], view[1][1], view[2][1]));
            
            // Generate billboarded vertex position
            vec3 vertPos = particlePos + camRight * aPos.x * particleSize + camUp * aPos.y * particleSize;
            
            // Use normal perpendicular to camera plane (for lighting)
            vec3 faceNormal = -normalize(vec3(view[0][2], view[1][2], view[2][2]));
            
            // We need to convert these to world space for our existing shader format
            FragPos = vec3(model * vec4(vertPos, 1.0));
            Normal = mat3(transpose(inverse(model))) * faceNormal;
            
            // Pass texture coordinates for the splat
            TexCoord = aPos.xy + 0.5; // Convert from [-0.5, 0.5] to [0, 1]
            
            // Calculate distance from center for the fragment shader
            DistFromCenter = length(aPos.xy);
            
            // Pass particle parameters to fragment shader for unique water effects
            ParticleParams = aParticleData;
            
            gl_Position = projection * view * model * vec4(vertPos, 1.0);
        }
    )";
    const char* fragmentShaderSource = R"(
        #version 330 core
        out vec4 FragColor;

        in vec3 FragPos;
        in vec3 Normal;
        in vec2 TexCoord;
        in float DistFromCenter;
        in vec4 ParticleParams; // xyz = position, w = size

        // Metaball parameters
        uniform int numParticles;         // Total number of particles
        uniform vec4 particlePositions[100]; // Array of particle positions (increase if needed)
        uniform float metaballThreshold = 1.0;  // Threshold for metaball effect
        uniform float metaballStrength = 0.5;   // Strength of metaball effect

        // Standard lighting parameters
        uniform vec3 lightPos;
        uniform vec3 viewPos;
        uniform vec3 lightColor;
        uniform vec3 objectColor;
        uniform float smoothing = 0.8;

        // Water droplet parameters
        uniform bool enableWaterDroplets = true;
        uniform float dropletScale = 15.0;     // Scale of water droplet pattern
        uniform float dropletIntensity = 0.5;  // Intensity of the water droplet effect
        uniform sampler2D waterTexture;        // Optional texture for water droplets

        // Random functions for noise
        float random(vec2 st) {
            return fract(sin(dot(st.xy, vec2(12.9898, 78.233))) * 43758.5453123);
        }

        float noise(vec2 st) {
            vec2 i = floor(st);
            vec2 f = fract(st);
            
            // Four corners in 2D of a tile
            float a = random(i);
            float b = random(i + vec2(1.0, 0.0));
            float c = random(i + vec2(0.0, 1.0));
            float d = random(i + vec2(1.0, 1.0));

            // Smooth interpolation
            vec2 u = f * f * (3.0 - 2.0 * f);

            return mix(a, b, u.x) + 
                (c - a) * u.y * (1.0 - u.x) + 
                (d - b) * u.x * u.y;
        }

        // Function to create water droplet pattern
        float waterDropletPattern(vec2 uv, float scale, float time) {
            float n = noise(uv * scale);
            
            // Create water droplet effect by adding multiple noise layers
            float droplet = smoothstep(0.4, 0.5, n);
            
            // Add some variation for more natural look
            droplet += 0.1 * noise(uv * scale * 2.0 + vec2(time * 0.1));
            
            return droplet;
        }

        void main() {
            // Calculate distance from center of splat
            float dist = length(TexCoord - vec2(0.5));
            
            // Basic mask with sharp cutoff
            float edge = 0.45;
            float mask = smoothstep(0.5, edge, dist);
            
            // Discard fragments outside the particle
            if (mask < 0.05)
                discard;
            
            // Standard lighting calculation
            vec3 norm = normalize(Normal);
            vec3 lightDir = normalize(lightPos - FragPos);
            vec3 viewDir = normalize(viewPos - FragPos);
            
            float ambientStrength = 0.4;
            vec3 ambient = ambientStrength * lightColor;
            
            float diff = max(dot(norm, lightDir), 0.0) * 0.4;
            vec3 diffuse = diff * lightColor;
            
            float specularStrength = 0.1;
            vec3 halfwayDir = normalize(lightDir + viewDir);
            float spec = pow(max(dot(norm, halfwayDir), 0.0), 16.0);
            vec3 specular = specularStrength * spec * lightColor;
            
            // Generate water droplet texture
            float droplet = 0.0;
            if (enableWaterDroplets) {
                // Use particle position as a unique seed for variation
                vec2 seed = vec2(ParticleParams.x * 10.0, ParticleParams.y * 10.0);
                
                // Create water droplet pattern
                droplet = waterDropletPattern(TexCoord * 3.0 + seed, dropletScale, 0.0);
                droplet *= dropletIntensity;
            }
            
            // Apply water droplet effect to normal and color
            vec3 adjustedNormal = normalize(norm + vec3(droplet * 0.1, droplet * 0.1, 0.0));
            
            // Metaball calculation - influences the alpha
            float metaballField = 0.0;
            // Self contribution
            metaballField += 1.0 / (dist * dist + 0.01) * metaballStrength;
            
            // Apply the base color with water droplet effect
            float dropletHighlight = droplet * 0.2;
            vec3 fluidColor = objectColor * (1.0 + (0.05 * (1.0 - dist) + dropletHighlight));
            
            // Enhanced fluid appearance with varying depth and color
            float depthFactor = 0.95 - droplet * 0.1;
            
            // Combine lighting
            vec3 result = (ambient + diffuse + specular) * fluidColor * depthFactor;
            
            // Enhanced alpha for metaball effect
            float alphaFactor = smoothstep(0.0, 0.25, mask);
            float baseAlpha = mix(0.9, 0.7, dist * 0.3) * alphaFactor;
            
            // Add subtle refraction at the edges
            float refractionFactor = 0.05 * (1.0 - mask) * (1.0 + droplet);
            result += refractionFactor * vec3(0.8, 0.9, 1.0);
            
            // Add subtle water caustics effect
            float causticEffect = noise(TexCoord * 20.0) * noise(TexCoord * 15.0 + vec2(0.2, 0.3));
            result += causticEffect * 0.1 * baseAlpha * vec3(0.8, 0.95, 1.0);
            
            // Final color
            FragColor = vec4(result, baseAlpha);
        }
    )";
    
    // Compile and link the shaders
    unsigned int vertexShader, fragmentShader;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    
    // Compile shaders and check for errors
    int success;
    char infoLog[512];
    
    glCompileShader(vertexShader);
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
        return 0;
    }
    
    glCompileShader(fragmentShader);
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
        return 0;
    }
    
    // Link shaders
    unsigned int program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);
    
    // Check for linking errors
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(program, 512, NULL, infoLog);
        std::cerr << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
        return 0;
    }
    
    // Shaders are linked, we can delete them
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    
    return program;
}
    // Fragment shader with metaball and water droplet effects
    