#include "particleSplatter.h"
#include <iostream>
#include <glm/gtc/type_ptr.hpp>
#include "shader_utils.h" // Import your shader utilities

ParticleSplatter::ParticleSplatter() : VAO(0), VBO(0), instanceVBO(0), 
                                       shaderProgram(0), numParticles(0),
                                       radius(0.05f), smoothing(0.8f) {}

ParticleSplatter::~ParticleSplatter() {
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &instanceVBO);
    if (shaderProgram != 0) {
        glDeleteProgram(shaderProgram);
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
    
    std::cout << "Fluid splatter initialized with radius: " << radius << std::endl;
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
    
    // Draw particles as instanced quads
    glBindVertexArray(VAO);
    glDrawArraysInstanced(GL_TRIANGLES, 0, 6, numParticles);
    glBindVertexArray(0);
    
    // Restore previous OpenGL state
    if (!blendEnabled) glDisable(GL_BLEND);
    if (depthMaskEnabled) glDepthMask(GL_TRUE);
    glBlendFunc(blendSrcRGB, blendDestRGB);
}

// Create shaders using embedded code (fallback if external shaders not found)
unsigned int ParticleSplatter::createEmbeddedShaderProgram() {
    // Vertex shader compatible with your existing shader format
    const char* vertexShaderSource = R"(
        #version 330 core
out vec4 FragColor;

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoord;
in float DistFromCenter;

uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;
uniform vec3 objectColor;
uniform float smoothing = 0.8;

void main() {
    // Calculate distance from center of splat
    float dist = length(TexCoord - vec2(0.5));
    
    // Sharp cutoff at edge with minimal softening
    // We want overlapping but distinct particles
    float edge = 0.45;
    float mask = smoothstep(0.5, edge, dist);
    
    // Discard fragments outside the particle
    if (mask < 0.05)
        discard;
    
    // Much more uniform lighting across the particle
    // This is key to reducing the "ball" appearance
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    vec3 viewDir = normalize(viewPos - FragPos);
    
    // Minimal ambient variation - more uniform color
    float ambientStrength = 0.4;
    vec3 ambient = ambientStrength * lightColor;
    
    // Reduced diffuse lighting effect to make particles look less rounded
    float diff = max(dot(norm, lightDir), 0.0) * 0.4; // Reducing diffuse impact by 60%
    vec3 diffuse = diff * lightColor;
    
    // Very minimal specular to avoid "water droplet" look
    float specularStrength = 0.1; // Greatly reduced from 0.8
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(norm, halfwayDir), 0.0), 16.0); // Less sharp highlights
    vec3 specular = specularStrength * spec * lightColor;
    
    // Almost no depth variation within each particle
    // This helps create a more uniform fluid appearance
    float depthFactor = 0.95; // Almost constant value
    
    // Very minimal color variation from center to edge
    vec3 fluidColor = objectColor * (1.0 + (0.05 * (1.0 - dist))); // Only 5% variation
    
    // Combine lighting with minimal variations
    vec3 result = (ambient + diffuse + specular) * fluidColor * depthFactor;
    
    // Sharp alpha falloff at edges when particles overlap
    // This creates a surface tension effect rather than blurry blending
    float alphaFactor = smoothstep(0.0, 0.25, mask);
    float alpha = mix(0.9, 0.7, dist * 0.3) * alphaFactor; // Less transparency variation
    
    FragColor = vec4(result, alpha);
}
    )";
    
    // Fragment shader compatible with your existing shader format

const char* fragmentShaderSource = R"(
    #version 330 core
    out vec4 FragColor;
    
    in vec3 FragPos;
    in vec3 Normal;
    in vec2 TexCoord;
    in float DistFromCenter;
    
    uniform vec3 lightPos;
    uniform vec3 viewPos;
    uniform vec3 lightColor;
    uniform vec3 objectColor;
    uniform float smoothing = 0.8;
    
    void main() {
        // Calculate distance from center of splat
        float dist = length(TexCoord - vec2(0.5));
        
        // Improved circular mask with much softer edges
        float edge = mix(0.5, 0.25, smoothing); 
        float mask = smoothstep(0.5, edge, dist);
        
        // Extended influence - allows blending even beyond the normal radius
        float extendedMask = smoothstep(0.7, edge, dist);
        
        // Discard only the very distant fragments
        if (extendedMask < 0.001)
            discard;
        
        // Normal lighting calculation
        vec3 norm = normalize(Normal);
        vec3 lightDir = normalize(lightPos - FragPos);
        vec3 viewDir = normalize(viewPos - FragPos);
        
        // Ambient light
        float ambientStrength = 0.3;
        vec3 ambient = ambientStrength * lightColor;
        
        // Diffuse
        float diff = max(dot(norm, lightDir), 0.0);
        vec3 diffuse = diff * lightColor;
        
        // Specular (Blinn-Phong) 
        float specularStrength = 0.8;
        vec3 halfwayDir = normalize(lightDir + viewDir);
        float spec = pow(max(dot(norm, halfwayDir), 0.0), 128.0);
        vec3 specular = specularStrength * spec * lightColor;
        
        // Enhanced water-like effect with improved depth perception
        float depthFactor = mix(0.7, 1.0, 1.0 - dist * 1.8);
        
        // Improved color gradient from center to edge
        vec3 coreColor = objectColor * 1.1;
        vec3 edgeColor = mix(objectColor, vec3(0.9, 0.95, 1.0), 0.3);
        vec3 waterColor = mix(edgeColor, coreColor, smoothstep(0.3, 0.0, dist));
        
        // Final lighting calculation
        vec3 result = (ambient + diffuse + specular) * waterColor * depthFactor;
        
        // Much more gradual transparency from center to edge
        float baseAlpha = mix(0.95, 0.0, pow(dist * 1.2, 1.5));
        float finalAlpha = mix(mask, extendedMask, 0.6) * baseAlpha;
        
        // Add a slight glow effect to enhance the fluid appearance
        vec3 glowColor = mix(result, objectColor * 1.5, 0.1);
        result = mix(result, glowColor, smoothstep(0.4, 0.0, dist));
        
        FragColor = vec4(result, finalAlpha);
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