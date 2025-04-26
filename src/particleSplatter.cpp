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
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aNormal;
        layout (location = 2) in vec4 aParticleData; // xyz = position, w = size
        
        out vec3 FragPos;
        out vec3 Normal;
        out vec2 TexCoord;
        out float DistFromCenter;
        
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
            
            gl_Position = projection * view * model * vec4(vertPos, 1.0);
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
            
            // Circular mask with smoothed edges
            float edge = mix(0.5, 0.35, smoothing); // Adjustable edge softness
            float mask = smoothstep(0.5, edge, dist);
            
            // Discard fragments outside the circle
            if (mask < 0.01)
                discard;
            
            // Normal lighting calculation (similar to your existing shader)
            vec3 norm = normalize(Normal);
            vec3 lightDir = normalize(lightPos - FragPos);
            vec3 viewDir = normalize(viewPos - FragPos);
            
            // Ambient
            float ambientStrength = 0.2;
            vec3 ambient = ambientStrength * lightColor;
            
            // Diffuse
            float diff = max(dot(norm, lightDir), 0.0);
            vec3 diffuse = diff * lightColor;
            
            // Specular (Blinn-Phong)
            float specularStrength = 0.6;
            vec3 halfwayDir = normalize(lightDir + viewDir);
            float spec = pow(max(dot(norm, halfwayDir), 0.0), 32.0);
            vec3 specular = specularStrength * spec * lightColor;
            
            // Apply water-like effect (refraction and color depth)
            float depthFactor = mix(0.7, 1.0, 1.0 - dist * 1.5);
            vec3 waterColor = mix(objectColor, objectColor * 1.3, 1.0 - dist);
            
            vec3 result = (ambient + diffuse + specular) * waterColor * depthFactor;
            
            // Apply transparency gradient from center to edge
            float alpha = mask * mix(1.0, 0.85, dist * 2.0);
            
            FragColor = vec4(result, alpha);
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