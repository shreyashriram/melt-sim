#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoords;
layout (location = 3) in vec4 aParticleData; // xyz = position, w = size
layout (location = 4) in vec4 aParticleVelocity; // xyz = velocity, w = unused

out vec3 FragPos;
out vec3 Normal;
out vec2 TexCoord;
out vec3 Tangent;
out vec3 Bitangent;
out vec3 ParticleVelocity;
out vec4 ParticleParams;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform float time;
uniform float rippleSpeed;
uniform float rippleStrength;

void main() {
    // Extract particle data
    vec3 particlePos = aParticleData.xyz;
    float particleSize = aParticleData.w;
    vec3 velocity = aParticleVelocity.xyz;
    
    // Calculate billboarding vectors - aligned with camera
    vec3 camRight = normalize(vec3(view[0][0], view[1][0], view[2][0]));
    vec3 camUp = normalize(vec3(view[0][1], view[1][1], view[2][1]));
    vec3 billboardNormal = normalize(vec3(view[0][2], view[1][2], view[2][2]));
    
    // Calculate velocity magnitude for deformation
    float velocityMag = length(velocity);
    
    // Add ripple effect based on time
    float waveOffset = 0.0;
    if (velocityMag < 0.1) { // Only add ripples to relatively stationary particles
        // Simple sine wave displacement based on distance from center
        float distFromCenter = length(aPos.xy);
        waveOffset = sin(distFromCenter * 10.0 + time * rippleSpeed) * rippleStrength * 
                     max(0.0, 1.0 - distFromCenter * 2.0); // Fade out toward edges
    }
    
    // Generate billboarded vertex position with ripple effect
    vec3 vertPos = particlePos + 
                  camRight * aPos.x * particleSize + 
                  camUp * aPos.y * particleSize + 
                  billboardNormal * waveOffset; // Apply wave offset along normal
    
    // Tangent and bitangent for normal mapping
    // In a billboard, tangent is along the right vector and bitangent along up
    Tangent = camRight;
    Bitangent = camUp;
    
    // Transform to world space
    FragPos = vec3(model * vec4(vertPos, 1.0));
    
    // Normal should point toward camera for the billboard
    Normal = mat3(transpose(inverse(model))) * billboardNormal;
    
    // Pass texture coordinates for the water surface
    TexCoord = aTexCoords;
    
    // Pass additional data to fragment shader
    ParticleVelocity = velocity;
    ParticleParams = vec4(particlePos, particleSize);
    
    gl_Position = projection * view * model * vec4(vertPos, 1.0);
}