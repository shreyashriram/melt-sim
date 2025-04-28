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

// Hash function for pseudo-random values
float hash(vec2 p) {
    p = fract(p * vec2(123.34, 456.21));
    p += dot(p, p + 45.32);
    return fract(p.x * p.y);
}

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
    
    // Calculate squash and stretch factors based on velocity
    float stretchFactor = 1.0;
    float squashFactor = 1.0;
    vec3 stretchDir = vec3(0.0, 0.0, 0.0);
    
    if (velocityMag > 0.1) {
        // Normalize velocity to get direction
        stretchDir = normalize(velocity);
        
        // Calculate stretch factor - increase along velocity direction
        // Enhanced to make stretching more pronounced
        stretchFactor = 1.0 + min(velocityMag * 0.7, 1.5);
        
        // Calculate squash factor - decrease perpendicular to velocity
        // Conserve volume by making squash = 1/sqrt(stretch)
        squashFactor = 1.0 / sqrt(stretchFactor);
        
        // Add subtle oscillation to the stretch factor for more dynamic movement
        stretchFactor *= 1.0 + sin(time * 5.0 * velocityMag) * 0.05;
    }
    
    // Project stretch direction onto billboard plane
    vec3 stretchProjRight = dot(stretchDir, camRight) * camRight;
    vec3 stretchProjUp = dot(stretchDir, camUp) * camUp;
    
    // Calculate modified right and up vectors for the billboard
    vec3 modifiedRight = camRight;
    vec3 modifiedUp = camUp;
    
    if (velocityMag > 0.1) {
        // Apply stretch along velocity direction and squash perpendicular to it
        modifiedRight = normalize(stretchProjRight) * (dot(stretchDir, camRight) > 0.0 ? stretchFactor : squashFactor);
        modifiedUp = normalize(stretchProjUp) * (dot(stretchDir, camUp) > 0.0 ? stretchFactor : squashFactor);
        
        // If projection is too small, apply uniform squash instead
        if (length(stretchProjRight) < 0.1) {
            modifiedRight = camRight * squashFactor;
        }
        if (length(stretchProjUp) < 0.1) {
            modifiedUp = camUp * squashFactor;
        }
    }
    
    // Add ripple effect based on time and velocity
    float waveOffset = 0.0;
    
    // Calculate distance from center for ripple effect
    float distFromCenter = length(aPos.xy);
    
    if (velocityMag < 0.5) { 
        // More complex ripple effect for stationary or slow particles
        float primaryWave = sin(distFromCenter * 10.0 + time * rippleSpeed) * 
                           max(0.0, 1.0 - distFromCenter * 2.0);
        
        float secondaryWave = cos(distFromCenter * 15.0 - time * rippleSpeed * 0.7) * 
                             max(0.0, 1.0 - distFromCenter * 2.5) * 0.3;
        
        waveOffset = (primaryWave + secondaryWave) * rippleStrength * (1.0 - velocityMag * 2.0);
    } else {
        // For fast-moving particles, add trailing wake effect
        // More pronounced at the trailing edge
        float angle = atan(aPos.y, aPos.x);
        float alignmentWithVelocity = dot(normalize(vec3(aPos.xy, 0.0)), normalize(vec3(-velocity.xy, 0.0)));
        
        if (alignmentWithVelocity > 0.0) {
            float wakeIntensity = alignmentWithVelocity * distFromCenter;
            waveOffset = sin(distFromCenter * 20.0 + time * rippleSpeed * 3.0) * 
                        wakeIntensity * rippleStrength * 0.5;
        }
    }
    
    // Add subtle random displacement for more natural water movement
    float randomOffset = hash(vec2(particlePos.xy) + time * 0.1) * 0.01;
    waveOffset += randomOffset;
    
    // Generate billboarded vertex position with ripple effect and deformation
    vec3 vertPos = particlePos + 
                  modifiedRight * aPos.x * particleSize + 
                  modifiedUp * aPos.y * particleSize + 
                  billboardNormal * waveOffset;
    
    // Add subtle vertex displacement based on velocity for turbulence
    if (velocityMag > 0.5) {
        float turbulence = sin(distFromCenter * 30.0 + time * 10.0) * cos(distFromCenter * 25.0 - time * 8.0);
        float turbAmt = 0.005 * min(velocityMag, 2.0);
        vertPos += (camRight * sin(time * 15.0) + camUp * cos(time * 12.0)) * turbulence * turbAmt;
    }
    
    // Tangent and bitangent for normal mapping
    // Update these to account for the deformation and make them perpendicular
    Tangent = normalize(modifiedRight);
    Bitangent = normalize(modifiedUp);
    
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