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