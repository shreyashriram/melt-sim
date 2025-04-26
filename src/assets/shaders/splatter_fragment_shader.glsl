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