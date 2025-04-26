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
    float specularStrength = 0.8;
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(norm, halfwayDir), 0.0), 64.0);
    vec3 specular = specularStrength * spec * lightColor;
    
    // Apply water-like effect (refraction and color depth)
    float depthFactor = mix(0.7, 1.0, 1.0 - dist * 1.5);
    vec3 waterColor = mix(objectColor, objectColor * 1.3, 1.0 - dist);
    
    vec3 result = (ambient + diffuse + specular) * waterColor * depthFactor;
    
    // Apply transparency gradient from center to edge
    float alpha = mask * mix(1.0, 0.7, dist * 2.0);
    
    FragColor = vec4(result, alpha);
}