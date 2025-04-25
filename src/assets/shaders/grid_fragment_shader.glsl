
// grid_fragment_shader.glsl
#version 330 core
out vec4 FragColor;

in float mass; // Receiving mass from vertex shader
in float arrowFlag; // Flag indicating if this is an arrow

uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;
uniform vec3 arrowColor = vec3(1.0, 1.0, 1.0); 

void main()
{
    // Check if we're rendering an arrow or a grid node
    if (arrowFlag > 0.5) {
        // This is an arrow - use the arrow color
        FragColor = vec4(arrowColor, 1.0);
    } else {
        // This is a grid node - use the mass-based coloring
        vec3 color;
        // Normalize mass to [0,1] range using a max mass reference
        float maxMass = 10.0; // Adjust based on expected mass range
        float normalizedMass = clamp(mass / maxMass, 0.0, 1.0);
        
        // Blue to Red gradient through green
        if (normalizedMass < 0.5) {
            // Blue to Green (low to medium mass)
            float t = normalizedMass * 2.0;
            color = vec3(0.0, t, 1.0 - t);
        } else {
            // Green to Red (medium to high mass)
            float t = (normalizedMass - 0.5) * 2.0;
            color = vec3(t, 1.0 - t, 0.0);
        }
        
        // Apply basic lighting
        vec3 ambient = 0.3 * color;
        // Final color with ambient lighting
        FragColor = vec4(ambient + 0.7 * color, 1.0);
        
        // Add a little glow effect for particles
        float distFromCenter = length(gl_PointCoord - vec2(0.5));
        if (distFromCenter > 0.5) {
            discard; // Create circular points
        }
        
        // Add glow based on mass
        float glow = smoothstep(0.5, 0.0, distFromCenter);
        FragColor.rgb += 0.3 * glow * vec3(normalizedMass, normalizedMass * 0.5, 0.0);
    }
}