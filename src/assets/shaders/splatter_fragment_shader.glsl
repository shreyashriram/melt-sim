#version 330 core
out vec4 FragColor;

in vec3 FragPos;
in vec3 Normal;
in vec2 TexCoord;
in vec3 Tangent;
in vec3 Bitangent;
in vec3 ParticleVelocity;
in vec4 ParticleParams;

// Lighting
uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;
uniform vec3 waterColor;

// Water parameters
uniform float reflectionStrength;
uniform float specularPower;
uniform float specularStrength;
uniform float fresnelBias;
uniform float fresnelScale;
uniform float fresnelPower;
uniform float rippleStrength;
uniform float rippleSpeed;
uniform float time;
uniform float surfaceTensionStrength;
uniform float blendDistance;

// Blending parameters
uniform int numParticles;
uniform vec4 particlePositions[100];

// Feature toggles
uniform bool reflectionsEnabled;
uniform bool refractionsEnabled;

// Textures
uniform sampler2D normalMap;
uniform sampler2D environmentMap;

// Calculate how much particles influence each other at a point
float calculateBlendFactor(vec3 position) {
    float blendFactor = 0.0;
    
    // Only calculate blend if we're not too far from any particle center
    vec2 texCoordFromCenter = TexCoord - vec2(0.5);
    float distFromCenter = length(texCoordFromCenter);
    
    // Only do the calculation near the edges of particles
    if (distFromCenter > 0.3) {
        for (int i = 0; i < numParticles; i++) {
            vec3 otherPos = particlePositions[i].xyz;
            float otherSize = particlePositions[i].w;
            
            // Skip if it's the same particle
            vec3 diff = otherPos - ParticleParams.xyz;
            float distSquared = dot(diff, diff);
            if (distSquared < 0.0001) continue;
            
            // Calculate how much this other particle influences this point
            // Use an exponential falloff
            float blendThreshold = blendDistance * blendDistance;
            if (distSquared < blendThreshold) {
                // Calculate influence based on distance
                float influence = exp(-distSquared / (blendThreshold * 0.5));
                
                // Scale influence based on how close we are to the edge
                influence *= smoothstep(0.0, 0.4, distFromCenter);
                
                // Add to total blend factor
                blendFactor += influence * surfaceTensionStrength;
            }
        }
    }
    
    // Clamp to reasonable range
    return min(blendFactor, 1.0);
}

void main() {
    // Calculate base particle shape (circular falloff)
    vec2 texCoordFromCenter = TexCoord - vec2(0.5);
    float distFromCenter = length(texCoordFromCenter);
    float circularMask = smoothstep(0.5, 0.4, distFromCenter);
    
    // Discard fragments outside the circular shape
    if (circularMask < 0.01)
        discard;
    
    // Animate texture coordinates for flowing water effect
    vec2 flowSpeed = vec2(0.03, 0.02) * rippleSpeed;
    vec2 flowDir = normalize(vec2(ParticleVelocity.xy) + vec2(0.01, 0.01));
    if (length(ParticleVelocity) > 0.1) {
        flowSpeed = flowSpeed + flowDir * length(ParticleVelocity) * 0.05;
    }
    vec2 animatedTexCoords = TexCoord * 3.0 + time * flowSpeed;
    
    // Get normal from normal map
    vec3 normalMapValue = texture(normalMap, animatedTexCoords).rgb;
    normalMapValue = normalMapValue * 2.0 - 1.0; // Convert from [0,1] to [-1,1]
    
    // Reduce normal map intensity at the edges for smoother transitions
    float edgeFactor = smoothstep(0.5, 0.35, distFromCenter);
    normalMapValue.xy *= mix(0.05, rippleStrength, edgeFactor);
    
    // Calculate blend factor between particles
    float blendFactor = calculateBlendFactor(FragPos);
    
    // Construct TBN matrix for normal mapping
    vec3 N = normalize(Normal);
    vec3 T = normalize(Tangent);
    vec3 B = normalize(Bitangent);
    mat3 TBN = mat3(T, B, N);
    
    // Apply normal map with TBN matrix
    vec3 mappedNormal = TBN * normalize(normalMapValue);
    
    // Smoothly blend between mapped normal and original normal near edges
    float normalBlendFactor = smoothstep(0.45, 0.3, distFromCenter);
    vec3 finalNormal = normalize(mix(N, mappedNormal, normalBlendFactor));
    
    // Apply particle blending to normal
    if (blendFactor > 0.0) {
        // Flatten normal at particle intersections for smoother transitions
        finalNormal = normalize(mix(finalNormal, N, blendFactor * 0.7));
    }
    
    // Calculate basic lighting vectors
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 lightDir = normalize(lightPos - FragPos);
    vec3 halfwayDir = normalize(lightDir + viewDir);
    
    // Calculate Fresnel effect (stronger reflections at glancing angles)
    float fresnel = fresnelBias + fresnelScale * pow(1.0 - max(0.0, dot(viewDir, finalNormal)), fresnelPower);
    
    // Sample environment map for reflections
    // We use a simple spherical mapping for the environment 
    vec3 reflectDir = reflect(-viewDir, finalNormal);
    
    // Convert 3D reflection direction to 2D texture coordinates
    // This is a simple spherical mapping - a real implementation would use cubemaps
    float m = 2.0 * sqrt(reflectDir.x*reflectDir.x + reflectDir.y*reflectDir.y + (reflectDir.z+1.0)*(reflectDir.z+1.0));
    vec2 envMapCoord = vec2(reflectDir.x/m + 0.5, reflectDir.y/m + 0.5);
    
    // Sample environment texture
    vec3 reflectionColor = texture(environmentMap, envMapCoord).rgb;
    
    // Ambient component
    float ambientStrength = 0.2;
    vec3 ambient = ambientStrength * lightColor;
    
    // Diffuse component
    float diff = max(dot(finalNormal, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;
    
    // Specular component (with blending)
    float spec = pow(max(dot(finalNormal, halfwayDir), 0.0), specularPower);
    
    // Enhance specular at particle-particle boundaries for cohesive look
    float enhancedSpec = spec;
    if (blendFactor > 0.0) {
        enhancedSpec = mix(spec, spec * 1.5, blendFactor);
    }
    
    vec3 specular = specularStrength * enhancedSpec * lightColor;
    
    // Apply water color with lighting
    vec3 baseColor = waterColor;
    vec3 litColor = (ambient + diffuse) * baseColor + specular;
    
    // Mix water color and reflection based on Fresnel term and setting
    vec3 finalColor = litColor;
    if (reflectionsEnabled) {
        finalColor = mix(litColor, reflectionColor, fresnel * reflectionStrength);
    }
    
    // Add extra highlights at the edges for depth
    float rimLight = pow(1.0 - max(0.0, dot(viewDir, finalNormal)), 4.0);
    finalColor += rimLight * lightColor * 0.2;
    
    // Add velocity-based highlight (brighter in direction of motion)
    if (length(ParticleVelocity) > 0.1) {
        vec3 velDir = normalize(vec3(ParticleVelocity));
        float velHighlight = pow(max(0.0, dot(viewDir, velDir)), 8.0) * length(ParticleVelocity);
        finalColor += velHighlight * vec3(0.7, 0.8, 1.0) * 0.3;
    }
    
    // Adjust alpha for smoother blending between particles
    float alpha = circularMask;
    
    // Make edges more transparent
    float edgeTransparency = smoothstep(0.5, 0.4, distFromCenter);
    alpha *= edgeTransparency;
    
    // Strengthen alpha in areas with particle overlap
    alpha = mix(alpha, min(alpha + 0.2, 1.0), blendFactor);
    
    // Output final color with alpha
    FragColor = vec4(finalColor, alpha);
}