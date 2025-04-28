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

// Generate procedural water texture
vec3 waterTexture(vec2 uv, float time) {
    // Layer several noise patterns for a complex water surface
    float scale1 = 8.0, scale2 = 15.0, scale3 = 25.0;
    
    // Animate UV coordinates based on flow direction and time
    vec2 flowDir = normalize(vec2(ParticleVelocity.xy) + vec2(0.01, 0.01));
    float flowSpeed = length(ParticleVelocity) * 0.5;
    vec2 animatedUV = uv + flowDir * flowSpeed * time * 0.1;
    
    // Generate multiple noise layers
    float noise1 = sin(animatedUV.x * scale1 + time * 0.1) * cos(animatedUV.y * scale1 + time * 0.15) * 0.5;
    float noise2 = sin(animatedUV.x * scale2 + time * 0.2 + 0.3) * cos(animatedUV.y * scale2 + time * 0.1 + 0.7) * 0.3;
    float noise3 = sin(animatedUV.x * scale3 + time * 0.15 + 1.2) * cos(animatedUV.y * scale3 + time * 0.2 + 2.3) * 0.2;
    
    // Combine for final texture
    float combinedNoise = noise1 + noise2 + noise3;
    
    // Create water patterns
    float waterPattern = abs(combinedNoise) * 0.3 + 0.7; // Range 0.7-1.0 for subtle effect
    
    // Mix with base water color
    return waterColor * waterPattern;
}

// Enhanced specular function with anisotropic highlights
float anisotropicSpecular(vec3 normal, vec3 viewDir, vec3 lightDir, vec3 flowDir, float roughness) {
    // Traditional Blinn-Phong
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float NdotH = max(dot(normal, halfwayDir), 0.0);
    float spec = pow(NdotH, specularPower * 2.0);
    
    // Anisotropic component - elongated highlights in flow direction
    vec3 anisotropicDir = normalize(flowDir);
    float anisotropy = length(ParticleVelocity) * 2.0; // More anisotropic when moving faster
    
    // Anisotropic highlight calculation
    vec3 anisotropicHalfway = normalize(halfwayDir - anisotropicDir * max(0.0, dot(halfwayDir, anisotropicDir)));
    float anisotropicDot = max(0.0, dot(normal, anisotropicHalfway));
    float anisotropicSpec = pow(anisotropicDot, specularPower * 0.5);
    
    // Blend between regular and anisotropic specular based on velocity
    return mix(spec, anisotropicSpec, min(1.0, anisotropy));
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
    
    // Use multiple texture coordinates at different scales for more detail
    vec2 animatedTexCoords1 = TexCoord * 3.0 + time * flowSpeed;
    vec2 animatedTexCoords2 = TexCoord * 6.0 - time * flowSpeed * 0.7;
    
    // Get normal from normal map (primary layer)
    vec3 normalMapValue1 = texture(normalMap, animatedTexCoords1).rgb;
    normalMapValue1 = normalMapValue1 * 2.0 - 1.0; // Convert from [0,1] to [-1,1]
    
    // Get secondary normal layer
    vec3 normalMapValue2 = texture(normalMap, animatedTexCoords2).rgb;
    normalMapValue2 = normalMapValue2 * 2.0 - 1.0;
    
    // Combine normal maps with different weights
    vec3 normalMapValue = normalize(normalMapValue1 * 0.7 + normalMapValue2 * 0.3);
    
    // Reduce normal map intensity at the edges for smoother transitions
    float edgeFactor = smoothstep(0.5, 0.35, distFromCenter);
    normalMapValue.xy *= mix(0.05, rippleStrength * 1.5, edgeFactor); // Increased normal strength for more detail
    
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
    
    // Calculate enhanced Fresnel effect (stronger reflections at glancing angles)
    float fresnel = fresnelBias + fresnelScale * pow(1.0 - max(0.0, dot(viewDir, finalNormal)), fresnelPower * 1.2);
    
    // Add secondary fresnel rim highlight for extra shininess
    float rimLight = pow(1.0 - max(0.0, dot(viewDir, finalNormal)), 8.0) * 0.5;
    fresnel += rimLight;
    
    // Sample environment map for reflections
    vec3 reflectDir = reflect(-viewDir, finalNormal);
    
    // Convert 3D reflection direction to 2D texture coordinates
    float m = 2.0 * sqrt(reflectDir.x*reflectDir.x + reflectDir.y*reflectDir.y + (reflectDir.z+1.0)*(reflectDir.z+1.0));
    vec2 envMapCoord = vec2(reflectDir.x/m + 0.5, reflectDir.y/m + 0.5);
    
    // Sample environment texture with slight distortion based on water movement
    vec2 distortedCoord = envMapCoord + normalMapValue.xy * 0.03;
    vec3 reflectionColor = texture(environmentMap, distortedCoord).rgb;
    
    // Generate procedural water texture
    vec3 waterTextureColor = waterTexture(TexCoord * 2.0, time);
    
    // Ambient component
    float ambientStrength = 0.2;
    vec3 ambient = ambientStrength * lightColor;
    
    // Diffuse component
    float diff = max(dot(finalNormal, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;
    
    // Enhanced specular component with anisotropic highlights
    float spec = anisotropicSpecular(finalNormal, viewDir, lightDir, vec3(flowDir, 0.0), 0.2);
    vec3 specular = specularStrength * 1.5 * spec * lightColor; // Increased specular strength
    
    // Add secondary specular highlight for extra shininess
    float secondarySpec = pow(max(dot(finalNormal, normalize(lightDir + vec3(0.1, 0.2, 0.1))), 0.0), specularPower * 0.5);
    specular += secondarySpec * lightColor * 0.3;
    
    // Enhance specular at particle-particle boundaries for cohesive look
    if (blendFactor > 0.0) {
        specular = mix(specular, specular * 1.8, blendFactor);
    }
    
    // Apply lighting to water texture
    vec3 litColor = (ambient + diffuse) * waterTextureColor + specular;
    
    // Mix water color and reflection based on Fresnel term and setting
    vec3 finalColor = litColor;
    if (reflectionsEnabled) {
        finalColor = mix(litColor, reflectionColor, fresnel * reflectionStrength * 1.3); // Increased reflection
    }
    
    // Add extra highlights at the edges for depth
    finalColor += rimLight * lightColor * 0.4;
    
    // Add velocity-based highlight (brighter in direction of motion)
    if (length(ParticleVelocity) > 0.1) {
        vec3 velDir = normalize(vec3(ParticleVelocity));
        float velHighlight = pow(max(0.0, dot(viewDir, velDir)), 8.0) * length(ParticleVelocity);
        finalColor += velHighlight * vec3(0.8, 0.9, 1.0) * 0.5;
    }
    
    // Add subtle blue/cyan tint to enhance water feel
    finalColor += vec3(0.0, 0.05, 0.1) * circularMask;
    
    // Add subtle iridescence effect
    float iridescence = sin(dot(viewDir, finalNormal) * 10.0 + time * 0.5) * 0.5 + 0.5;
    finalColor += vec3(0.0, iridescence * 0.05, iridescence * 0.1) * spec;
    
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