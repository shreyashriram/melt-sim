#version 330 core
        out vec4 FragColor;

        in vec3 FragPos;
        in vec3 Normal;
        in vec2 TexCoord;
        in float DistFromCenter;
        in vec4 ParticleParams; // xyz = position, w = size
        in vec3 ParticleVelocity; // Velocity for trail effects

        // Metaball parameters
        uniform int numParticles;         // Total number of particles
        uniform vec4 particlePositions[100]; // Array of particle positions (increase if needed)
        uniform float metaballThreshold = 1.0;  // Threshold for metaball effect
        uniform float metaballStrength = 0.5;   // Strength of metaball effect

        // Standard lighting parameters
        uniform vec3 lightPos;
        uniform vec3 viewPos;
        uniform vec3 lightColor;
        uniform vec3 objectColor;
        uniform float smoothing = 0.8;

        // Water droplet parameters
        uniform bool enableWaterDroplets = true;
        uniform float dropletScale = 15.0;     // Scale of water droplet pattern
        uniform float dropletIntensity = 0.5;  // Intensity of the water droplet effect
        uniform sampler2D waterTexture;        // Optional texture for water droplets

        // Random functions for noise
        float random(vec2 st) {
            return fract(sin(dot(st.xy, vec2(12.9898, 78.233))) * 43758.5453123);
        }

        float noise(vec2 st) {
            vec2 i = floor(st);
            vec2 f = fract(st);
            
            // Four corners in 2D of a tile
            float a = random(i);
            float b = random(i + vec2(1.0, 0.0));
            float c = random(i + vec2(0.0, 1.0));
            float d = random(i + vec2(1.0, 1.0));

            // Smooth interpolation
            vec2 u = f * f * (3.0 - 2.0 * f);

            return mix(a, b, u.x) + 
                (c - a) * u.y * (1.0 - u.x) + 
                (d - b) * u.x * u.y;
        }

        // Function to create water droplet pattern
        float waterDropletPattern(vec2 uv, float scale, float time) {
            float n = noise(uv * scale);
            
            // Create water droplet effect by adding multiple noise layers
            float droplet = smoothstep(0.4, 0.5, n);
            
            // Add some variation for more natural look
            droplet += 0.1 * noise(uv * scale * 2.0 + vec2(time * 0.1));
            
            return droplet;
        }

        void main() {
            // Calculate velocity effect
            float velocityMag = length(ParticleVelocity);
            float velocityEffect = clamp(velocityMag * 0.5, 0.0, 1.0);
            
            // Calculate distance from center of splat, adjusted for velocity
            float dist = length(TexCoord - vec2(0.5));
            
            // Adjust edge for velocity - more elongated droplets have sharper edges
            float edge = mix(0.45, 0.3, velocityEffect);
            float mask = smoothstep(0.5, edge, dist);
            
            // Trail effect for faster particles
            float trailEffect = smoothstep(0.0, 1.0, velocityMag * 0.3);
            
            // Apply trail to mask - make particles with higher velocity more transparent at the edges
            mask = mix(mask, mask * (1.0 - dist), trailEffect * 0.5);
            
            // Discard fragments outside the particle
            if (mask < 0.05)
                discard;
            
            // Standard lighting calculation
            vec3 norm = normalize(Normal);
            vec3 lightDir = normalize(lightPos - FragPos);
            vec3 viewDir = normalize(viewPos - FragPos);
            
            float ambientStrength = 0.4;
            vec3 ambient = ambientStrength * lightColor;
            
            float diff = max(dot(norm, lightDir), 0.0) * 0.4;
            vec3 diffuse = diff * lightColor;
            
            float specularStrength = 0.1;
            vec3 halfwayDir = normalize(lightDir + viewDir);
            float spec = pow(max(dot(norm, halfwayDir), 0.0), 16.0);
            vec3 specular = specularStrength * spec * lightColor;
            
            // Generate water droplet texture
            float droplet = 0.0;
            if (enableWaterDroplets) {
                // Use particle position as a unique seed for variation
                vec2 seed = vec2(ParticleParams.x * 10.0, ParticleParams.y * 10.0);
                
                // Create water droplet pattern
                droplet = waterDropletPattern(TexCoord * 3.0 + seed, dropletScale, 0.0);
                droplet *= dropletIntensity;
            }
            
            // Enhanced droplet effects based on velocity
            if (velocityMag > 0.2) {
                // Add streaks in the direction of movement for fast particles
                vec2 normalizedVelocity2D = normalize(vec2(ParticleVelocity.x, ParticleVelocity.y));
                float streakEffect = pow(abs(dot(normalize(TexCoord - vec2(0.5)), normalizedVelocity2D)), 3.0);
                droplet += streakEffect * velocityEffect * 0.3;
            }
            
            // Apply water droplet effect to normal and color
            vec3 adjustedNormal = normalize(norm + vec3(droplet * 0.1, droplet * 0.1, 0.0));
            
            // Metaball calculation - influences the alpha
            float metaballField = 0.0;
            // Self contribution
            metaballField += 1.0 / (dist * dist + 0.01) * metaballStrength;
            
            // Apply the base color with water droplet effect
            float dropletHighlight = droplet * 0.2;
            vec3 fluidColor = objectColor * (1.0 + (0.05 * (1.0 - dist) + dropletHighlight));
            
            // Enhanced fluid appearance with varying depth and color
            float depthFactor = 0.95 - droplet * 0.1;
            
            // Add velocity-based color tinting (blue to white as velocity increases)
            vec3 velocityColor = mix(fluidColor, vec3(1.0, 1.0, 1.0), velocityEffect * 0.3);
            
            // Combine lighting
            vec3 result = (ambient + diffuse + specular) * velocityColor * depthFactor;
            
            // Enhanced alpha for metaball effect
            float alphaFactor = smoothstep(0.0, 0.25, mask);
            float baseAlpha = mix(0.9, 0.7, dist * 0.3) * alphaFactor;
            
            // Add subtle refraction at the edges
            float refractionFactor = 0.05 * (1.0 - mask) * (1.0 + droplet);
            result += refractionFactor * vec3(0.8, 0.9, 1.0);
            
            // Add subtle water caustics effect
            float causticEffect = noise(TexCoord * 20.0) * noise(TexCoord * 15.0 + vec2(0.2, 0.3));
            result += causticEffect * 0.1 * baseAlpha * vec3(0.8, 0.95, 1.0);
            
            // Final color with trail transparency for fast-moving particles
            float finalAlpha = baseAlpha * (1.0 - trailEffect * 0.2);
            FragColor = vec4(result, finalAlpha);
        }