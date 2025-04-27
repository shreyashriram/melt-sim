#version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aNormal;
        layout (location = 2) in vec4 aParticleData; // xyz = position, w = size
        layout (location = 3) in vec4 aParticleVelocity; // xyz = velocity, w = unused

        out vec3 FragPos;
        out vec3 Normal;
        out vec2 TexCoord;
        out float DistFromCenter;
        out vec4 ParticleParams; // Pass particle data to fragment shader
        out vec3 ParticleVelocity; // Pass velocity for fragment shader effects

        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        uniform float dropletDeformFactor = 0.7; // Controls how much to stretch particles based on velocity

        void main() {
            // Extract particle data
            vec3 particlePos = aParticleData.xyz;
            float particleSize = aParticleData.w;
            vec3 velocity = aParticleVelocity.xyz;
            
            // Calculate velocity magnitude and direction
            float velocityMag = length(velocity);
            vec3 velocityDir = velocityMag > 0.01 ? normalize(velocity) : vec3(0.0, -1.0, 0.0);
            
            // Calculate billboarding vectors - aligned with camera but stretched in velocity direction
            vec3 camRight = normalize(vec3(view[0][0], view[1][0], view[2][0]));
            vec3 camUp = normalize(vec3(view[0][1], view[1][1], view[2][1]));
            
            // Transform these basis vectors to align one with velocity direction for stretching
            // We want to stretch in the velocity direction and compress perpendicular to it
            
            // Calculate the billboard right and up directions
            vec3 billboardNormal = normalize(vec3(view[0][2], view[1][2], view[2][2]));
            
            // Find the best alignment between velocity and the camera basis
            // Project velocity onto camera plane
            vec3 projVelocity = velocity - dot(velocity, billboardNormal) * billboardNormal;
            vec3 projVelocityDir = length(projVelocity) > 0.01 ? normalize(projVelocity) : normalize(camUp);

            // Calculate perpendicular vector to both the projected velocity and billboard normal
            vec3 perpVelocity = normalize(cross(projVelocityDir, billboardNormal));
            
            // Deformation factors - stretch in velocity direction, compress in perpendicular
            float stretchFactor = 1.0 + velocityMag * dropletDeformFactor;
            float compressFactor = 1.0 / (1.0 + velocityMag * 0.3);
            
            // Scale in velocity direction
            vec3 adjustedRight, adjustedUp;
            
            if (velocityMag > 0.01) {
                // Use the projected velocity as one axis and the perpendicular as the other
                adjustedRight = perpVelocity * compressFactor;
                adjustedUp = projVelocityDir * stretchFactor;
            } else {
                // When nearly stationary, use regular billboarding
                adjustedRight = camRight;
                adjustedUp = camUp;
            }
            
            // Generate billboarded vertex position with velocity-based deformation
            vec3 vertPos = particlePos + 
                          adjustedRight * aPos.x * particleSize + 
                          adjustedUp * aPos.y * particleSize;
            
            // Create a normal that gives good lighting with the deformed shape
            vec3 faceNormal;
            if (velocityMag > 0.01) {
                // Blend between velocity direction and camera-facing normal for dynamic billboarding
                faceNormal = normalize(mix(billboardNormal, -velocityDir, 0.3));
            } else {
                // Standard normal when nearly stationary
                faceNormal = billboardNormal;
            }
            
            // We need to convert these to world space for our existing shader format
            FragPos = vec3(model * vec4(vertPos, 1.0));
            Normal = mat3(transpose(inverse(model))) * faceNormal;
            
            // Pass texture coordinates for the splat - adjusted to account for stretching
            TexCoord = aPos.xy + 0.5; // Convert from [-0.5, 0.5] to [0, 1]
            
            // Calculate distance from center for the fragment shader
            // Account for the elongated shape by scaling the distance calculation
            vec2 scaledPos = vec2(aPos.x / compressFactor, aPos.y / stretchFactor);
            DistFromCenter = length(scaledPos);
            
            // Pass particle parameters to fragment shader for unique water effects
            ParticleParams = aParticleData;
            ParticleVelocity = velocity;
            
            gl_Position = projection * view * model * vec4(vertPos, 1.0);
        }