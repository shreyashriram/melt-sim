#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in float aMass; // Mass attribute

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

// Add a uniform to indicate whether we're rendering nodes or arrows
uniform float isArrow = 0.0; // 0.0 for nodes, 1.0 for arrows

out float mass; // Pass mass to fragment shader
out float arrowFlag; // Flag to indicate if this is an arrow

void main()
{
    gl_Position = projection * view * model * vec4(aPos, 1.0);
    mass = aMass; // Pass the mass to the fragment shader
    arrowFlag = isArrow; // Pass the arrow flag to fragment shader
}

