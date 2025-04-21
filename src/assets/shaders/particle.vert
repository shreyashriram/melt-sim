#version 330 core
layout (location = 0) in vec2 aPos;
uniform mat4 projection;
uniform mat4 view;

void main() {
    gl_PointSize = 5.0;
    gl_Position = projection * view * vec4(aPos, 0.0, 1.0);
}