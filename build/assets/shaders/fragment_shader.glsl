#version 330 core
out vec4 FragColor;

in vec3 FragPos;    // from vertex shader
in vec3 Normal;     // from vertex shader

uniform vec3 lightPos;   // light position in world space
uniform vec3 viewPos;    // camera position
uniform vec3 lightColor; // usually white
uniform vec3 objectColor; // water base color, like vec3(0.0, 0.4, 0.7)

void main()
{
    // Normalize inputs
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
    float spec = pow(max(dot(norm, halfwayDir), 0.0), 64.0); // Shininess
    vec3 specular = specularStrength * spec * lightColor;

    vec3 result = (ambient + diffuse + specular) * objectColor;

    FragColor = vec4(result, 1.0);
}