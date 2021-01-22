#version 330 core
out vec4 FragColor;

in vec3 Normal;
in vec3 FragPos;

uniform vec3 LightPos;
uniform vec3 LightColor;
uniform vec3 ObjectColor;

void main()
{
    float AmbientStrength = 0.1;
    vec3 Ambient = AmbientStrength * LightColor;

    vec3 Norm = normalize(Normal);
    vec3 LightDir = normalize(LightPos - FragPos);
    float Diff = max(dot(Norm, LightDir), 0.0);
    vec3 Diffuse = Diff * LightColor;

    vec3 Result = Diff * LightColor;
    FragColor = vec4(Result, 1.0);
}