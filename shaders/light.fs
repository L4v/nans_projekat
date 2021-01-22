#version 330 core
out vec4 FragColor;

struct material
{
    sampler2D Diffuse;
    sampler2D Specular;
    float Shininess;
};

struct light
{
    vec3 Position;

    vec3 Ambient;
    vec3 Diffuse;
    vec3 Specular;
};

in vec3 Normal;
in vec3 FragPos;
in vec2 TexCoord;

uniform vec3 ViewPos;
uniform light Light;
uniform material Material;

void main()
{
    // NOTE(Jovan): Ambient
    vec3 Ambient = Light.Ambient * vec3(texture(Material.Diffuse, TexCoord));

    // NOTE(Jovan): Diffuse
    vec3 Norm = normalize(Normal);
    vec3 LightDir = normalize(Light.Position - FragPos);
    float Diff = max(dot(Norm, LightDir), 0.0);
    vec3 Diffuse = Light.Diffuse * Diff * vec3(texture(Material.Diffuse, TexCoord));

    // NOTE(Jovan): Specular
    vec3 ViewDir = normalize(ViewPos - FragPos);
    vec3 ReflectDir = reflect(-LightDir, Norm);
    float Spec = pow(max(dot(ViewDir, ReflectDir), 0.0), Material.Shininess);
    vec3 Specular = Light.Specular * Spec * vec3(texture(Material.Specular, TexCoord));

    // NOTE(Jovan): Phong
    FragColor = vec4(Ambient + Diffuse + Specular, 1.0);
}