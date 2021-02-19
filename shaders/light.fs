#version 330 core

#define POINT_LIGHT_COUNT 5

out vec4 FragColor;

struct material
{
    sampler2D Diffuse;
    sampler2D Specular;
    float Shininess;
};

struct dirLight
{
    vec3 Direction;

    vec3 Ambient;
    vec3 Diffuse;
    vec3 Specular;
};

struct pointLight
{
    vec3 Position;

    float Constant;
    float Linear;
    float Quadratic;

    vec3 Ambient;
    vec3 Diffuse;
    vec3 Specular;
};

in vec3 Normal;
in vec3 FragPos;
in vec2 TexCoord;

uniform vec3 ViewPos;
uniform dirLight DirLight;
uniform pointLight PointLights[POINT_LIGHT_COUNT];
uniform material Material;
uniform float TexScale;

// NOTE(Jovan): Scaled TexCoord
vec2 TexCoords;

vec3
CalcDirLight(dirLight light, vec3 normal, vec3 viewDir);
vec3
CalcPointLight(pointLight, vec3 normal, vec3 fragPos, vec3 viewDir);

void
main()
{
    TexCoords = TexScale * (TexCoord + 0.5) + 0.5 * TexScale; 

    // TODO(Jovan): Optimize, remove duplicate code (reflection, etc...)
    vec3 Norm = normalize(Normal);
    vec3 ViewDir = normalize(ViewPos - FragPos);

    //NOTE(Jovan): Directional lighting
    vec3 Result = CalcDirLight(DirLight, Norm, ViewDir);

    // NOTE(Jovan): Point lighting
    for(int PtLightIndex = 0;
        PtLightIndex < POINT_LIGHT_COUNT;
        ++PtLightIndex)
    {
        Result += CalcPointLight(PointLights[PtLightIndex], Norm, FragPos, ViewDir);
    }

    // TODO(Jovan): Spotlight

    FragColor = vec4(Result, 1.0);
}

vec3
CalcDirLight(dirLight light, vec3 normal, vec3 viewDir)
{
    vec3 LightDir = normalize(-light.Direction);
    // NOTE(Jovan): Diffuse
    float Diff = max(dot(normal, LightDir), 0.0);
    // NOTE(Jovan): Specular
    vec3 ReflectDir = reflect(-LightDir, normal);
    float Spec = pow(max(dot(viewDir, ReflectDir), 0.0), Material.Shininess);
    // NOTE(Jovan): Ambient and phong
    vec3 Ambient  = light.Ambient  * vec3(texture(Material.Diffuse, TexCoords));
    vec3 Diffuse  = light.Diffuse  * Diff * vec3(texture(Material.Diffuse, TexCoords));
    vec3 Specular = light.Specular * Spec * vec3(texture(Material.Specular, TexCoords));
    return (Ambient + Diffuse + Specular);
}

vec3
CalcPointLight(pointLight light, vec3 normal, vec3 fragPos, vec3 viewDir)
{
    vec3 LightDir = normalize(light.Position - fragPos);
    // NOTE(Jovan): Diffuse
    float Diff = max(dot(normal, LightDir), 0.0);
    // NOTE(Jovan): Specular
    vec3 ReflectDir = reflect(-LightDir, normal);
    float Spec = pow(max(dot(viewDir, ReflectDir), 0.0), Material.Shininess);
    // NOTE(Jovan): Attenuation
    float Distance    = length(light.Position - fragPos);
    float Attenuation = 1.0 / (light.Constant + light.Linear * Distance + 
        light.Quadratic * (Distance * Distance));    
    // NOTE(Jovan): Ambient and phong
    vec3 Ambient  = light.Ambient  * vec3(texture(Material.Diffuse, TexCoords));
    vec3 Diffuse  = light.Diffuse  * Diff * vec3(texture(Material.Diffuse, TexCoords));
    vec3 Specular = light.Specular * Spec * vec3(texture(Material.Specular, TexCoords));
    Ambient  *= Attenuation;
    Diffuse  *= Attenuation;
    Specular *= Attenuation;
    return (Ambient + Diffuse + Specular);
}