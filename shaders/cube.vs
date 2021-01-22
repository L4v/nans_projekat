#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec2 aTexCoord;
layout (location = 2) in vec3 aNormal;

uniform mat4 Model;
uniform mat4 View;
uniform mat4 Projection;

out vec2 TexCoord;
out vec3 FragPos;
out vec3 Normal;

void main()
{
	gl_Position = Projection * View * Model * vec4(aPos, 1.0);
	TexCoord = aTexCoord;
	Normal = aNormal;
	FragPos = vec3(Model * vec4(aPos, 1.0));
}