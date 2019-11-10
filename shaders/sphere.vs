#version 330 core
layout (location = 0) in vec3 aPosa;
layout (location = 1) in vec2 aTexCoorda;

uniform mat4 Model;
uniform mat4 View;
uniform mat4 Projection;

out vec2 TexCoorda;

void main()
{
	gl_Position = Projection * View * Model * vec4(aPosa, 1.0);
	TexCoorda = aTexCoorda;
}