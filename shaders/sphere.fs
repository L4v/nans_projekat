#version 330 core

out vec4 FragColora;

in vec2 TexCoorda;

uniform sampler2D Texture;

void main()
{
	FragColora = texture(Texture, TexCoorda);
}