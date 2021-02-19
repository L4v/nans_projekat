#version 330 core

out vec4 FragColor;

in vec4 vertexColor;

uniform sampler2D Texture;

void main()
{
	FragColor = vertexColor;
}