#version 330 core
out vec4 FragColor;

uniform vec3 ObjectColor;

void main()
{
	FragColor = vec4(ObjectColor, 1.0);
}