#version 330
uniform sampler2D Image;

out vec4 FragColor;
in vec2 texCoords;
void main()
{
	FragColor = vec4(texture(Image,texCoords).xyz,1.0);
}
