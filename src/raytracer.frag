#version 330
uniform sampler2D Image;
uniform int pass;

out vec4 FragColor;
in vec2 texCoords;
void main()
{
	vec3 color = texture(Image,texCoords).xyz;
	FragColor = vec4(color.x/pass,color.y/pass,color.z/pass,1.0);
}
