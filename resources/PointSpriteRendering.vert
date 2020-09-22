#version 120
									 
uniform float pointRadius;
uniform float pointScale;
uniform float densityScale;
uniform float densityOffset;
varying vec4 eyeSpaceCoord;

void main()
{
    // calculate window-space point size
    vec3 posEye = vec3(gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0));
    float dist = length(posEye);
    gl_PointSize = pointRadius * (pointScale / dist);
	
    gl_TexCoord[0] = gl_MultiTexCoord0;
    gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);
	
    gl_FrontColor = gl_Color;
    eyeSpaceCoord = vec4(gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0));
}