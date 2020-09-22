/*
 * FluidRenderer.cpp
 *
 *  Created on: 23/07/2011
 *      Author: josericardo
 */

#include "FluidRenderer.h"

namespace Simulation
{
	namespace Fluid
	{
		void FluidRenderer::initializeVBO(float4* _positions, int count)
		{
			positions.genBuffer();
			positions.bindBuffer();
			positions.setData(sizeof(float4) * count, _positions, GL_DYNAMIC_DRAW);
			positions.unbindBuffer();
            
            // Criar texturas
            m_ColorTex = createTexture(GL_TEXTURE_2D, 800, 600, GL_RGBA8, GL_RGBA);
            m_DepthTex = createTexture(GL_TEXTURE_2D, 800, 600, GL_DEPTH_COMPONENT24, GL_DEPTH_COMPONENT);
            
            // Criar FBO
            m_FBO.AttachTexture(GL_TEXTURE_2D, m_ColorTex, GL_COLOR_ATTACHMENT0);
            m_FBO.AttachTexture(GL_TEXTURE_2D, m_DepthTex, GL_DEPTH_ATTACHMENT);
            m_FBO.IsValid();
            
		}
        
        GLuint FluidRenderer::createTexture(GLenum target, int w, int h, GLint internalFormat, GLenum format)
        {
            GLuint texId;
            glGenTextures(1, &texId);
            glBindTexture(target, texId);
            glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexImage2D(target, 0, internalFormat, w, h, 0, format, GL_UNSIGNED_BYTE, NULL);
            return texId;
        }
	}
}
