/*
 * FluidRenderer.h
 *
 *  Created on: 22/07/2011
 *      Author: josericardo
 */

#ifndef FLUIDRENDERER_H_
#define FLUIDRENDERER_H_

#include "jrfxgl.h"
#include "Gl/glew.h"
#include <vector_types.h>

namespace Simulation
{
	namespace Fluid
	{
		class FluidRenderer : public JRFXGL::Graphics::Mesh
		{
		public:
			JRFXGL::Graphics::GLVertexBuffer&   getVBOPosition() {return positions; }

		protected:
			void initializeVBO(float4* _positions, int count);
            
        private:
            GLuint createTexture(GLenum target, int w, int h, GLint internalFormat, GLenum format);

		private:
			JRFXGL::Graphics::GLVertexBuffer     positions;
            JRFXGL::Graphics::FramebufferObject  m_FBO;
            GLuint                               m_ColorTex;
            GLuint                               m_DepthTex;
		};
	}
}

#endif /* FLUIDRENDERER_H_ */
