/*
 * FluidParticleRenderMaterial.cpp
 *
 *  Created on: 25/07/2011
 *      Author: josericardo
 */

#include "FluidParticleRenderMaterial.h"


namespace Simulation
{
	namespace Fluid
	{
		FluidParticleMaterial::FluidParticleMaterial()
		{
			program.reset(new JRFXGL::Graphics::Program(
					"resources/PointSpriteRendering.vert",
					"resources/PointSpriteRendering.frag"));

		}

		void FluidParticleMaterial::apply()
		{
			program->use();
			program->setUniform("pointRadius", pointRadius);
			program->setUniform("pointScale", pointScale);
		}

		void FluidParticleMaterial::disable()
		{
			program->disable();
		}
	}
}
