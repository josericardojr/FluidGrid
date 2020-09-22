/*
 * FluidParticleRenderMaterial.h
 *
 *  Created on: 25/07/2011
 *      Author: josericardo
 */

#ifndef FLUIDPARTICLERENDERMATERIAL_H_
#define FLUIDPARTICLERENDERMATERIAL_H_

#include "JRFXGL.h"
#include <boost/shared_ptr.hpp>

namespace Simulation
{
	namespace Fluid
	{
		class FluidParticleMaterial : public JRFXGL::Graphics::Material
		{
		public:
			FluidParticleMaterial();

			void setPointRadius(float _pointradius){ pointRadius = _pointradius; }
			void setPointScale(float _pointScale){ pointScale = _pointScale; }

			void apply();
			void disable();

		private:
			JRFXGL::Graphics::ProgramPtr program;
			float pointRadius;
			float pointScale;
		};

		typedef boost::shared_ptr<FluidParticleMaterial> FluidParticleMaterialPtr;

	}
}


#endif /* FLUIDPARTICLERENDERMATERIAL_H_ */
