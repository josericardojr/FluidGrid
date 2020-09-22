#ifndef __Emitter_H__

#define __Emitter_H__

#include <vector>
#include "GPUSPHFluid.h"
#include "JRFXGL.h"

#ifdef __linux__
	#include <stdint.h>
	typedef uint32_t uint32;
#else
	//typedef uint32 uint32;
#endif

/* 
 * O tempo de vida da partícula é armazenado no vetor posição, no campo w (lifetime < -100 = infinito)
*/

namespace Simulation
{	
	namespace Fluid { class GPUSPHFluid; }

	namespace Emitter
	{
		class IEmitter
		{
		public:
			virtual void update(Simulation::Fluid::GPUSPHFluid* fluid) = 0;

			virtual void reset(float4* positions, int count) = 0;
					
		};


		class BlockEmitter : public IEmitter
		{
		public:
			/// <summary>
			/// Initializes a new instance of the BlockEmitter class.
			/// </summary>
			BlockEmitter(mc::math::Vec3f _block_size_min, mc::math::Vec3f _block_size_max, float _particle_spacing);

			void reset(float4* positions, int count);

			void update(Simulation::Fluid::GPUSPHFluid* fluid){}

			void UpdateData(){}

		private:
			mc::math::Vec3f  mBlockSizeMin;
			mc::math::Vec3f  mBlockSizeMax;
			float            mParticleSpacing;

		};
	}
}
#endif
