#include "Emitter.h"

namespace Simulation
{
	namespace Emitter
	{

		BlockEmitter::BlockEmitter(mc::math::Vec3f _block_size_min,
				mc::math::Vec3f _block_size_max, float _particle_spacing)
		{
			mParticleSpacing = _particle_spacing;
			mBlockSizeMin = _block_size_min;
			mBlockSizeMax = _block_size_max;
		}

		void BlockEmitter::reset(float4* positions, int count)
		{
			int num_particles_used = 0;
			float posY = mBlockSizeMin.y;
			while (num_particles_used < count)
			{
				for (float z = mBlockSizeMin.z; z < mBlockSizeMax.z; z += mParticleSpacing)
				{
					for (float x = mBlockSizeMin.x; x < mBlockSizeMax.x; x += mParticleSpacing)
					{
						// Verificar se atingiu o limite maximo de particulas
						if (num_particles_used < count)
						{
							positions[num_particles_used] = make_float4(x, posY, z, 1); // Criar part�culas com lifetime infinito
							num_particles_used++;
						}
					}
				}
				posY += mParticleSpacing;
			}

			// Atualizar particulas restantes como inutilizadas
			for (;num_particles_used < count; num_particles_used++)
			{
				positions[num_particles_used] = make_float4(0.0f, 0.0f, 0.0f, -1.0f);
			}

		}
	}
}
