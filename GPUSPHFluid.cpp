#include "GPUSPHFluid.h"
#include <cutil_inline.h>


namespace Simulation
{
	namespace Fluid
	{
		GPUSPHFluid::GPUSPHFluid(Simulation::Application::FluidParams _fluidParams, Emitter::IEmitter* _emitter)
			: mEmitter(_emitter),
			  mFluidParams(_fluidParams)
		{
			Application::CommonFluidParams& commonFluidParams =
					Application::ApplicationData::instance()->getCommonFluidParams();

			// Inicializar estrutura de dados
			cutilSafeCall(cudaMalloc((void**)&d_particleId, sizeof(unsigned int) * mFluidParams.max_particles));
			cutilSafeCall(cudaMalloc((void**)&d_hashKey, sizeof(unsigned int) * mFluidParams.max_particles));
			cutilSafeCall(cudaMalloc((void**)&d_gridPosStart, sizeof(unsigned int) * commonFluidParams.getNumBuckets()));
			cutilSafeCall(cudaMalloc((void**)&d_gridPosEnd, sizeof(unsigned int) * commonFluidParams.getNumBuckets()));
			cutilSafeCall(cudaMalloc((void**)&d_ordered_position, sizeof(float4) * mFluidParams.max_particles));
			cutilSafeCall(cudaMalloc((void**)&d_forces, sizeof(float4) * mFluidParams.max_particles));
			cutilSafeCall(cudaMalloc((void**)&d_mass_pressure, sizeof(float2) * mFluidParams.max_particles));
			cutilSafeCall(cudaMalloc((void**)&d_velocity, sizeof(float4) * mFluidParams.max_particles));
			cutilSafeCall(cudaMalloc((void**)&d_ordered_vel, sizeof(float4) * mFluidParams.max_particles));

			// Resetar valores
			cutilSafeCall(cudaMemset(d_velocity, 0, sizeof(float4) * mFluidParams.max_particles));
			cutilSafeCall(cudaMemset(d_ordered_vel, 0, sizeof(float4) * mFluidParams.max_particles));

			JRFXGL::Util::LogManager::instance()->log("Dados criados");

			// Inicializar emissor
			float4* positions = new float4[mFluidParams.max_particles];
			mEmitter -> reset(positions, mFluidParams.max_particles);

			initializeVBO(positions, mFluidParams.max_particles);
			delete positions;

			// Inicializar interop
			cudaBuffer.setBuffer(&getVBOPosition());
		}

		GPUSPHFluid::~GPUSPHFluid()
		{
			cudaFree(d_particleId);
			cudaFree(d_hashKey);
			cudaFree(d_gridPosEnd);
			cudaFree(d_gridPosStart);
			cudaFree(d_ordered_position);
			cudaFree(d_forces);
			cudaFree(d_mass_pressure);
			cudaFree(d_velocity);
			cudaFree(d_ordered_vel);
		}


		void GPUSPHFluid::beginUpdate()
		{
			cudaBuffer.mapBuffer();

		}

		void GPUSPHFluid::endUpdate()
		{
			cudaBuffer.unmapBuffer();
		}

		float4* GPUSPHFluid::getPosition()
		{
			return (float4*)(cudaBuffer.getData());
		}

		void GPUSPHFluid::render()
		{
			glEnableClientState(GL_VERTEX_ARRAY);
			glEnable(GL_POINT_SPRITE);
			glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
			glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
			getVBOPosition().bindBuffer();
			glVertexPointer(4, GL_FLOAT, 0, 0);
			mMaterial->apply();
			glColor3f(1, 0, 0);
			glDrawArrays(GL_POINTS, 0, mFluidParams.max_particles);
			getVBOPosition().unbindBuffer();
			mMaterial->disable();
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisable(GL_POINT_SPRITE);
		}
	}
}
