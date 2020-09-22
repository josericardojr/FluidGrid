#ifndef __GPUSPHFluidSim_H__

#define __GPUSPHFluidSim_H__

/*
 * Classe que representa a simulaÁ„o de um tipo de fluido 
 */
//#include "cudpp.h"
#include <thrust/device_vector.h>
#include "ApplicationData.cuh"
#include "Emitter.h"
#include "FluidRenderer.h"
#include "GLCudaBuffer.h"



namespace Simulation
{

	namespace Emitter{ class IEmitter;}

	namespace Fluid
	{
		class GPUSPH_Simulator;

		class GPUSPHFluid : public FluidRenderer
		{
		public:
			// Construtor
			GPUSPHFluid(Application::FluidParams _fluidParams, Emitter::IEmitter* _emitter);

			// Destrutor
			~GPUSPHFluid();

			/*
			 * Método que deve ser chamado antes dos updates
			 */
			void beginUpdate();

			/*
			 * Metodo que dever ser chamado apos updates
			 */
			void endUpdate();

			/*
			 * MÈtodo para calcular a localizaÁ„o de vizinhanÁa
			 */
//			void Find_Neighbourhood();

			/*
			 * MÈtodo para calcular a densidade/pressao
			 */
//			void Calculate_Density();


			/*
			 * MÈtodo para calcular forÁas internas e externas em CPU
			 */
		//	void Calculate_Forces(uint _time_stamp, float _phys_particle_rad, CollisionStructures::HashTable* _srb_hash_table,
		//		float3 _in_sphere_collection_cell_size, float4* _in_srb_pos, float _in_phys_kSpring, float _in_phys_kDamping, int3 _in_cell_size_relation,
		//		CollisionStructures::HashTable* _rb_hash_table, float4* _in_rb_pos, float4* _in_rb_vel, float3 _global_external_force)
		//	{
		//		mSPHSimulator -> Calculate_Forces(mFluidData, mFluidParams.viscosity, _phys_particle_rad, _in_sphere_collection_cell_size,
		//			_srb_hash_table, _in_srb_pos, _in_phys_kSpring, _in_phys_kDamping, _in_cell_size_relation, _rb_hash_table,
		//			_in_rb_pos, _in_rb_vel, _global_external_force, _time_stamp);
		//	}

			/*
			 * MÈtodo para calcular forÁas internas e externas em GPU
			 */
//			void Calculate_Forces();
	//

			/*
			 * MÈtodo para calcular forÁas internas e externas em GPU
			 */
	//		void Calculate_Forces(float _phys_particle_rad, float _in_phys_kSpring, float _in_phys_kDamping,
	//			int3 _in_cell_size_relation, float3 _in_phys_cell_size, float4* _in_srb_pos, uint* _in_srb_start_offset,
	//			uint* _in_srb_end_offset, uint _in_srb_num_buckets, float3 _global_external_force, float4* _in_rb_pos,
	//			float4* _in_rb_vel, uint* _in_rb_start_offset, uint* _in_rb_end_offset, uint _in_rb_num_buckets)
	//		{
	//			mSPHSimulator -> Calculate_Forces(mFluidData, mFluidParams.viscosity, _phys_particle_rad, _in_phys_kSpring, _in_phys_kDamping,
	//				_in_phys_cell_size, _in_srb_pos, _in_cell_size_relation, _in_srb_start_offset, _in_srb_end_offset,
	//				_in_srb_num_buckets, _global_external_force, _in_rb_pos, _in_rb_vel, _in_rb_start_offset,
	//				_in_rb_end_offset, _in_rb_num_buckets);
	//		}

			// Calcular forÁas somente considerando corpos rÌgidos estaticos
	//		void Calculate_Forces(float _phys_particle_rad, float _in_phys_kSpring, float _in_phys_kDamping,
	///			int3 _in_cell_size_relation, float3 _in_phys_cell_size, float4* _in_srb_pos, uint* _in_srb_start_offset,
		//		uint* _in_srb_end_offset, uint _in_srb_num_buckets, float3 _global_external_force)
	//		{
	///			mSPHSimulator -> Calculate_Forces(mFluidData, mFluidParams.viscosity, _phys_particle_rad, _in_phys_kSpring, _in_phys_kDamping,
	//				_in_phys_cell_size, _in_srb_pos, _in_cell_size_relation, _in_srb_start_offset, _in_srb_end_offset,
	//				_in_srb_num_buckets, _global_external_force);
	//		}


			/*
			 * MÈtodo para calcular forÁas internas e externas em GPU
			 */
	//		void Calculate_Forces(float3 _global_external_force)
	//		{
	//			mSPHSimulator -> Calculate_Forces(mFluidData, mFluidParams.viscosity, _global_external_force);
	//		}

			/*
			 * MÈtodo para IntegraÁ„o temporal
			 */
//			void Integrate();


			//Engine::BoundingBox& GetBoundingBox(){ return mBoundingBox; }


			/* Métodos auxiliares */
			unsigned int* getHashKey(){ return d_hashKey; }
			unsigned int* getParticleID(){ return d_particleId; }
			unsigned int* getHashStart(){ return d_gridPosStart; }
			unsigned int* getHashEnd(){ return d_gridPosEnd; }
			float4*       getPosition();
			float4*       getOrderedPosition(){ return d_ordered_position; }
			float4*       getOrderedVelocities(){ return d_ordered_vel; }
			float4*       getVelocity(){ return d_velocity; }
			float4*       getForces(){ return d_forces; }
			float2*       getDensityField(){ return d_mass_pressure; }
			int           getNumParticles(){ return mFluidParams.max_particles; }
			float         getGasConstant(){ return mFluidParams.gasConstant; }
			float         getRestDensity(){ return mFluidParams.restDensity; }
			float         getMass(){ return mFluidParams.mass; }
			float         getViscosity(){ return mFluidParams.viscosity; }
			//vl::AABB      getAABB(){ return vl::AABB(vl::fvec3(-1000, -1000, -1000), vl::fvec3(1000, 1000, 1000)); }

			//CUDPPHandle   getPlanHandle(){ return planHandle; }

			virtual void render();

		public:
			int                               mId;
			Application::FluidParams          mFluidParams;
			Emitter::IEmitter*                mEmitter;

			// Informações do Fluido
			unsigned int*    d_particleId;
			unsigned int*    d_hashKey;
			unsigned int*    d_gridPosStart;
			unsigned int*    d_gridPosEnd;
			float4*          d_ordered_position;
			float4*          d_forces;
			float2*          d_mass_pressure;
			float4*          d_velocity;
			float4*          d_ordered_vel;

			//CUDPPHandle      planHandle;
			Simulation::GLCudaBufferVec4 cudaBuffer;
		};
	}
}

#endif
