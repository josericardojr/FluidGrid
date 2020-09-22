/*
 * ApplicationData.h
 *
 *  Created on: 31/01/2011
 *      Author: josericardodasilvajunior
 */

#ifndef __ApplicationData_H__
#define __ApplicationData_H__

#include "vector_types.h"
#include <vector_functions.h>

namespace Simulation
{
	namespace Application
	{
		static const float PI = 3.14159265f;
		
		/*
		 * Parametros utilizados para a simulacao
		 */
		struct SimulationParams
		{
			float              dampingForce;
			float              springForce;
		};

		/*
		 * Parametros da simulacao
		 */
		struct FluidParams
		{
			float               mass;          // Massa de cada particula
			float               restDensity;
			float               viscosity;
			float               stiffness;
			float               particle_radius; // Utilizado para colisao com corpos rigidos
			float               gasConstant;
			int                 max_particles;
		};

		/*
		 * Armazena os parametros cumuns a simulacao de fluidos
		 */
		struct CommonFluidParams
		{
			float3              cell_size;
			int3                bucket_size;
			unsigned int        keyBits;
			float               kernel_radius;
			float               kernel_radius_squared;
			float               pointRadius;
			float               pointScale;
			float               velocity_limit;
			float3              volume_min;
			float3              volume_max;
			float3              global_external_force;
			float               simulation_scale;
			float               deltaTimeIntegration;
			int3                fluid_to_srb_relation;

			__device__ __host__ int getNumBuckets()
			{
				return bucket_size.x * bucket_size.y * bucket_size.z;
			}
		};


		/*
		 * Armazena os parametros cumuns a simulacao de corpos rigidos
		 */
		struct CommonRigidBodyParams
		{
			int3                rb_to_fluid_relation;
			int3                st_bucket_size;
			int3                dn_bucket_size;
			unsigned int        st_keyBits;
			float               velocity_limit;
			float3              volume_min;
			float3              volume_max;
			float4              global_external_force;
			float               simulation_scale;
			float               deltaTimeIntegration;
			float               particleRadius;


			__device__ __host__ int getStNumBuckets()
			{
				return st_bucket_size.x * st_bucket_size.y * st_bucket_size.z;
			}

			__device__ __host__ int getDnNumBuckets()
			{
				return dn_bucket_size.x * dn_bucket_size.y * dn_bucket_size.z;
			}
			
			__device__ __host__ float3 getCellSize()
			{
				return make_float3(particleRadius, particleRadius, particleRadius);
			}
		};

		/*
		 * Armazena informacoes pre-calculadas dos kernels utilizados
		 */
		struct KernelsCache
		{
			float kPoly6;
			float kSpikyGradient;
			float kViscosityLaplacian;
		};


		/*
		 * Classe singleton utilizada para armazenamento dos parametros da simulacao
		 */
		class ApplicationData
		{
		public:
			static ApplicationData* instance();

		public:
			CommonFluidParams& getCommonFluidParams(){ return mCommonFluidParams; }
			KernelsCache&      getKernelsCache(){ return mKernelsCache; }
			CommonRigidBodyParams &getCommonRigidBodyParams(){ return mCommonRigidBodyParams; }
			SimulationParams&     getSimulationParams(){ return mSimulationParams; }

			void setSimulationParams(SimulationParams _simulationParams)
			{
				mSimulationParams = _simulationParams;
			}
			
			void setCommonFluidParams(CommonFluidParams _commonFluidParams)
			{
				mCommonFluidParams = _commonFluidParams;
			}

			void setCommonRigidBodyParams(CommonRigidBodyParams _commonRigidBodyParams)
			{
				mCommonRigidBodyParams = _commonRigidBodyParams;
			}

			void setKernelsCache(KernelsCache _kernelsCache)
			{
				mKernelsCache = _kernelsCache;
			}


		private:
			ApplicationData(){}

		private:
			static ApplicationData*  mInstance;
			CommonFluidParams        mCommonFluidParams;
			CommonRigidBodyParams    mCommonRigidBodyParams;
			SimulationParams         mSimulationParams;
			KernelsCache             mKernelsCache;
		};
	}
}


#endif /* APPLICATIONDATA_H_ */
