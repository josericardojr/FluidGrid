#ifndef __K_FLUIDSIM_CU__

#define __K_FLUIDSIM_CU__

#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_inline.h>
#include <thrust/sort.h>
#include "K_SimulationLib.cuh"
#include "K_DeviceFunctions.cu"
#include <cutil_math.h>
#include "K_Boundaries_Walls.cu"
#include "../ApplicationData.cuh"

	
/*
 * Texture mapping
 */
#if USE_TEX
texture<float4, 1, cudaReadModeElementType> oldPosTex;
texture<float4, 1, cudaReadModeElementType> oldVelTex;
texture<float2, 1, cudaReadModeElementType> densityFieldTex;
texture<float4, 1, cudaReadModeElementType> forceTex;
texture<uint, 1, cudaReadModeElementType> gridParticleHashTex;
texture<uint, 1, cudaReadModeElementType> cellStartTex;
texture<uint, 1, cudaReadModeElementType> cellEndTex;
#endif

// Constantes
namespace SimulationLib
{

	namespace Fluid
	{

		/* Armazena kernels pre-calculados */
		__device__ __constant__ Simulation::Application::KernelsCache dKernelsCache;
		
		/* Armazena parametros comuns da simulacao de fluidos */
		__device__ __constant__ Simulation::Application::CommonFluidParams dCommonFluidParams;


		/* 
 		 * Processamento da densidade-pressao em uma celula
 		 */
		__device__ float ProcessMassDensityGridD( float4*     oldPos,
												unsigned int  gridHash,
												unsigned int* cellStart,
												unsigned int* cellEnd,
												float4        pPos,
												unsigned int  particleIdx,
												float         mass)
		{
			float mass_density = 0;
			uint startIndex = FETCH(cellStart, gridHash);
		
			if (startIndex != 0xffffffff)
			{
				float3 pos1_real = make_float3(pPos) * dCommonFluidParams.simulation_scale;

				uint endIndex = FETCH(cellEnd, gridHash);
		
				for (int i = startIndex; i < endIndex; i++)
				{
					float3 pos2_real = make_float3(FETCH(oldPos, i)) *
							dCommonFluidParams.simulation_scale;
		
					// Calcular diferenca
					float3 diff = pos1_real - pos2_real;
					float l = length(diff);
		
					if (l < dCommonFluidParams.kernel_radius)
						mass_density += mass * dKernelsCache.kPoly6 *
							powf( (dCommonFluidParams.kernel_radius_squared - (l * l)), 3);
				}	
			}

			return mass_density;
		}
		
		
		/*
		 * Metodo para o calculo da forca de viscosidade
		 */
		__device__ float3 SPHCalculateViscosityForceD(float3 pos1, float3 pos2, float4 v1, float4 v2,
													  float mass, float2 dp1, float2 dp2, float kViscosity, float _in_kernel_radius)
		{
			float3 diff = pos1 - pos2;

			float l = length(diff);

			if (l > _in_kernel_radius) return make_float3(0.0f);


			float lap = dKernelsCache.kViscosityLaplacian * (_in_kernel_radius - l);

			float3 f = make_float3( ((v2.x - v1.x) / (dp1.x * dp2.x)) * lap * kViscosity,
				    ((v2.y - v1.y) / (dp1.x * dp2.x)) * lap * kViscosity,
				    ((v2.z - v1.z) / (dp1.x * dp2.x)) * lap * kViscosity);

			return f;
		}

		/*
		 * Metodo para o calculo da pressao
		 */
		__device__ float3 SPHCalculatePressureForceD(float3 p1, float3 p2, float2 dp1, float2 dp2, float mass, float _in_kernel_radius)
		{
			float p =  ((dp1.y + dp2.y) / (dp1.x * dp2.x));

			float3 diff = p1 - p2;
			float l = length(diff);

			if (l > _in_kernel_radius) return make_float3(0.0f);

			float r = dKernelsCache.kSpikyGradient * ( (_in_kernel_radius - l) * (_in_kernel_radius - l));
			diff = (diff / l) * r * p * -1.0;

			return diff;
		}



		/* 
		 * Metodo para o processamento das forcas internas entre particulas da grid
		 */
		__device__ float3 SPHCalculateInternalForcesD( unsigned int  particleIdx,
													   unsigned int  gridHash,
													   float         scale,
													   float         _in_kernel_radius,
													   float         mass,
													   float         kViscosity,
													   float2*       densityField,
													   float4*       oldPos,
													   float4*       oldVel,
													   unsigned int* cellStart,
													   unsigned int* cellEnd,
													   float4        pos)
		{
			float3 total_forces = make_float3(0.0f);
			uint startIndex = FETCH(cellStart, gridHash);

			// Primeiramente deve ser verificado se a cï¿½lula possui partï¿½culas
			if (startIndex!= 0xffffffff)
			{
				float3 pos1_real = make_float3(pos) * scale; // posiï¿½ï¿½o real da partï¿½cula
				float2 dp1 = FETCH(densityField, particleIdx); // Densidade/pressï¿½o da partï¿½cula
				float4 v1 = FETCH(oldVel, particleIdx);          // Velocidade da partï¿½cula

				uint endIndex = FETCH(cellEnd, gridHash);

				// Percorrer todos as partï¿½culas nesta cï¿½lula
				for (int i = startIndex; i < endIndex; i++)
				{
					// Para o calculo das forcas de pressao e viscosidade, a particula em questao nao pode ser processada com ela mesma
					if ( i != particleIdx )
					{
						float3 pos2_real = make_float3(FETCH(oldPos, i)) * scale;
						float2 dp2 = FETCH(densityField, i);
						float4 v2 = FETCH(oldVel, i);

						// Processamento da forca de pressao
						total_forces += SPHCalculatePressureForceD(pos1_real, pos2_real, dp1, dp2, mass, _in_kernel_radius);

						// Processamento da forca de viscosidade
						total_forces += SPHCalculateViscosityForceD(pos1_real, pos2_real, v1, v2, mass, dp1, dp2, kViscosity, _in_kernel_radius);
					}
				}
			}

			return total_forces;
		}

		/*
		 * Metodo auxiliar responsavel pelo calculo das forcas internas em uma particula
		 * (Pressao e Viscosidade)
		 */
		__device__ float3 CalculateFluidInternalForcesD(float4& pos, int particleIndex, float mass,
				float viscosity, float2* densityField, float4* oldPos, float4* oldVel,
				unsigned int* cellStart, unsigned int* cellEnd)
		{
			float3 total_forces = make_float3(0.0);

			/* Calculo das forca internas (pressao e viscosidade) */
			for (int z = -1; z <= 1; z++)
			{
				for ( int y = -1; y <= 1; y++)
				{
					for (int x = -1; x <= 1; x++)
					{
						float3 offset = dCommonFluidParams.cell_size;
						offset.x *= x; offset.y *=y; offset.z *= z;

						unsigned int hashCell = SimulationLib::Common::HashKeyD(make_float3(pos) +
								offset, dCommonFluidParams.cell_size, dCommonFluidParams.getNumBuckets() );

						 total_forces += SPHCalculateInternalForcesD(particleIndex, hashCell,
								 dCommonFluidParams.simulation_scale, dCommonFluidParams.kernel_radius,
								 mass, viscosity, densityField, oldPos, oldVel, cellStart, cellEnd, pos);
					}
				}
			}

			return total_forces;
		}

		/*
		 * Calculo da densidade/pressao 
		 */
		__global__ void CalculateMassDensityK( float4*       oldPos,
											   float2*       mass_pressure,
											   unsigned int* hashTable,
											   unsigned int* gridStart,
											   unsigned int* gridEnd,
											   unsigned int  numParticles,
											   float         gasConstant,
											   float         restDensity,
											   float         mass)
		{
			unsigned int particleIndex = blockDim.x * blockIdx.x + threadIdx.x;	
		
			if (particleIndex < numParticles)
			{
				float4 pos = FETCH(oldPos, particleIndex);
				
				// Verificar se esta particula esta ativa
				if (pos.w > 0.0)
				{
					float mass_density = 0;
		
					for (int z = -1; z <= 1; z++)
					{
						for ( int y = -1; y <= 1; y++)
						{
							for (int x = -1; x <= 1; x++)
							{
								float3 offset = dCommonFluidParams.cell_size;
								offset.x *= x; offset.y *= y; offset.z *= z;
								unsigned int hashCell = SimulationLib::Common::HashKeyD(make_float3(pos) + offset,
										dCommonFluidParams.cell_size, dCommonFluidParams.getNumBuckets() );
													
								 mass_density += ProcessMassDensityGridD(oldPos,
										 hashCell, gridStart, gridEnd, pos, particleIndex, mass);
							}
						}
					}	
		
					float2 mp = make_float2(max(1.0, mass_density),
							gasConstant * (mass_density - restDensity));
		
					mass_pressure[particleIndex] = mp;
				}
			}
		}	


		/*
		 * Metodo para o calculo de forcas sem corpos rigidos
		 */
		__global__ void SPHCalculateForcesK(float4*       oldPos,
											float4*       oldVel,
											float4*       forces,
											float2*       densityField,
											unsigned int* cellStart,
											unsigned int* cellEnd,
											float4*       newVel,
											unsigned int* particleId,
											float         mass,
										    unsigned int  numParticles,
										    float         viscosity)
		{
			unsigned int particleIndex = blockDim.x * blockIdx.x + threadIdx.x;

			if (particleIndex < numParticles)
			{
				float4 pos = FETCH(oldPos, particleIndex);
				float4 vel = FETCH(oldVel, particleIndex);
				float2 dp = FETCH(densityField, particleIndex);
				float3 total_forces = make_float3(0.0);

				// Verificar se esta particula esta ativa
				if (pos.w > 0.0)
				{
					total_forces = CalculateFluidInternalForcesD(pos, particleIndex,
						mass, viscosity, densityField, oldPos, oldVel, cellStart, cellEnd);

					total_forces = total_forces * mass + dCommonFluidParams.global_external_force;

					float3 wallAvoid = calculateWallsNoPenetrationForce(make_float3(pos), make_float3(vel), dCommonFluidParams.volume_min,
						dCommonFluidParams.volume_max, 0.003043174, 20000, 256, dCommonFluidParams.simulation_scale);

					wallAvoid += calculateWallsNoSlipForce(make_float3(pos), make_float3(vel), total_forces + wallAvoid,
							make_float3(-11.0f, -40.0f, -11.0f), make_float3(11.0f, 200.0f, 11.0f), 0.003043174, 0, 0, dCommonFluidParams.simulation_scale);

					total_forces += wallAvoid;

					float speed = length(total_forces);

					if (speed > dCommonFluidParams.velocity_limit)
					{
						float t = dCommonFluidParams.velocity_limit / speed;
						total_forces *= t;
					}
				}

				unsigned int originalIndex = particleId[particleIndex];
				newVel[originalIndex] = vel + make_float4((total_forces  * dCommonFluidParams.deltaTimeIntegration));
			}
		}


		/*
		 * Metodo para o calculo das forcas com corpos rigidos
		 */
		__global__ void SPHCalculateForcesWithSRBK(
				/* Parametros do fluido */
				float4*       oldPos,
				float4*       oldVel,
				float4*       newVel,
				float4*       forces,
				float2*       density_pressure,
				unsigned int* particleId,
				unsigned int* cellStart,
				unsigned int* cellEnd,
				unsigned int  numParticles,
				float         mass,
				float         kViscosity,
				float         fluidParticleRadius,

				/* Parametros do corpo rigido estatico */
				float4*       _in_srb_positions,       // Posição das partículas dos corpos rígidos staticos (w = raio partícula)
				uint*         _in_srb_start_offset,     // Offset inicial de cada bucket
				uint*         _in_srb_end_offset,       // Offset final de cada bucket
				uint          _in_srb_num_buckets,      // Número de Buckets
				float         _in_srb_particle_radius,
				float         _in_srb_kspring,
				float         _in_srb_kdamping,
				float3          _in_srb_cell_size,
				int3          _in_srb_cell_size_relation,

				float4*       _in_rb_positions,
				float4*       _in_rb_velocities,
				uint*         _in_rb_start_offset,
				uint*         _in_rb_end_offset,
				uint          _in_rb_num_buckets,

				bool useStaticRB,
				bool useDynamicRB)
		{
			unsigned int particleIndex = blockDim.x * blockIdx.x + threadIdx.x;

			if (particleIndex < numParticles)
			{
				float4 pos = oldPos[particleIndex];
				float4 vel = oldVel[particleIndex];
				float2 dp = density_pressure[particleIndex];

				float3 total_forces = make_float3(0.0f);


				// Verificar se esta particula esta ativa
				if (pos.w > 0.0)
				{
					/* Calculo das forca internas (pressao e viscosidade) */
					float3 internal_forces = CalculateFluidInternalForcesD(pos, particleIndex,
							mass, kViscosity, density_pressure, oldPos, oldVel, cellStart, cellEnd);


					float3 external_force = make_float3(0.0);

					/* Calculo da forca gerada pela colisao de fluido com corpos rigidos estaticos */
					if (useStaticRB)
						external_force += SimulationLib::Common::ProcessCollisionForcesD(pos,
							vel, fluidParticleRadius, _in_srb_positions, NULL,
							_in_srb_particle_radius, _in_srb_kspring, _in_srb_kdamping,
							_in_srb_cell_size_relation, _in_srb_start_offset,
							_in_srb_end_offset, _in_srb_cell_size, _in_srb_num_buckets,
							false, false, 0, 0.003); //todo colocar epsilon como variavel global

					/* C‡lculo da fora gerada pela colis‹o de fluido com corpos r’gidos dinamicos */
					if (useDynamicRB)
						external_force += SimulationLib::Common::ProcessCollisionForcesD(pos,
							vel, fluidParticleRadius, _in_rb_positions, _in_rb_velocities,
							_in_srb_particle_radius, _in_srb_kspring, _in_srb_kdamping,
							_in_srb_cell_size_relation, _in_rb_start_offset, _in_rb_end_offset,
							_in_srb_cell_size, _in_rb_num_buckets, true, false, 0, 0.003);


					total_forces = internal_forces * mass +
						dCommonFluidParams.global_external_force + external_force;

					float3 wallAvoid = calculateWallsNoPenetrationForce(make_float3(pos), make_float3(vel), dCommonFluidParams.volume_min,
						dCommonFluidParams.volume_max, 0.003043174, 20000, 256, dCommonFluidParams.simulation_scale);

					wallAvoid += calculateWallsNoSlipForce(make_float3(pos), make_float3(vel), total_forces + wallAvoid,
						make_float3(-11.0f, -40.0f, -11.0f), make_float3(11.0f, 200.0f, 11.0f), 0.003043174, 0, 0, dCommonFluidParams.simulation_scale);

					total_forces += wallAvoid;

					float speed = length(total_forces);

					if (speed > dCommonFluidParams.velocity_limit)
					{
						float t = dCommonFluidParams.velocity_limit / speed;
						total_forces *= t;
					}
				}

				unsigned int originalIndex = particleId[particleIndex];
				newVel[originalIndex] = vel + make_float4((total_forces  *
					dCommonFluidParams.deltaTimeIntegration));
			}
		}

		///*
		// * Kernel para executar a integraï¿½ï¿½o temporal
		// */
		__global__ void SPHIntegrateK(float4* _in_sph_positions,   // SPH Posiï¿½ï¿½o partï¿½culas
									  float4* _in_sph_velocity,                   // SPH Velocidade partï¿½culas
									  float4* _in_sph_forces,             // SPH Forï¿½as internas
									  float2*  _in_density_press,
									  unsigned int* _in_grid_start,
									  unsigned int* _in_grid_end,
									  unsigned int  numParticles)
		{
			unsigned int particleIdx = blockDim.x * blockIdx.x + threadIdx.x;


			if (particleIdx < numParticles)
			{

				// Verificar se particula esta ativa
				if (_in_sph_positions[particleIdx].w > 0.0)
				{
//					float3 accel = make_float3( _in_sph_forces[particleIdx].x, _in_sph_forces[particleIdx].y, _in_sph_forces[particleIdx].z);

	//				_in_sph_velocity[particleIdx] += make_float4( accel.x * _in_elapsed_time, accel.y * _in_elapsed_time, accel.z * _in_elapsed_time, 0.0);
					_in_sph_positions[particleIdx] +=  _in_sph_velocity[particleIdx] * ( dCommonFluidParams.deltaTimeIntegration / dCommonFluidParams.simulation_scale);

					//XSPHD(_in_sph_positions, _in_sph_velocity, _in_density_press, _in_grid_start, _in_grid_end,
						//_in_hash_table_size, _in_sph_num_particles, _in_cell_size, _in_mass, _in_scale, _in_kernel_radius);

					/* Test check */
					/*if (_in_sph_positions[particleIdx].y < _in_min_volume.y)
					{
						_in_sph_positions[particleIdx].y = _in_min_volume.y;
						_in_sph_velocity[particleIdx] .y *= -1;
					} else if (s_positions[threadIdx.x].y >= _in_max_volume.y)
					{
						s_positions[threadIdx.x].y = _in_max_volume.y;
						s_velocities[threadIdx.x].y *= -1;
					}

					if (_in_sph_positions[particleIdx].x < _in_min_volume.x)
					{
						_in_sph_positions[particleIdx].x = _in_min_volume.x + 0.05;
						_in_sph_velocity[particleIdx] .x *= -1;
					}
					else if (_in_sph_positions[particleIdx].x > _in_max_volume.x)
					{
						_in_sph_positions[particleIdx].x = _in_max_volume.x - 0.05;
						_in_sph_velocity[particleIdx].x *= -1;
					}

					if (_in_sph_positions[particleIdx].z < _in_min_volume.z)
					{
						_in_sph_positions[particleIdx].z = _in_min_volume.z + 0.05;
						_in_sph_velocity[particleIdx] .z *= -1;
					}
					else if (_in_sph_positions[particleIdx].z > _in_max_volume.z)
					{
						_in_sph_positions[particleIdx].z = _in_max_volume.z - 0.05;
						_in_sph_velocity[particleIdx] .z *= -1;
					}*/
				}
			}
		}

	}
}




/*
 * Funcao externa para o calculo da densidade
 */		
extern "C"		
	void gCalcMassDensity(float4*       orderedPositions,
						  float2*       mass_pressure,
						  unsigned int* hashTable,
						  unsigned int* gridStart,
						  unsigned int* gridEnd,
						  int           numParticles,
						  float         gasConstant,
						  float         restDensity,
						  float         mass)
	{
		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)numParticles / blockSize.x)), 1, 1);

	 	Simulation::Application::CommonFluidParams& commonFluidParams =
	 			Simulation::Application::ApplicationData::instance()->getCommonFluidParams();
#if USE_TEX
	    cutilSafeCall(cudaBindTexture(0, oldPosTex, orderedPositions, numParticles*sizeof(float4)));
	    cutilSafeCall(cudaBindTexture(0, cellStartTex, gridStart, commonFluidParams.getNumBuckets()*sizeof(uint)));
	    cutilSafeCall(cudaBindTexture(0, cellEndTex, gridEnd, commonFluidParams.getNumBuckets()*sizeof(uint)));
#endif

		SimulationLib::Fluid::CalculateMassDensityK<<<gridSize, blockSize>>>(orderedPositions, mass_pressure,
			hashTable, gridStart, gridEnd, numParticles, gasConstant, restDensity, mass);

#if USE_TEX
	    cutilSafeCall(cudaUnbindTexture(oldPosTex));
	    cutilSafeCall(cudaUnbindTexture(cellStartTex));
	    cutilSafeCall(cudaUnbindTexture(cellEndTex));
#endif
	}	
	
	



/*
 * Calcular forcas sem considerar corpos rigidos
 */
extern "C"
	void gSPHCalculateForces(float4*       orderedPositions,
				  			 float4*       orderedVelocity,
							 float4*       forces,
							 float2*       density_pressure,
							 unsigned int* gridStart,
							 unsigned int* gridEnd,
							 float4*       vel,
							 unsigned int* particleId,
							 unsigned int  numParticles,
							 float         mass,
							 float         viscosity)

	{
		dim3 blockSize(128, 1, 1);
	 	dim3 gridSize((int)(ceil((float)numParticles / blockSize.x)), 1, 1);

	 	Simulation::Application::CommonFluidParams& commonParams =
	 			Simulation::Application::ApplicationData::instance()->getCommonFluidParams();

#if USE_TEX
    	cutilSafeCall(cudaBindTexture(0, oldPosTex, orderedPositions, numParticles*sizeof(float4)));
    	cutilSafeCall(cudaBindTexture(0, oldVelTex, orderedVelocity, numParticles*sizeof(float4)));
    	cutilSafeCall(cudaBindTexture(0, densityFieldTex, density_pressure, numParticles*sizeof(float2)));
    	cutilSafeCall(cudaBindTexture(0, cellStartTex, gridStart, commonParams.getNumBuckets()*sizeof(uint)));
    	cutilSafeCall(cudaBindTexture(0, cellEndTex, gridEnd, commonParams.getNumBuckets()*sizeof(uint)));
#endif

		SimulationLib::Fluid::SPHCalculateForcesK<<<gridSize, blockSize>>>(
				orderedPositions, orderedVelocity, forces, density_pressure, gridStart,
				gridEnd, vel, particleId, mass, numParticles, viscosity);

#if USE_TEX
    	cutilSafeCall(cudaUnbindTexture(oldPosTex));
    	cutilSafeCall(cudaUnbindTexture(oldVelTex));
    	cutilSafeCall(cudaUnbindTexture(densityFieldTex));
    	cutilSafeCall(cudaUnbindTexture(cellStartTex));
    	cutilSafeCall(cudaUnbindTexture(cellEndTex));
#endif
	}

/*
 * Integracao temporal das particulas do SPH
 */
extern "C"
	void gSPHIntegrate(float4* _in_sph_positions,   // SPH Posiï¿½ï¿½o partï¿½culas
					   float4* _in_sph_velocity,                   // SPH Velocidade partï¿½culas
					   float4* _in_sph_forces,             // SPH Forï¿½as internas
					   uint    _in_num_particles,
					   float2*  _in_density_press,
					   unsigned int* _in_grid_start,
					   unsigned int* _in_grid_end)
	{
		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)_in_num_particles / blockSize.x)), 1, 1);
		size_t s_mem = sizeof(float4) * 3 * blockSize.x;

		SimulationLib::Fluid::SPHIntegrateK<<<gridSize, blockSize>>>(_in_sph_positions,
			_in_sph_velocity, _in_sph_forces, _in_density_press, _in_grid_start, _in_grid_end, _in_num_particles);
	}

/*
 * Funcao externa para inicializacao dos parametros comuns do fluido
 */
extern "C"
	void gInitializeCommonFluidParams(Simulation::Application::CommonFluidParams hFluidParams)
	{
		cutilSafeCall(cudaMemcpyToSymbol(SimulationLib::Fluid::dCommonFluidParams,
				&hFluidParams, sizeof(Simulation::Application::CommonFluidParams)));

	}

/*
 * Funcao externa para inicializar kernels
 */
extern "C"
	void gInitializeKernelParams(Simulation::Application::KernelsCache _hKernelsCache)
	{
		//copiar valores do kernel pre calculado
		cutilSafeCall(
			cudaMemcpyToSymbol(SimulationLib::Fluid::dKernelsCache, &_hKernelsCache,
					sizeof(Simulation::Application::KernelsCache)));
	}

/*
 * Funcao para o calculo das forcas internas do fluido e tambem considerando os corpos estaticos
 * na cena
 */
extern "C"
	void gSPHCalculateForcesWithSRB(float4*       oldPositions,
				  			        float4*       oldVelocity,
							 		float4*       forces,
							 		float4*       newVel,
									float2*       density_pressure,
									unsigned int* particleId,
									unsigned int* gridStart,
									unsigned int* gridEnd,
									unsigned int  numParticles,
									float         mass,
							        float         kViscosity,
									float         _in_sph_particle_radius,
									int3          _in_cell_size_relation,

									float         _in_srb_particle_radius,
									float3        _in_srb_cell_size,
									float         _in_kspring,
									float         _in_kdamping,

									float4*       _in_srb_positions,       // Posição das partículas dos corpos rígidos staticos (w = raio partícula)
									uint*         _in_srb_start_offset,     // Offset inicial de cada bucket
									uint*         _in_srb_end_offset,       // Offset final de cada bucket
									uint          _in_srb_num_buckets,

									float4*       _in_rb_positions,
									float4*       _in_rb_velocities,
									uint*         _in_rb_start_offset,
									uint*         _in_rb_end_offset,
									uint          _in_rb_num_buckets,

									bool useStaticRB,
									bool useDynamicRB)

	{

		dim3 blockSize(128, 1, 1);
	 	dim3 gridSize((int)(ceil((float)numParticles / blockSize.x)), 1, 1);

	 	Simulation::Application::CommonFluidParams& commonParams =
	 			Simulation::Application::ApplicationData::instance()->getCommonFluidParams();

#if USE_TEX
    	cutilSafeCall(cudaBindTexture(0, oldPosTex, oldPositions, numParticles*sizeof(float4)));
    	cutilSafeCall(cudaBindTexture(0, oldVelTex, oldVelocity, numParticles*sizeof(float4)));
    	cutilSafeCall(cudaBindTexture(0, densityFieldTex, density_pressure, numParticles*sizeof(float2)));
    	cutilSafeCall(cudaBindTexture(0, cellStartTex, gridStart, commonParams.getNumBuckets()*sizeof(uint)));
    	cutilSafeCall(cudaBindTexture(0, cellEndTex, gridEnd, commonParams.getNumBuckets()*sizeof(uint)));
#endif

		SimulationLib::Fluid::SPHCalculateForcesWithSRBK<<<gridSize, blockSize>>>(oldPositions, oldVelocity, newVel,
			forces, density_pressure, particleId, gridStart, gridEnd, numParticles, mass,
			kViscosity, _in_sph_particle_radius, _in_srb_positions, _in_srb_start_offset,
			_in_srb_end_offset, _in_srb_num_buckets, _in_srb_particle_radius, _in_kspring,
			_in_kdamping, _in_srb_cell_size, _in_cell_size_relation, _in_rb_positions,
			_in_rb_velocities, _in_rb_start_offset, _in_rb_end_offset, _in_rb_num_buckets,
			useStaticRB, useDynamicRB);

#if USE_TEX
    	cutilSafeCall(cudaUnbindTexture(oldPosTex));
    	cutilSafeCall(cudaUnbindTexture(oldVelTex));
    	cutilSafeCall(cudaUnbindTexture(densityFieldTex));
    	cutilSafeCall(cudaUnbindTexture(cellStartTex));
    	cutilSafeCall(cudaUnbindTexture(cellEndTex));
#endif
	}
	
/*
 * Funcao para realizacao de sort
 */
extern "C" 
	void gSort(unsigned int* _key, unsigned int* _values, int _numParticles)
	{
		thrust::sort_by_key(thrust::device_ptr<unsigned int>(_key),
				thrust::device_ptr<unsigned int>(_key + _numParticles),
				thrust::device_ptr<unsigned int>(_values));
	}



#endif
