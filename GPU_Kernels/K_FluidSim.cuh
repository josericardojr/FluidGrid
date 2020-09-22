#ifndef __K_FLUIDSIM_CUH__

#define __K_FLUIDSIM_CUH__

#include "../ApplicationData.cuh"
	
extern "C"		
	void gCalcMassDensity(float4*       orderedPositions,
						  float2*       mass_pressure,
						  unsigned int* hashTable,
						  unsigned int* gridStart,
						  unsigned int* gridEnd,
						  int           numParticles,
						  float         gasConstant,
						  float         restDensity,
						  float         mass);

extern "C"
	void gInitializeKernelParams(Simulation::Application::KernelsCache _hKernelsCache);

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
							 float         viscosity);

extern "C"
	void gSPHIntegrate(float4* _in_sph_positions,   // SPH Posição partículas
					   float4* _in_sph_velocity,                   // SPH Velocidade partículas
					   float4* _in_sph_forces,             // SPH Forças internas
					   unsigned int    _in_num_particles,
					   float2*  _in_density_press,
					   unsigned int* _in_grid_start,
					   unsigned int* _in_grid_end);

extern "C"
	void gInitializeCommonFluidParams(Simulation::Application::CommonFluidParams hFluidParams);

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
								bool useDynamicRB);
								
								
extern "C"
	void gSort(unsigned int* _key, unsigned int* _values, int _numParticles);

#endif
