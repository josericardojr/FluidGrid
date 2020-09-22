#ifndef __K_RIGIDBODYSIM_H__

#define __K_RIGIDBODYSIM_H__

extern "C"
	void gCalculateRigidBodyForcesFromRB( float4*       _in_RBPositions,    // Posi��es de cada RB part�cula
									float4*       _in_RBVelocites,
									float4*       _in_RBCMPosition,
									float4*       _out_torque,
									float4*       _out_Forces,        // Total das for�as em cada part�cula
									unsigned int* _in_RBGridStart,    // Offset inicial de cada c�lula (corpos dinamicos)
									unsigned int* _in_RBGridEnd,      // Offset final de cada c�lula (corpos dinamicos)
									unsigned int  _in_NumParticles,   // N�mero total de part�culas
									float3        _in_CellSize,       // Tamanho de cada c�lula
									float         _in_Radius,         // Raio de cada part�cula
									int           _in_RBHashTableSize,// Total de Buckets para RB
									int           _in_SRBHashTableSize, // Total de Buckets para SRB
									float4*       _in_SRBPositions,   // Posi��o de cada SRB part�cula
									unsigned int* _in_SRBGridStart,   // Offset inicial de cada c�lula (corpos est�ticos)
									unsigned int* _in_SRBGridEnd,     // Offset final de cada c�lula (corpos est�ticos)
									float         _in_RBkSpring,      // Constante de Spring
									float         _in_RBkDamping,     // Constante de damping
									int           _in_numSRB
									);

// M�todo para o calculo de for�a em cada part�cula de um corpo r�gido com os fluidos
extern "C"
	void gCalculateRigidBodyForcesFromFluids(float4*       _in_RBPositions,    // Posi��es de cada RB part�cula
										float4*       _in_RBVelocites,
										float4*       _in_RBCMPos,
										float4*       _out_Forces,        // Total das for�as em cada part�cula
										float4*       _out_torque,
										unsigned int* _in_RBGridStart,    // Offset inicial de cada c�lula (corpos dinamicos)
										unsigned int* _in_RBGridEnd,      // Offset final de cada c�lula (corpos dinamicos)
										unsigned int  _in_NumParticles,   // N�mero total de part�culas
										float3        _in_CellSize,       // Tamanho de cada c�lula
										float         _in_Radius,         // Raio de cada part�cula
										int           _in_RBHashTableSize,// Total de Buckets para RB

										float4*      _in_sph_particles,
										float4*      _in_sph_velocities,
										unsigned int*        _in_sph_grid_start,
										unsigned int*        _in_sph_grid_end,
										float        _in_sph_particle_radius,
										float3       _in_sph_cell_size,
										unsigned int         _in_sph_num_buckets,
										int3         _in_cell_size_relation,
					    				float        _in_sph_kSpring,
					    				float        _in_sph_kDamping
								);

extern "C"
	void gCopyRBDataToParticle(unsigned int* _out_force_rb,
			unsigned int* _out_force_index, float4* _in_rb_particle_forces,
			size_t num_particles);

extern "C"
	void gReorderForceTorque(float4* _out_force, float4* _out_torque, float4* _in_force, float4* _in_torque,
			unsigned int *_in_ordered_index, unsigned int _num_values);

extern "C"
	void gAddForceTorque(float4* _in_sum_force, float4* _in_sum_torque, float4* _in_ordened_force,
			float4* _out_to_force, float4* _out_to_torque, int* _in_sum_offset, size_t num_rigid_body);
extern "C"
	void gIntegrateRB(float4* _in_RBPosition, float* _in_mass, float4* _in_RBVelocities, float4* _in_rbforces, size_t num_rb,
			float _elapsed_time_in_seconds, float4 _in_external_force, float3 min_volume, float3 max_volume, float4* _angular_momentum,
			float4* _torque, float4* _inertia_tensor, float4* _orientation, float4* _inv_inertia_tensor, float velocity_limit);


#endif
