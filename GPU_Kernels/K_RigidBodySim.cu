#ifndef __K_RIGIDBODYSIM_CU__

#define __K_RIGIDBODYSIM_CU__

#include <cutil_math.h>
#include "K_DeviceFunctions.cu"
#include "K_MathUtils.cu"

// Constantes
namespace SimulationLib
{

	namespace RigidBody
	{
		/*
		 * Metodo responsavel pelo calculo de forcas entre particulas de corpos rigidos
		 */
		__global__ void CalculateRigidBodyForcesFromRBK(float4* _in_RBPositions, float4* _in_RBVelocites,
			float4* _in_RBCMPosition, float4* _out_torque, float4* _out_Forces,
			unsigned int* _in_RBGridStart,  unsigned int* _in_RBGridEnd,
			unsigned int  _in_NumParticles, float3 _in_CellSize, float  _in_Radius,
			int _in_RBHashTableSize, int _in_SRBHashTableSize, float4* _in_SRBPositions,
			unsigned int* _in_SRBGridStart, unsigned int* _in_SRBGridEnd,
			float _in_RBkSpring, float _in_RBkDamping, int _in_numSRB)

		{
			unsigned int particleIndex = blockDim.x * blockIdx.x + threadIdx.x;

			if (particleIndex < _in_NumParticles)
			{
				float4 pos = _in_RBPositions[particleIndex];
				float3 total_forces = make_float3(0.0f);
				float4 vel = _in_RBVelocites[(int)pos.w];
				float4 rb_cm_pos = _in_RBCMPosition[(int)pos.w];

				extern __shared__ unsigned int s_flags[];

				// Forca com corpos dinamicos
				total_forces += SimulationLib::Common::ProcessCollisionForcesD(pos, vel,
					_in_Radius, _in_RBPositions, _in_RBVelocites, _in_Radius, _in_RBkSpring,
					_in_RBkDamping, make_int3(1, 1, 1), _in_RBGridStart, _in_RBGridEnd,
					_in_CellSize, _in_RBHashTableSize, true, true, (int)pos.w, 0.001f);


				// Forca com corpos estaticos
				if (_in_numSRB > 0)
					total_forces += SimulationLib::Common::ProcessCollisionForcesD(pos, vel,
						_in_Radius, _in_SRBPositions, NULL, _in_Radius, _in_RBkSpring,
						_in_RBkDamping, make_int3(1, 1, 1), _in_SRBGridStart, _in_SRBGridEnd,
						_in_CellSize, _in_SRBHashTableSize, false, false, 0, 0.001f);


				_out_Forces[particleIndex] = make_float4(total_forces.x, total_forces.y, total_forces.z, 0);
				_out_Forces[particleIndex].w = (int)pos.w;

				// Torque
				float4 dir = pos - rb_cm_pos;
				_out_torque[particleIndex] = MathUtils::Vec3CrossD(dir, make_float4(total_forces.x, total_forces.y, total_forces.z, 0));
				_out_torque[particleIndex].w = (int)pos.w;
			}
		}

		/*
		 * Calcular forca do fluido nos corpos rigidos dinamicos
		 */
		__global__ void CalculateRBForcesFromFluidsK(float4* _in_RBPositions, float4* _in_RBVelocites,
			float4* _in_RBCMPos, float4* _out_Forces, float4* _out_torque,  unsigned int* _in_RBGridStart,
			unsigned int* _in_RBGridEnd, unsigned int  _in_NumParticles, float3 _in_CellSize, float _in_Radius,
			int _in_RBHashTableSize, float4* _in_sph_particles, float4* _in_sph_velocites, uint* _in_sph_grid_start,
			uint* _in_sph_grid_end, float _in_sph_particle_radius, float3 _in_sph_cell_size, uint _in_sph_num_buckets,
			int3 _in_cell_size_relation, float _in_sph_kSpring, float _in_sph_kDamping)

		{
			unsigned int particleIndex = blockDim.x * blockIdx.x + threadIdx.x;

			if (particleIndex < _in_NumParticles)
			{
				float4 pos = _in_RBPositions[particleIndex];
				float3 total_forces = make_float3(0.0f);
				float4 vel = _in_RBVelocites[(int)pos.w];
				float4 rb_cm_pos = _in_RBCMPos[(int)pos.w];

				// Força da ação dos fluidos
				total_forces += SimulationLib::Common::ProcessCollisionForcesD(pos, vel,
					_in_Radius, _in_sph_particles, _in_sph_velocites, _in_sph_particle_radius,
					_in_sph_kSpring, _in_sph_kDamping, _in_cell_size_relation, _in_sph_grid_start,
				    _in_sph_grid_end, _in_sph_cell_size, _in_sph_num_buckets, false, false, 0, 0.0001);


				_out_Forces[particleIndex] = make_float4(total_forces.x, total_forces.y, total_forces.z, 0.0f);
				_out_Forces[particleIndex].w = (int)pos.w;

				// Torque
				float4 dir = pos - rb_cm_pos;
				_out_torque[particleIndex] = SimulationLib::MathUtils::Vec3CrossD(
					dir, make_float4(total_forces.x, total_forces.y, total_forces.z, 0));
				_out_torque[particleIndex].w = (int)pos.w;
			}
		}

		//todo verificar funcao
		__global__ void CopyRBDataToParticleK(uint* _out_force_rb, uint* _out_force_index, float4* _in_rb_particle_forces, size_t num_particles)
		{
			unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

			if (idx <  num_particles)
			{
				_out_force_rb[idx] = (unsigned int)_in_rb_particle_forces[idx].w;
				_out_force_index[idx] = idx;
			}
		}

		//todo verificar funcao
		__global__ void ReorderForceTorqueK(float4* _out_force, float4* _out_torque, float4* _in_force, float4* _in_torque,
					uint *_in_ordered_index, uint _num_values)
		{
			extern __shared__ float4 s_forces[];

			unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

			if (idx < _num_values)
			{
				unsigned int idxvalue = _in_ordered_index[idx];

				//s_forces[threadIdx.x] = _in_values[idxvalue];
				//__syncthreads();


				//_out_values[idx] = s_forces[threadIdx.x];
				_out_force[idx] =  _in_force[idxvalue];
				_out_torque[idx] = _in_torque[idxvalue];
			}

		}

		/*
		 * Metodo responsavel pela adicao de torque e forca no centro de massa
		 * do corpo rigido
		 */
		__global__ void AddForceTorqueK(float4* _in_sum_force, float4* _in_sum_torque, float4* _in_ordened_force,
								  float4* _out_to_force, float4* _out_to_torque, int* _in_sum_offset, size_t num_rigid_body)
		{
			unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;

			if (idx <  num_rigid_body)
			{
				int idxoffset = _in_sum_offset[idx];
				int idxrb = (int) _in_ordened_force[idxoffset].w;
				float4 force_amount = _in_sum_force[idxoffset];
				float4 torque_amount = _in_sum_torque[idxoffset];

				_out_to_force[idxrb] += force_amount;
				_out_to_torque[idxrb] += torque_amount;
			}
		}

		/*
		 * Metodo utilizado para integracao temporal das posicoes de corpos rigidos
		 */
		__global__ void IntegrateRBK(float4* _in_RBPosition, float* _in_mass,
			float4* _in_RBVelocities, float4* _in_rbforces, size_t num_rb,
			float _elapsed_time_in_seconds, float4 _in_external_force,
			float3 min_volume, float3 max_volume, float4* _angular_momentum,
			float4* _torque, float4* _inertia_tensor, float4* _orientation,
			float4* _inv_inertia_tensor, float velocity_limit)
		{

			int idx = (blockIdx.x * blockDim.x) + threadIdx.x;

			if (idx < num_rb)
			{
				float4 force_amount = _in_rbforces[idx];

				float spd = length(force_amount);

				if (spd > velocity_limit)
				{
					float t = velocity_limit / spd;
					force_amount *= t;
				}

				float4 pos = _in_RBPosition[idx];
				float mass = _in_mass[idx];

				float4 accel = ((force_amount + _in_external_force) / mass) * _elapsed_time_in_seconds;
				float4 vel = _in_RBVelocities[idx] + (accel * _elapsed_time_in_seconds);

				pos.x += vel.x * _elapsed_time_in_seconds;
				pos.y += vel.y * _elapsed_time_in_seconds;
				pos.z += vel.z * _elapsed_time_in_seconds;

				//// Calcular momento angular
				_angular_momentum[idx] += _torque[idx] * _elapsed_time_in_seconds;
				_angular_momentum[idx] *= 0.9;

				//// Calcular velocidade angular
				float4 angularVelocity;
				angularVelocity = SimulationLib::MathUtils::AngularVelocityD(
						&_inv_inertia_tensor[idx * 3], _angular_momentum[idx],
						_orientation[idx]);

				angularVelocity.w = 0;
				angularVelocity =SimulationLib::MathUtils::Quaternion_ProductD(
					angularVelocity, _orientation[idx]) * 0.5;


				_orientation[idx] += SimulationLib::MathUtils::Quaternion_ProductD(
					angularVelocity * _elapsed_time_in_seconds * 0.1f, _orientation[idx]);
				_orientation[idx] = SimulationLib::MathUtils::Quaternion_NormalizeD(
					_orientation[idx]);

				//// Derivar orientação
//				float4 star[3];
//				float4 orientIntegration[3];
//				float4 newOrientation[3];
//				SimulationLib::MathUtils::Star(angularVelocity, star);
//
//				SimulationLib::MathUtils::MatrixMultD(star, quaternionToMatrix,
//						orientIntegration);
//				orientIntegration[0] *= _elapsed_time_in_seconds;
//				orientIntegration[1] *= _elapsed_time_in_seconds;
//				orientIntegration[2] *= _elapsed_time_in_seconds;
//				SimulationLib::MathUtils::MatrixAddD(orientIntegration,
//						quaternionToMatrix, newOrientation);
//
//				_orientation[idx] = SimulationLib::MathUtils::MatrixToQuaternionD(
//						newOrientation);


				if (pos.y < min_volume.y)
				{
					pos.y = min_volume.y;
					vel.y *= -0.7f;
				}

				if (pos.x <= min_volume.x)
				{
					pos.x = min_volume.x;
					vel.x *= -0.7f;
				}
				else if (pos.x >= max_volume.x)
				{
					pos.x = max_volume.x;
					vel.x *= -0.7f;
				}

				if (pos.z <= min_volume.z)
				{
					pos.z = min_volume.z;
					vel.z *= -0.7f;
				}
				else if (pos.z >= max_volume.z)
				{
					pos.z = max_volume.z;
					vel.z *= -0.7f;
				}

				_in_RBPosition[idx] = pos;
				_in_RBVelocities[idx] = vel;

			}
		}



	}
}

/*
 * Funcao utilizada para o calculo de forcas geradas na colisao entre corpos rigidos
 * estaticos e dinamicos
 */
extern "C"
	void gCalculateRigidBodyForcesFromRB( float4*       _in_RBPositions,    // Posições de cada RB partícula
									float4*       _in_RBVelocites,
									float4*       _in_RBCMPosition,
									float4*       _out_torque,
									float4*       _out_Forces,        // Total das forças em cada partícula
									unsigned int* _in_RBGridStart,    // Offset inicial de cada célula (corpos dinamicos)
									unsigned int* _in_RBGridEnd,      // Offset final de cada célula (corpos dinamicos)
									unsigned int  _in_NumParticles,   // Número total de partículas
									float3        _in_CellSize,       // Tamanho de cada célula
									float         _in_Radius,         // Raio de cada partícula
									int           _in_RBHashTableSize,// Total de Buckets para RB
									int           _in_SRBHashTableSize, // Total de Buckets para SRB
									float4*       _in_SRBPositions,   // Posição de cada SRB partícula
									unsigned int* _in_SRBGridStart,   // Offset inicial de cada célula (corpos estáticos)
									unsigned int* _in_SRBGridEnd,     // Offset final de cada célula (corpos estáticos)
									float         _in_RBkSpring,      // Constante de Spring
									float         _in_RBkDamping,    // Constante de damping
									int           _in_numSRB
									)
	{
		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)_in_NumParticles / blockSize.x)), 1, 1);
		size_t smem = sizeof(unsigned int) * blockSize.x;

		SimulationLib::RigidBody::CalculateRigidBodyForcesFromRBK<<<gridSize, blockSize, smem>>>(
			_in_RBPositions, _in_RBVelocites, _in_RBCMPosition, _out_torque, _out_Forces,
			_in_RBGridStart, _in_RBGridEnd, _in_NumParticles, _in_CellSize, _in_Radius,
			_in_RBHashTableSize, _in_SRBHashTableSize,  _in_SRBPositions,
			_in_SRBGridStart, _in_SRBGridEnd, _in_RBkSpring, _in_RBkDamping, _in_numSRB);
	}


/*
 * Funcao utilizada para o calculo das forcas dos corpos rigidos nos corpos rigidos
 */
extern "C"
	void gCalculateRigidBodyForcesFromFluids(float4*       _in_RBPositions,
											float4*       _in_RBVelocites,
											float4*       _in_RBCMPos,
											float4*       _out_Forces,        // Total das forças em cada partícula
											float4*       _out_torque,
											unsigned int* _in_RBGridStart,    // Offset inicial de cada célula (corpos dinamicos)
											unsigned int* _in_RBGridEnd,      // Offset final de cada célula (corpos dinamicos)
											unsigned int  _in_NumParticles,   // Número total de partículas
											float3        _in_CellSize,       // Tamanho de cada célula
											float         _in_Radius,         // Raio de cada partícula
											int           _in_RBHashTableSize,// Total de Buckets para RB

											float4*      _in_sph_particles,
											float4*      _in_sph_velocities,
											uint*        _in_sph_grid_start,
											uint*        _in_sph_grid_end,
											float        _in_sph_particle_radius,
											float3       _in_sph_cell_size,
											uint         _in_sph_num_buckets,
											int3         _in_cell_size_relation,
						    				float        _in_sph_kSpring,
						    				float        _in_sph_kDamping)
	{
		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)_in_NumParticles / blockSize.x)), 1, 1);

		SimulationLib::RigidBody::CalculateRBForcesFromFluidsK<<<gridSize, blockSize>>>(
			_in_RBPositions, _in_RBVelocites, _in_RBCMPos, _out_Forces, _out_torque,
			_in_RBGridStart, _in_RBGridEnd, _in_NumParticles, _in_CellSize, _in_Radius,
			_in_RBHashTableSize, _in_sph_particles, _in_sph_velocities, _in_sph_grid_start,
			_in_sph_grid_end, _in_sph_particle_radius, _in_sph_cell_size, _in_sph_num_buckets,
			_in_cell_size_relation, _in_sph_kSpring, _in_sph_kDamping);
	}

extern "C"
	void gCopyRBDataToParticle(uint* _out_force_rb, uint* _out_force_index, float4* _in_rb_particle_forces, size_t num_particles)
	{
		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)num_particles / blockSize.x)), 1, 1);

		SimulationLib::RigidBody::CopyRBDataToParticleK<<<gridSize, blockSize>>>(
			_out_force_rb,  _out_force_index, _in_rb_particle_forces, num_particles);
	}

extern "C"
	void gReorderForceTorque(float4* _out_force, float4* _out_torque, float4* _in_force, float4* _in_torque,
			uint *_in_ordered_index, uint _num_values)
	{
		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)_num_values / blockSize.x)), 1, 1);

		size_t smem = sizeof(float4) * blockSize.x;

		SimulationLib::RigidBody::ReorderForceTorqueK<<<gridSize, blockSize, smem>>>(
			_out_force, _out_torque, _in_force, _in_torque, _in_ordered_index,
			_num_values);
	}

extern "C"
	void gAddForceTorque(float4* _in_sum_force, float4* _in_sum_torque, float4* _in_ordened_force,
		float4* _out_to_force, float4* _out_to_torque, int* _in_sum_offset, size_t num_rigid_body)
	{
		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)num_rigid_body / blockSize.x)), 1, 1);

		SimulationLib::RigidBody::AddForceTorqueK<<<gridSize, blockSize>>>(_in_sum_force,  _in_sum_torque, _in_ordened_force,
			_out_to_force, _out_to_torque, _in_sum_offset, num_rigid_body);
	}



extern "C"
	void gIntegrateRB(float4* _in_RBPosition, float* _in_mass, float4* _in_RBVelocities, float4* _in_rbforces, size_t num_rb,
		float _elapsed_time_in_seconds, float4 _in_external_force, float3 min_volume, float3 max_volume, float4* _angular_momentum,
		float4* _torque, float4* _inertia_tensor, float4* _orientation, float4* _inv_inertia_tensor, float velocity_limit)
	{
		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)num_rb / blockSize.x)), 1, 1);

		SimulationLib::RigidBody::IntegrateRBK<<<gridSize, blockSize>>>(_in_RBPosition, _in_mass, _in_RBVelocities, _in_rbforces, num_rb,
			_elapsed_time_in_seconds, _in_external_force, min_volume, max_volume, _angular_momentum, _torque, _inertia_tensor,
			_orientation, _inv_inertia_tensor, velocity_limit);
	}

#endif
