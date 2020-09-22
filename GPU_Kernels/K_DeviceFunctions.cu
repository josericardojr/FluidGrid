#ifndef __K_DeviceFunctions_CU__

#define __K_DeviceFunctions_CU__

#include <vector_types.h>
#include <vector_functions.h>
#include "K_Boundaries_Common.cu"
#include <cutil_math.h>

namespace SimulationLib {
	namespace Common {
		
		/*
		 * Metodo utilizado para geracao de codigo de hash considerando o problema do 
		 * espelhamento do dominio
		 */
		__device__ inline unsigned int HashKeyD(float3 pos, float3 cellSize, int hashTableSize)
		{
			int x = (int) ((pos.x + 5000) / cellSize.x);
			int y = (int) ((pos.y + 5000) / cellSize.y);
			int z = (int) ((pos.z + 5000) / cellSize.z);
			
			return (unsigned int) ((x * 73856093) ^
				(y * 19349663) ^
				(z * 83492791)) % hashTableSize;
		}
		
		/*
		 * Metodo para verificar colisao entre uma esfera passada como argumento e as
		 * esferas pertencentes a uma celula
		 */
		__device__ inline float3 ProcessCellForcesD(
				/* Dados da esfera */
				float4        _in_sphere_position,
				float4        _in_sphere_velocity,
				float         _in_sphere_radius,

				/* Dados do conjunto de esferas a serem testadas */
				float         _in_sphere_collection_radius,
				float4*       _in_sphere_collection_positions,
				float4*       _in_sphere_collection_velocities,

				unsigned int  _in_hash_key,
				unsigned int* _in_hash_grid_start_offset,
				unsigned int* _in_hash_grid_end_offset,
				float         _in_kSpring,
				float         _in_kDamping,
				bool          _in_w_rb_index, // Indica se o valor em w deve ser usado para localizar
											  // informacoes no array de RB
				bool          _in_check_grouping, // Indica se e necessario checar se a particula
				                                  // pertence ao mesmo grupo das outras
				int           _in_w_group,
				float         _epsilon)

		{
			float3 _out_total_forces = make_float3(0.0, 0.0, 0.0);

			// Distancia minima
			float minimun_distance = _in_sphere_radius + _in_sphere_collection_radius + _epsilon;

			// Primeiramente deve ser verificado se a celula possui particulas
			if (_in_hash_grid_start_offset[_in_hash_key] != 0xffffffff)
			{
				// Percorrer todos as particulas nesta celula
				for (int i = _in_hash_grid_start_offset[_in_hash_key];
						i < _in_hash_grid_end_offset[_in_hash_key]; i++)
				{
					float4 pos = _in_sphere_collection_positions[i];
					float4 vel = make_float4(0.0, 0.0, 0.0, 0.0);

					if (_in_check_grouping && ( ((int) pos.w ) == _in_w_group)) continue;

					// Verificar se o array de distancia foi fornecido (rb estatico nao possuem
					// velocidade)
					if (_in_sphere_collection_velocities != NULL)
						if (_in_w_rb_index)
							vel = _in_sphere_collection_velocities[(int)pos.w];
						else
							vel = _in_sphere_collection_velocities[i];

					// Verificar distância
					float3 relative_distance = make_float3(pos) - make_float3(_in_sphere_position);
					float lg = length(relative_distance);

					if (lg < minimun_distance)
					{
						// Vetor normalizado com a direcao da normal onde ocorreu a colisao
						float3 normalDir = relative_distance / lg;

						float3 relative_vel = make_float3(vel) - make_float3(_in_sphere_velocity);


					/*	float3 fSpring = -_in_kSpring * (minimun_distance - lenght) * nDir;
						float3 fDamping = make_float3( _in_kDamping * relative_vel.x, _in_kDamping * relative_vel.y, _in_kDamping * relative_vel.z);
						float3 fTangent = 2 * (relative_vel -
							((relative_vel.x * nDir.x +relative_vel.y * nDir.y +relative_vel.z * nDir.z) * nDir));*/

						_out_total_forces += calculateRepulsionForce(relative_vel, normalDir, lg,
								_in_kDamping, _in_kSpring);

					}
				}
			}

			return _out_total_forces;
		}


		/*
		 * Metodo que processa forcas decorrente de colisao entre uma particula com um
		 * conjunto de particulas em uma celula e suas vizinhas
		 */
		__device__ inline float3 ProcessCollisionForcesD(
				/* Informacoes das particulas */
				float4        _in_sphere_position,
				float4        _in_sphere_velocity,
				float         _in_sphere_radius,

				/* Informacoes da grid que contem as particulas a serem verificadas */
				float4*       _in_sphere_collection_positions,
				float4*       _in_sphere_collection_velocities,
				float         _in_sphere_collection_radius,
				float         _in_kSpring,
				float         _in_kDamping,
				int3          _in_cell_size_relation,

				unsigned int* _in_hash_grid_start_offset,
				unsigned int* _in_hash_grid_end_offset,
				float3         _in_sphere_collection_cell_size,
				unsigned int  _in_sphere_collection_num_buckets,
				bool          _in_w_rb_index,
				bool          _in_check_grouping,
				int           _in_w_group,
				float         _epsilon)
		{

			float3 _out_force = make_float3(0.0f);

			for (int z = -_in_cell_size_relation.z; z <= _in_cell_size_relation.z; z++)
			{
				for ( int y = -_in_cell_size_relation.y; y <= _in_cell_size_relation.y; y++)
				{
					for (int x = -_in_cell_size_relation.x; x <= _in_cell_size_relation.x; x++)
					{
						// Calcular a posicao de Hash
						float3 offset = _in_sphere_collection_cell_size;
						offset.x *= x; offset.y *= y; offset.z *= z;
						unsigned int hashCell = HashKeyD(make_float3(_in_sphere_position) + offset,
								_in_sphere_collection_cell_size, _in_sphere_collection_num_buckets );

						_out_force += ProcessCellForcesD(_in_sphere_position, _in_sphere_velocity, _in_sphere_radius, _in_sphere_collection_radius,
							 _in_sphere_collection_positions, _in_sphere_collection_velocities, hashCell, _in_hash_grid_start_offset, _in_hash_grid_end_offset,
							 _in_kSpring, _in_kDamping, _in_w_rb_index, _in_check_grouping, _in_w_group, _epsilon);
					}
				}
			}

			return _out_force;
		}
		
	}
}

#endif
