#ifndef __K_HASHGRID_CU__

#define __K_HASHGRID_CU__

#include "K_DeviceFunctions.cu"
#include "K_SimulationLib.cuh"
#include "K_MathUtils.cu"
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_math.h>
#include <cutil_inline.h>

#if USE_TEX
texture<float4, 1, cudaReadModeElementType> oldPosTex;
texture<float4, 1, cudaReadModeElementType> oldVelTex;
#endif



namespace SimulationLib {
	namespace DataStructure {

		/*
		 * Calcular codigo de hash
		 */
		__global__ void CalcHashKeyK(unsigned int* hashKey, unsigned int* particleID, float4* position, 
				unsigned int numParticles, float3 cellSize, int hashTableSize)
		{
			int idx = blockDim.x * blockIdx.x + threadIdx.x;
		
			if (idx < numParticles)
			{	
				float4 pos = position[idx];
				particleID[idx] = idx;
				hashKey[idx] = SimulationLib::Common::HashKeyD(
						make_float3(pos.x, pos.y, pos.z), cellSize, hashTableSize);
			}
		}
		
		/*
		 * Calculo do codigo de hash com atualizacao de posicao atraves de um offset
		 */
		__global__ void CalcHashKeyWithOffsetK( unsigned int* _out_hashKey,
			unsigned int* _out_particleID, float4* _in_vectorOffset, float4* _in_orientation,
			float4* _in_vectorPos, float4* _out_positions, unsigned int  _in_numParticles,
			float3 _in_cellSize, int _in_hashTableSize )
		{
			int idx = blockDim.x * blockIdx.x + threadIdx.x;

			if (idx < _in_numParticles)
			{
				float4 pos = _in_vectorPos[idx];
				int offsetIndex = (int) pos.w;

				float4 absPos = SimulationLib::MathUtils::Quaternion_TransformPointD(pos,
					_in_orientation[offsetIndex]) + _in_vectorOffset[offsetIndex];
				absPos.w = pos.w;
				_out_particleID[idx] = idx;
				_out_positions[idx] = absPos;
				_out_hashKey[idx] = SimulationLib::Common::HashKeyD(
						make_float3(absPos.x, absPos.y, absPos.z), _in_cellSize, _in_hashTableSize);
			}
		}


		/*
		 * Kernel para calculo do offset das celulas, ordenando posicao e velocidade,
		 * garantindo acesso coalescido
		 */
		__global__ void FindCellOffsetOrderinPosVelK(unsigned int*   gridPosStart,
													 unsigned int*   gridPosEnd,
													 unsigned int*   hashTable,
													 unsigned int*   particleID,
													 float4*         oldPos,
													 float4*         particlePositionsOrdered,
													 float4*         oldVel,
													 float4*         particleOrderedVel,
													 unsigned int    numParticles)

		{
			extern __shared__ unsigned int sharedHash[];    // blockSize + 1 elements
			unsigned int index = __umul24(blockIdx.x,blockDim.x) + threadIdx.x;

			unsigned int hash;

			// handle case when no. of particles not multiple of block size
			if (index < numParticles)
			{
				hash = hashTable[index];

				// Load hash data into shared memory so that we can look
				// at neighboring particle's hash value without loading
				// two hash values per thread
				sharedHash[threadIdx.x+1] = hash;

				if (index > 0 && threadIdx.x == 0)
				{
					// first thread in block must load neighbor particle hash
					sharedHash[0] = hashTable[index-1];
				}
			}

			__syncthreads();

			if (index < numParticles)
			{
				// If this particle has a different cell index to the previous
				// particle then it must be the first particle in the cell,
				// so store the index of this particle in the cell.
				// As it isn't the first particle, it must also be the cell end of
				// the previous particle's cell

			    if (index == 0 || hash != sharedHash[threadIdx.x])
			    {
				    gridPosStart[hash] = index;

					if (index > 0)
						gridPosEnd[sharedHash[threadIdx.x]] = index;
			    }

				if (index == numParticles - 1)
					gridPosEnd[hash] = index + 1;

			    // Now use the sorted index to reorder the pos and vel data
			    unsigned int sortedIndex = particleID[index];
			    float4 pos = FETCH(oldPos, sortedIndex);
			    float4 vel = FETCH(oldVel, sortedIndex);

			    particlePositionsOrdered[index] = pos;
			    particleOrderedVel[index] = vel;
			}

		}


		/*
		 * Kernel para calculo do offset das celulas, ordenando posicao,
		 * garantindo acesso coalescido
		 */
		__global__ void FindCellOffsetOrderingPosK(unsigned int*   gridPosStart,
																		unsigned int*   gridPosEnd,
																		unsigned int*   hashTable,
																		unsigned int*   particleID,
																		float4*         oldPos,
																		float4*         particlePositionsOrdered,
																		unsigned int    numParticles)
		{
			extern __shared__ unsigned int sharedHash[];    // blockSize + 1 elements
			unsigned int index = (blockIdx.x * blockDim.x) + threadIdx.x;

			unsigned int hash;

			// handle case when no. of particles not multiple of block size
			if (index < numParticles)
			{
				hash = hashTable[index];

				// Load hash data into shared memory so that we can look
				// at neighboring particle's hash value without loading
				// two hash values per thread
				sharedHash[threadIdx.x+1] = hash;

				if (index > 0 && threadIdx.x == 0)
				{
					// first thread in block must load neighbor particle hash
					sharedHash[0] = hashTable[index-1];
				}
			}

			__syncthreads();

			if (index < numParticles)
			{
				// If this particle has a different cell index to the previous
				// particle then it must be the first particle in the cell,
				// so store the index of this particle in the cell.
				// As it isn't the first particle, it must also be the cell end of
				// the previous particle's cell

			    if (index == 0 || hash != sharedHash[threadIdx.x])
			    {
				    gridPosStart[hash] = index;

					if (index > 0)
						gridPosEnd[sharedHash[threadIdx.x]] = index;
			    }

				if (index == numParticles - 1)
					gridPosEnd[hash] = index + 1;

			    // Now use the sorted index to reorder the pos and vel data
			    unsigned int sortedIndex = particleID[index];
			    float4 pos = FETCH(oldPos, sortedIndex);       // macro does either global read or texture fetch

			    particlePositionsOrdered[index] = pos;
			}
		}



	}
}

/* 
 * Metodo externo para calcular codigo de hash
 */
extern "C" 
   void gCalcHashKey(unsigned int* hashKey, unsigned int* particleID, float4 *positions,
		   unsigned int numParticles, unsigned int numBuckets, float3 cell_size)
   {
   		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)numParticles / blockSize.x)), 1, 1);
		
		SimulationLib::DataStructure::CalcHashKeyK<<<gridSize, blockSize>>>(hashKey, particleID, positions,
				numParticles, cell_size, numBuckets);
   }

/*
 * Metodo externdo para calculo de codigo de hash utilizando um deslocamento de posicao
 */
extern "C"
	void gCalcHashKeyWithOffset( unsigned int* _out_hashKey, unsigned int* _out_particleID,
			float4* _in_vectorOffset, float4* _in_orientation,
			float4* _in_vectorPos, float4* _out_positions, unsigned int  _in_numParticles,
			float3 _in_cellSize, int _in_hashTableSize)
	{
		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)_in_numParticles / blockSize.x)), 1, 1);

		SimulationLib::DataStructure::CalcHashKeyWithOffsetK<<<gridSize, blockSize>>>(_out_hashKey, _out_particleID,
			_in_vectorOffset, _in_orientation, _in_vectorPos, _out_positions, _in_numParticles,
			_in_cellSize, _in_hashTableSize);
	}

/*
 * Funcao global para ordenacao de posicao e velocidade, gerando uma tabela de hash
 */
extern "C"
	void gFindCellOffsetOrderingPosVel(unsigned int* gridPosStart,
									   unsigned int* gridPosEnd,
									   unsigned int* hashTable,
									   unsigned int* particleID,
									   float4*       particlePositions,
									   float4*       orderedParticlePositions,
									   float4*       particleVel,
									   float4*       particleOrderedVel,
									   unsigned int  numParticles,
									   unsigned int  numBuckets)
	{
		dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)numParticles / blockSize.x)), 1, 1);

		// set all cells to empty
		cutilSafeCall(cudaMemset(gridPosStart, 0xffffffff,
				numBuckets*sizeof(unsigned int)));

		unsigned int smemSize = sizeof(unsigned int)*(blockSize.x + 1);

		#if USE_TEX
		// Mapear memoria de textura
	    cutilSafeCall(cudaBindTexture(0, oldPosTex, particlePositions, numParticles*sizeof(float4)));
	    cutilSafeCall(cudaBindTexture(0, oldVelTex, particleVel, numParticles*sizeof(float4)));
		#endif

		SimulationLib::DataStructure::FindCellOffsetOrderinPosVelK<<< gridSize, blockSize, smemSize>>>(
			gridPosStart,
			gridPosEnd,
			hashTable,
			particleID,
			particlePositions,
			orderedParticlePositions,
			particleVel,
			particleOrderedVel,
			numParticles);

		// Remover mapeamento
		#if USE_TEX
	    cutilSafeCall(cudaUnbindTexture(oldPosTex));
	    cutilSafeCall(cudaUnbindTexture(oldVelTex));
		#endif
	}


/*
 * Funcao global para ordenacao de posicao, gerando uma tabela de hash
 */
extern "C"
	void gFindCellOffsetOrderingPos(unsigned int* gridPosStart,
									unsigned int* gridPosEnd,
									unsigned int* hashTable,
									unsigned int* particleID,
									float4*       particlePositions,
									float4*       orderedParticlePositions,
									unsigned int  numParticles,
									unsigned int  numBuckets)
	{
    	dim3 blockSize(128, 1, 1);
		dim3 gridSize((int)(ceil((float)numParticles / blockSize.x)), 1, 1);

		// set all cells to empty
		cutilSafeCall(cudaMemset(gridPosStart, 0xffffffff, numBuckets*sizeof(unsigned int)));

		unsigned int smemSize = sizeof(unsigned int)*(128+1);

		#if USE_TEX
		// Mapear memoria de textura
		cutilSafeCall(cudaBindTexture(0, oldPosTex, particlePositions, numParticles*sizeof(float4)));
		#endif

		SimulationLib::DataStructure::FindCellOffsetOrderingPosK<<< gridSize, blockSize, smemSize>>>(
				gridPosStart,
				gridPosEnd,
				hashTable,
				particleID,
				particlePositions,
				orderedParticlePositions,
				numParticles);
	}




#endif
