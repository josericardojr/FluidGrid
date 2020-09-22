#ifndef __K_HashGrid_CUH__

#define __K_HashGrid_CUH__

extern "C"
	void 
		gCalcHashKey(unsigned int* hashKey, unsigned int* particleID, float4 *positions, 
			unsigned int numParticles, unsigned int numBuckets, float3 cell_size);
	
extern "C"
	void gCalcHashKeyWithOffset(unsigned int* _out_hashKey,    // Vetor de hash key para cada elemento
								unsigned int* _out_particleID, // Vetor que conterá o ID de cada posição
								float4*       _in_vectorOffset,// Vetor com os offsets a serem adicionados
								float4*       _in_orientation, // Orientacao
								float4*       _in_vectorPos,   // Vetor com as posições
								float4*       _out_positions,  // Armazena o resultado da posição com o offset
								unsigned int  _in_numParticles,// Número de partículas
								float3        _in_cellSize,    // Tamanho de cada célula da grid regular
								int           _in_hashTableSize// Número de buckets
								);					  
			
			
extern "C"
	void gFindCellOffsetOrderingPosVel( unsigned int* gridPosStart,             // Vetor com o início (offset) de cada célula
		unsigned int* gridPosEnd,               // Vetor com a posição final de cada célula
		unsigned int* hashTable,                // Vetor com a tabela de hash keys
		unsigned int* particleID,               // Vetor com o ID de cada partícula
		float4*       particlePositions,        // Vetor de posições
		float4*       orderedParticlePositions, // Vetor de posições ordenadas
		float4*       particleVel,              // Velocidade de cada partícula
		float4*       particleOrderedVel,       // Velocidade ordenada de cada partícula
		unsigned int  numParticles,             // Total de partículas
		unsigned int  numBuckets
		);
			

extern "C"
	void gFindCellOffsetOrderingPos( unsigned int* gridPosStart,             // Vetor com o início (offset) de cada célula
						 unsigned int* gridPosEnd,               // Vetor com a posição final de cada célula
						 unsigned int* hashTable,                // Vetor com a tabela de hash keys
						 unsigned int* particleID,               // Vetor com o ID de cada partícula
						 float4*       particlePositions,        // Vetor de posições
						 float4*       orderedParticlePositions, // Vetor de posições ordenadas
						 unsigned int  numParticles,             // Total de partículas
						 unsigned int  numBuckets                // Número de buckets
						 );
						 
extern "C"
	void gFindCellOffsetOrderingPos(unsigned int* gridPosStart,
									unsigned int* gridPosEnd,
									unsigned int* hashTable,
									unsigned int* particleID,
									float4*       particlePositions,
									float4*       orderedParticlePositions,
									unsigned int  numParticles,
									unsigned int  numBuckets);						 
		
#endif