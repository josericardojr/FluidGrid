#ifndef __K_HashGrid_CUH__

#define __K_HashGrid_CUH__

extern "C"
	void 
		gCalcHashKey(unsigned int* hashKey, unsigned int* particleID, float4 *positions, 
			unsigned int numParticles, unsigned int numBuckets, float3 cell_size);
	
extern "C"
	void gCalcHashKeyWithOffset(unsigned int* _out_hashKey,    // Vetor de hash key para cada elemento
								unsigned int* _out_particleID, // Vetor que conter� o ID de cada posi��o
								float4*       _in_vectorOffset,// Vetor com os offsets a serem adicionados
								float4*       _in_orientation, // Orientacao
								float4*       _in_vectorPos,   // Vetor com as posi��es
								float4*       _out_positions,  // Armazena o resultado da posi��o com o offset
								unsigned int  _in_numParticles,// N�mero de part�culas
								float3        _in_cellSize,    // Tamanho de cada c�lula da grid regular
								int           _in_hashTableSize// N�mero de buckets
								);					  
			
			
extern "C"
	void gFindCellOffsetOrderingPosVel( unsigned int* gridPosStart,             // Vetor com o in�cio (offset) de cada c�lula
		unsigned int* gridPosEnd,               // Vetor com a posi��o final de cada c�lula
		unsigned int* hashTable,                // Vetor com a tabela de hash keys
		unsigned int* particleID,               // Vetor com o ID de cada part�cula
		float4*       particlePositions,        // Vetor de posi��es
		float4*       orderedParticlePositions, // Vetor de posi��es ordenadas
		float4*       particleVel,              // Velocidade de cada part�cula
		float4*       particleOrderedVel,       // Velocidade ordenada de cada part�cula
		unsigned int  numParticles,             // Total de part�culas
		unsigned int  numBuckets
		);
			

extern "C"
	void gFindCellOffsetOrderingPos( unsigned int* gridPosStart,             // Vetor com o in�cio (offset) de cada c�lula
						 unsigned int* gridPosEnd,               // Vetor com a posi��o final de cada c�lula
						 unsigned int* hashTable,                // Vetor com a tabela de hash keys
						 unsigned int* particleID,               // Vetor com o ID de cada part�cula
						 float4*       particlePositions,        // Vetor de posi��es
						 float4*       orderedParticlePositions, // Vetor de posi��es ordenadas
						 unsigned int  numParticles,             // Total de part�culas
						 unsigned int  numBuckets                // N�mero de buckets
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