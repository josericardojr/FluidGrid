/*
 * RigidBodySimulator.h
 *
 *  Created on: 03/02/2011
 *      Author: josericardodasilvajunior
 */

#ifndef RIGIDBODYSIMULATOR_H_
#define RIGIDBODYSIMULATOR_H_

#include <cudpp.h>
#include <cutil_inline.h>
#include "RigidBody.h"
#include "ApplicationData.cuh"
#include <vector>
#include "jrfxgl.h"
//#include "AuxFunctions.h"

typedef unsigned int uint;

namespace JRFXGL
{
	namespace RigidBodySim
	{

		/*
		 * Estrutura utilizada para o armazenamento das informacoes de colisao
		 */
		struct CollisionData
		{
			/* Informacoes de Colisao em CPU */
			float4*                h_CollidedParticleRelativePos;
			uint*                  h_CollidedParticleFlag;
			int*                   h_CollidedParticleStartOffset;

			/* Informacoes de colisao em GPU */
			uint*              d_CollidedParticleFlag;
			int*               d_CollidedParticleStartOffset;
			uint*              d_CollidedParticlesHashKey;
			uint*              d_CollidedParticleId;
			float4*            d_CollidedParticlesPos;
			float4*            d_CollidedParticlesForces;
			float4*            d_CollidedParticleTorque;
			uint*              d_CollidedBucketStart;
			uint*              d_CollidedBucketEnd;
			uint*              d_CollidedHashTable;

			CollisionData()
			{
				h_CollidedParticleRelativePos = NULL;
				h_CollidedParticleFlag = NULL;
				h_CollidedParticleStartOffset = NULL;

				d_CollidedParticleFlag = NULL;
				d_CollidedParticleStartOffset = NULL;
				d_CollidedParticlesHashKey = NULL;
				d_CollidedParticleId = NULL;
				d_CollidedParticlesPos = NULL;
				d_CollidedBucketStart = NULL;
				d_CollidedBucketEnd = NULL;
				d_CollidedHashTable = NULL;
				d_CollidedParticlesForces = NULL;
				d_CollidedParticleTorque = NULL;
			}

			void initialize(int numParticles, int numEntities, int numBuckets)
			{
				// Armazena o array de partículas que colidiram com o corpo rígido
				h_CollidedParticleRelativePos = new float4[numParticles];

				// Array de flags de partículas de corpos rígidos com colisão. Utilizando na etapa de soma prefixada em gpu
				h_CollidedParticleFlag = new uint[numParticles];

				// Offset com a localização da primeira partícula de cada corpo rígido que colidiu com corpo rígido
				h_CollidedParticleStartOffset = new int[numEntities];

				// Memoria de GPU
				cutilSafeCall(cudaMalloc((void**)&d_CollidedParticleFlag,
						sizeof(uint) * numParticles));
				cutilSafeCall(cudaMalloc((void**)&d_CollidedParticleStartOffset,
						sizeof(int) * numEntities));
				cutilSafeCall(cudaMalloc((void**)&d_CollidedParticlesHashKey,
						sizeof(uint) * numParticles));
				cutilSafeCall(cudaMalloc((void**)&d_CollidedParticleId,
						sizeof(uint) * numParticles));
				cutilSafeCall(cudaMalloc((void**)&d_CollidedParticlesPos,
						sizeof(float4) * numParticles));
				cutilSafeCall(cudaMalloc((void**)&d_CollidedParticlesForces,
						sizeof(float4) * numParticles));
				cutilSafeCall(cudaMalloc((void**)&d_CollidedBucketStart,
						sizeof(uint) * numBuckets));
				cutilSafeCall(cudaMalloc((void**)&d_CollidedBucketEnd,
						sizeof(uint) * numBuckets));
				cutilSafeCall(cudaMalloc((void**)&d_CollidedHashTable,
						sizeof(uint) * numBuckets));
				cutilSafeCall(cudaMalloc((void**)&d_CollidedParticleTorque,
						sizeof(float4) * numParticles));

			}

			void clear()
			{
				if (h_CollidedParticleRelativePos != NULL) free(h_CollidedParticleRelativePos);
				if (h_CollidedParticleFlag != NULL) free(h_CollidedParticleFlag);
				if (h_CollidedParticleStartOffset != NULL) free(h_CollidedParticleStartOffset);

				if (d_CollidedParticleFlag != NULL)
					cutilSafeCall(cudaFree(d_CollidedParticleFlag));
				if (d_CollidedParticleStartOffset != NULL)
					cutilSafeCall(cudaFree(d_CollidedParticleStartOffset));
				if (d_CollidedParticlesHashKey != NULL)
					cutilSafeCall(cudaFree(d_CollidedParticlesHashKey));
				if (d_CollidedParticleId != NULL)
					cutilSafeCall(cudaFree(d_CollidedParticleId));
				if (d_CollidedParticlesPos != NULL)
					cutilSafeCall(cudaFree(d_CollidedParticlesPos));
				if (d_CollidedParticlesForces != NULL)
					cutilSafeCall(cudaFree(d_CollidedParticlesForces));
				if (d_CollidedBucketEnd != NULL)
					cutilSafeCall(cudaFree(d_CollidedBucketEnd));
				if (d_CollidedBucketStart != NULL)
					cutilSafeCall(cudaFree(d_CollidedBucketStart));
				if (d_CollidedHashTable != NULL)
					cutilSafeCall(cudaFree(d_CollidedHashTable));
				if (d_CollidedParticleTorque != NULL)
					cutilSafeCall(cudaFree(d_CollidedParticleTorque));
			}
		};

		class RigidBodySimulator
		{
		private:
			RigidBodySimulator();

		public:
			static RigidBodySimulator* instance();

			void addElement(RigidBody* rbElement);
			void reset();
			void release();

			/* Metodos acessores para recuperar informacoes dos corpos rigidos */
			unsigned int* getStaticRBGridStartOffset(){ return m_dStaticRBGridOffsetStart; }
			unsigned int* getStaticRBGridEndOffset(){ return m_dStaticRBGridOffsetEnd; }
			float4*       getStaticRBParticlePositions(){ return m_dStaticRBRelativeParticlePositions; }
			int           getNumStaticParticles(){ return mStaticParticlePositions.size(); }

			//void generateCollisionData(vl::AABB fluidBB);
			void generateCollisionData();
			void startInteraction();
			void endInteraction();

			void calculateRigidBodyForcesFromRB();
			void calculateRigidBodyForcesFromFluid(float4* dfluidParticlePos,
					float4* dfluidParticleVel, uint* dfluidHashStart, uint* dfluidHashEnd);
			void applyForces();
			void integrate();

			// Recuperar informacao do numero de colisoes que ocorreu na ultima atualizacao
			int getNumCollisionsWithFluid(){ return mNumCollisionRBFluid; }

			int getNumStaticRigidBodies(){ return mStaticEntities.size(); }

			CollisionData* getFluidCollisionData(){ return &mFluidCollisionData; }
			float4* getDynamicVelocities(){ return m_dDynamicRBCenterMassVelocity; }

		private:
			void clearGPUData();
			void createHashTable();

			void updateStaticElements();
			void updateDynamicRigidBody();


		private:
			/* Corpos rigidos estaticos */
			std::vector<RigidBody*>     mStaticEntities;
			std::vector<float4>         mStaticParticlePositions;

			/* Corpos rigidos dinamicos */
			std::vector<RigidBody*>     mDynamicEntities;
			std::vector<uint>      mDynamicNumParticles;    // Armazena o número de partícula de cada corpo rígido
			std::vector<float4>    mDynamicParticlePosRel;  // Armazena a posição relativa de cada partícula do RB (só atualiza a GPU 1 vez)
			std::vector<uint>      mDynamicParticleEndOffset; // Armazena o offset da partícula inicial de cada RB
			std::vector<uint>      mDynamicParticleStartFlag;    // Armazena os flags para utilização na soma em GPU
			std::vector<float4>    mDynamicOrientation;     // Orientação do corpo rígido em CPU
			std::vector<float4>    mDynamicCenterMass;      // Posicao do centro de massa do corpo rigido
			std::vector<float>     mDynamicMass;            // Massa do corpo rigido

			uint2*                 mRBCollisions;          // Flags que indicam o corpo rigido que colidiu e com quem colidiu
			float4*                mRBPositions;
			float4*                mRBOrientations;
			std::vector<sInertiaTensor> mInertiaTensor;
			std::vector<sInertiaTensor> mInvInertiaTensor;

			CollisionData          mRBCollisionData;
			CollisionData          mFluidCollisionData;

			int                    mNumCollisionRBRB;      // Armazena o total de corpos rigidos que colidiram
			int                    mNumCollisionRBFluid;   // Armazena o total de corpos rigidos que colidiram com o fluido
			int                    mNumParticlesRBRBCol;   // Armazena o total de particulas de corpos rigidos que colidiram
			int                    mNumParticlesRBFluidCol; // Armazena o total de particulas de corpos rigidos que colidiram com fluido


			// GPU Data
			/* Corpos rigidos estaticos */
			float4*                     m_dStaticRBRelativeParticlePositions;
			unsigned int*               m_dStaticRBGridOffsetStart;
			unsigned int*               m_dStaticRBGridOffsetEnd;

			/* Corpos rigidos dinamicos */
			CUDPPHandle                 mSortHandle;
			CUDPPHandle                 mScanHandle;
			float4*	           m_dDynamicRBCenterMassPositions;
			float4*	           m_dDynamicParticleRelPosition;
			float4*            m_dDynamicRBCenterMassVelocity;
			float4*            m_dDynamicRBForces;
			float4*            m_dDynamicRBTorque;
			float4*            m_dDynamicRBOrientation;
			float*             m_dDynamicRBMass;
			float4*            m_dDynamicInertiaTensor;
			float4*            m_dDynamicInvInertiaTensor;
			float4*            m_dDynamicAngularMomentum;
			float4*            m_dDynamicParticlesSizedTmp;
			float4*            m_dDynamicParticlesSizedTmp2;


			static RigidBodySimulator* mInstance;
			bool mInteractionStarted;
		};
	}
}

#endif /* RIGIDBODYSIMULATOR_H_ */
