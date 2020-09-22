/*
 * RigidBodySimulator.cpp
 *
 *  Created on: 03/02/2011
 *      Author: josericardodasilvajunior
 */

#include <cutil_inline.h>
#include <vector_types.h>
#include <vector_functions.h>
#include "cudpp.h"
#include "../GPU_Kernels/K_HashGrid.cuh"
#include "../GPU_Kernels/K_RigidBodySim.cuh"
#include <omp.h>

#include "RigidBodySimulator.h"

void printGPUData(float4* gpuData, int num)
{
/*	float4* host = (float4*) malloc(sizeof(float4) * num);
	cutilSafeCall(cudaMemcpy(host, gpuData, sizeof(float4)*num, cudaMemcpyDeviceToHost));

	for (int i = 0; i < num; i++)
	{
		vl::Log::print("X " + vl::String::fromDouble(host[i].x) +
			" Y " + vl::String::fromDouble(host[i].y) +
			" Z " + vl::String::fromDouble(host[i].z) +
			" W " + vl::String::fromDouble(host[i].w) + "\n");
	}

	free(host);*/
}
//Todo Colocar a classificacao dos elementos das particulas de corpos rigidos dinamicos que colidiram ao mesmo tempo com a classificacao das particulas de fluido utilizando a FERMI

JRFXGL::RigidBodySim::RigidBodySimulator*
	JRFXGL::RigidBodySim::RigidBodySimulator::mInstance = NULL;

namespace JRFXGL
{
	namespace RigidBodySim
	{
		RigidBodySimulator* RigidBodySimulator::instance()
		{
			if (mInstance == NULL)
				mInstance = new RigidBodySimulator();

			return mInstance;
		}

		RigidBodySimulator::RigidBodySimulator()
		{
			m_dStaticRBRelativeParticlePositions = NULL;
			m_dStaticRBGridOffsetStart = NULL;
			m_dStaticRBGridOffsetEnd = NULL;

			m_dDynamicRBCenterMassPositions = NULL;
			m_dDynamicParticleRelPosition = NULL;
			m_dDynamicRBCenterMassVelocity = NULL;
			m_dDynamicRBOrientation = NULL;
			m_dDynamicRBMass = NULL;
			m_dDynamicParticlesSizedTmp = NULL;
			m_dDynamicParticlesSizedTmp2 = NULL;
			m_dDynamicRBForces = NULL;
			m_dDynamicRBTorque = NULL;
			m_dDynamicInertiaTensor = NULL;
			m_dDynamicInvInertiaTensor = NULL;
			m_dDynamicAngularMomentum = NULL;

			mRBCollisions = NULL;
			mRBPositions = NULL;
			mRBOrientations = NULL;

			mInteractionStarted = false;
		}

		void RigidBodySimulator::addElement(RigidBody* rbElement)
		{
			mc::math::Vec3f pos = rbElement->getPosition();
			float4* relativeParticlePositions = rbElement->getRelativeParticlePositions();

			if (rbElement -> isStatic())
			{
				mStaticEntities.push_back(rbElement);

				for (int i = 0; i < rbElement->getNumParticles(); i++)
				{
					mStaticParticlePositions.push_back(make_float4(
						relativeParticlePositions[i].x + pos.x,
						relativeParticlePositions[i].y + pos.y,
						relativeParticlePositions[i].z + pos.z, 0.0f));
				}
			}
			else
			{
				mc::math::Quatf& ori = rbElement->getOrientation();

				// Objeto a ser simulado
				mDynamicEntities.push_back(rbElement);

				// Número de partículas do corpo rigido
				mDynamicNumParticles.push_back(rbElement->getNumParticles());

				// Offset para localização da partícula final de um corpo rígido dinamico
				mDynamicParticleEndOffset.push_back(mDynamicParticlePosRel.size() +
						rbElement->getNumParticles() - 1);

				// Posição do centro de massa
				mDynamicCenterMass.push_back(make_float4(pos.x, pos.y, pos.z, 0));

				// Orientacao do corpo rígido
				mDynamicOrientation.push_back(make_float4(ori.x, ori.y, ori.z, ori.w));

				// Tensor de inércia
				mInertiaTensor.push_back(rbElement->getInertiaTensor());
				mInvInertiaTensor.push_back(rbElement->getInvInertiaTensor());

				// Massa
				mDynamicMass.push_back(rbElement->getMass());

				for (int i = 0; i < rbElement->getNumParticles(); i++)
				{
					(i == 0) ? mDynamicParticleStartFlag.push_back(1) :
							mDynamicParticleStartFlag.push_back(0);

					// Campo w da posicao relativa armazena o indice do corpo rigido no vetor mDynamicCenterMass
					// da qual esta particula pertence, a fim de calcular a posicao absoluta posteriormente
					mDynamicParticlePosRel.push_back(make_float4(relativeParticlePositions[i].x,
						relativeParticlePositions[i].y, relativeParticlePositions[i].z, mDynamicCenterMass.size()-1));

				}
			}
		}

		void RigidBodySimulator::createHashTable()
		{
			Simulation::Application::CommonRigidBodyParams& commonRBParameters =
					Simulation::Application::ApplicationData::instance()->getCommonRigidBodyParams();

			// hash table dos corpos rigidos estaticos
			cutilSafeCall(cudaMalloc((void**) &m_dStaticRBGridOffsetStart,
					sizeof(unsigned int) * commonRBParameters.getStNumBuckets()));
			cutilSafeCall(cudaMalloc((void**) &m_dStaticRBGridOffsetEnd,
					sizeof(unsigned int) * commonRBParameters.getStNumBuckets()));
		}

		void RigidBodySimulator::updateStaticElements()
		{
			// Verificar se existem corpos rigidos estaticos na simulacao
			if (mStaticEntities.size() > 0)
			{
				Simulation::Application::CommonRigidBodyParams& commonRBParameters =
						Simulation::Application::ApplicationData::instance()->getCommonRigidBodyParams();


				// Alocar memoria GPU para elementos estaticos
				cutilSafeCall(
					cudaMalloc((void**) &m_dStaticRBRelativeParticlePositions,
							sizeof(float4) * mStaticParticlePositions.size()));
				// Copiar informacoes
				cutilSafeCall(cudaMemcpy(m_dStaticRBRelativeParticlePositions, &mStaticParticlePositions[0],
					sizeof(float4) * mStaticParticlePositions.size(), cudaMemcpyHostToDevice));


				// Memoria temporaria utilizada para classificacao na tabela de hash
				unsigned int* dParticleID;
				unsigned int* dParticleHashKey;
				float4* dParticleOrderedPos;

				cutilSafeCall(
					cudaMalloc((void**) &dParticleOrderedPos,
							sizeof(float4) * mStaticParticlePositions.size()));

				cutilSafeCall(cudaMalloc((void**) &dParticleID,
						sizeof(unsigned int) * mStaticParticlePositions.size()));

				cutilSafeCall( cudaMalloc((void**) &dParticleHashKey,
						sizeof(unsigned int) * mStaticParticlePositions.size()));

				// Classificacao
				CUDPPConfiguration config;
				config.algorithm = CUDPP_SORT_RADIX;
				config.datatype = CUDPP_UINT;
				config.op = CUDPP_MAX;
				config.options = CUDPP_OPTION_KEY_VALUE_PAIRS;

				CUDPPHandle handle = 0;
				CUDPPResult res = cudppPlan(&handle, config, mStaticParticlePositions.size(), 1, 0);

				//todo assert(handle);
				gCalcHashKey(dParticleHashKey, dParticleID, m_dStaticRBRelativeParticlePositions,
						mStaticParticlePositions.size(), commonRBParameters.getStNumBuckets(),
					commonRBParameters.getCellSize());

				// Ordenar
				cudppSort(handle, dParticleHashKey, dParticleID, commonRBParameters.st_keyBits,
						mStaticParticlePositions.size());

				// Reclassificar
				gFindCellOffsetOrderingPos(m_dStaticRBGridOffsetStart, m_dStaticRBGridOffsetEnd,
						dParticleHashKey, dParticleID, m_dStaticRBRelativeParticlePositions,
						dParticleOrderedPos, mStaticParticlePositions.size(),
						commonRBParameters.getStNumBuckets());

				// Copiar posicoes ordenadas
				cutilSafeCall(cudaMemcpy(m_dStaticRBRelativeParticlePositions,
						dParticleOrderedPos, sizeof(float4) * mStaticParticlePositions.size(),
						cudaMemcpyDeviceToDevice));

				cudaFree(dParticleOrderedPos);
				cudaFree(dParticleHashKey);
				cudaFree(dParticleID);
				cudppDestroyPlan(handle);
			}
		}

		void RigidBodySimulator::updateDynamicRigidBody()
		{
			if (mDynamicEntities.size() == 0)
				return;

			// Calcular número de corpos rígidos por thread a ser processado na cpu
			int num_procs = omp_get_num_procs();
			int num_elem_thread = (int) (ceil((float)mDynamicEntities.size() / (float) num_procs));

			// Array de posicoes e orientacoes
			mRBPositions = new float4[mDynamicEntities.size()];
			mRBOrientations = new float4[mDynamicEntities.size()];

			// Criar um array de colisão utilizado pelos núcleos da cpu. Cada núcleo é responsável pela atualização
			// de uma parte desse array. Cada núcleo armazena em seu primeiro elemento o total de elementos inseridos
			mRBCollisions = new uint2[num_procs * (num_elem_thread + 1)];

			// Alocar memoria para colisao de corpos rigidos
			Simulation::Application::CommonRigidBodyParams& commonRBParameters =
					Simulation::Application::ApplicationData::instance()->getCommonRigidBodyParams();

			mRBCollisionData.initialize(mDynamicParticlePosRel.size(), mDynamicEntities.size(),
					commonRBParameters.getDnNumBuckets());
			mFluidCollisionData.initialize(mDynamicParticlePosRel.size(), mDynamicEntities.size(),
					commonRBParameters.getDnNumBuckets());


			// Alocar memória
			cutilSafeCall(cudaMalloc((void**) &m_dDynamicRBCenterMassPositions,
					sizeof(float4) * mDynamicEntities.size()));
			//cudaMalloc((void**) &d_DynamicRBOrientation, sizeof(float4) * dynamic_elemens);
			cutilSafeCall(cudaMalloc((void**) &m_dDynamicParticleRelPosition,
					sizeof(float4) * mDynamicParticlePosRel.size()));
			cutilSafeCall(cudaMalloc((void**) &m_dDynamicRBCenterMassVelocity,
					sizeof(float4) * mDynamicEntities.size()));
			cutilSafeCall(cudaMalloc((void**) &m_dDynamicRBMass,
					sizeof(float) * mDynamicEntities.size()));
			cutilSafeCall(cudaMalloc((void**)&m_dDynamicParticlesSizedTmp,
					sizeof(float4) * mDynamicParticlePosRel.size()));
			cutilSafeCall(cudaMalloc((void**)&m_dDynamicParticlesSizedTmp2,
					sizeof(float4) * mDynamicParticlePosRel.size()));
			cutilSafeCall(cudaMalloc((void**)&m_dDynamicRBForces,
					sizeof(float4) * mDynamicEntities.size()));
			cutilSafeCall(cudaMalloc((void**)&m_dDynamicRBTorque,
					sizeof(float4) * mDynamicEntities.size()));
			cutilSafeCall(cudaMalloc((void**)&m_dDynamicRBOrientation,
					sizeof(float4) * mDynamicEntities.size()));
			cutilSafeCall(cudaMalloc((void**)&m_dDynamicAngularMomentum,
					sizeof(float4) * mDynamicEntities.size()));
			cutilSafeCall(cudaMalloc((void**)&m_dDynamicInertiaTensor,
					sizeof(float4) * 3 * mDynamicEntities.size()));
			cutilSafeCall(cudaMalloc((void**)&m_dDynamicInvInertiaTensor,
					sizeof(float4) * 3 * mDynamicEntities.size()));

			// Criar a estrutura de ordenação
			CUDPPConfiguration sortConfig;
			sortConfig.algorithm = CUDPP_SORT_RADIX;
			sortConfig.datatype = CUDPP_UINT;
			sortConfig.op = CUDPP_MIN;
			sortConfig.options = CUDPP_OPTION_KEY_VALUE_PAIRS;
			cudppPlan(&mSortHandle, sortConfig, mDynamicParticlePosRel.size(), 1, 0);


			// Criar estrutura de scan segmentado
			CUDPPConfiguration scanConfig;
			scanConfig.algorithm = CUDPP_SEGMENTED_SCAN;
			scanConfig.datatype = CUDPP_FLOAT4;
			scanConfig.op = CUDPP_ADD;
			scanConfig.options = CUDPP_OPTION_FORWARD | CUDPP_OPTION_INCLUSIVE;
			cudppPlan(&mScanHandle, scanConfig, mDynamicParticlePosRel.size(), 1, 0);


			// Copiar posicoes do centro de massa, orientacao e massa
			cutilSafeCall(cudaMemcpy(m_dDynamicRBCenterMassPositions, &mDynamicCenterMass[0],
					sizeof(float4) * mDynamicCenterMass.size(), cudaMemcpyHostToDevice));

			cutilSafeCall(cudaMemcpy(m_dDynamicRBMass, &mDynamicMass[0],
								sizeof(float) * mDynamicMass.size(), cudaMemcpyHostToDevice));

			cutilSafeCall(cudaMemcpy(m_dDynamicRBOrientation, &mDynamicOrientation[0],
								sizeof(float4) * mDynamicOrientation.size(), cudaMemcpyHostToDevice));

			cutilSafeCall(cudaMemcpy(m_dDynamicInertiaTensor, reinterpret_cast<float4*>(&mInertiaTensor[0]),
								sizeof(float4) * mInertiaTensor.size() * 3, cudaMemcpyHostToDevice));

			cutilSafeCall(cudaMemcpy(m_dDynamicInvInertiaTensor, reinterpret_cast<float4*>(&mInvInertiaTensor[0]),
								sizeof(float4) * mInvInertiaTensor.size() * 3, cudaMemcpyHostToDevice));


			cutilSafeCall(cudaMemset(
					m_dDynamicRBForces, 0, sizeof(float4) * mDynamicEntities.size()));
			cutilSafeCall(cudaMemset(
					m_dDynamicRBCenterMassVelocity, 0, sizeof(float4) * mDynamicEntities.size()));
			cutilSafeCall(cudaMemset(
					m_dDynamicAngularMomentum, 0, sizeof(float4) * mDynamicEntities.size()));
		}

		void RigidBodySimulator::startInteraction()
		{
			assert(mInteractionStarted != true && "Started already called!");

			mInteractionStarted = true;
			mNumCollisionRBRB = 0;
			mNumCollisionRBFluid = 0;
			mNumParticlesRBRBCol = 0;
			mNumParticlesRBFluidCol = 0;
		}

		void RigidBodySimulator::endInteraction()
		{
			assert(mInteractionStarted != false && "Started not called!");
			mInteractionStarted = false;
		}


		/*
		 * Usado somente para corpos rigidos
		 */
		void RigidBodySimulator::generateCollisionData()
		{
			assert(mInteractionStarted != false && "Interaction not started!");

			// Calcular total de corpos rigidos por cpu thread
			int numthreads = omp_get_num_procs();
			int num_elem_thread = (int) (ceil((float)mDynamicEntities.size() /
					(float) numthreads));

			// Resetar flag de colisão
			memset(mRBCollisions, 0,
					sizeof(uint2) * numthreads * (num_elem_thread+1));

			// Inicio do passo paralelo
			#pragma omp parallel default(shared)
			{
				int data_per_threads = (int) floor((float)mDynamicEntities.size() /
						(float)numthreads);
				int current_thread = omp_get_thread_num();
				int start_element = current_thread * data_per_threads;
				int end_element = (current_thread == numthreads - 1) ?
						mDynamicEntities.size() : start_element + data_per_threads;
				int array_offset = current_thread * (num_elem_thread + 1);
				int array_step = 1;

				// Primeiro elemento do array armazena o numero de colisoes existentes
				// no array
				mRBCollisions[array_offset].x = 0;

				for (int i = start_element; i < end_element; i++)
				{
					int jTo = mDynamicEntities.size() - 1;
					bool collided = false;

					// Verificar colisao com outros corpos rigidos dinamicos
					for (int j = 0; j <= jTo; j++)
					{
						if (i == j) continue;

						//if (mDynamicEntities[i]->BBCollide(mDynamicEntities[j]))
						{
							mRBCollisions[array_offset + array_step].x = i;
							mRBCollisions[array_offset + array_step].y = 1;
							collided = true;
							break;
						}
					}


					// Se não colidiu com nenhum corpo rigido dinamico, então verificar
					// se colidiu com algum corpo rigido estatico
					if (!collided && mStaticEntities.size() > 0)
					{
						for (uint w = 0; w < mStaticEntities.size(); w++)
						{
					//		if (mDynamicEntities[i]->BBCollide(mStaticEntities[w]))
							{
								mRBCollisions[array_offset + array_step].x = i;
								mRBCollisions[array_offset + array_step].y = 1;
								collided = true;
								break;
							}
						}
					}

					// Incrementar o contador de colisão
					if (collided)
					{
						mRBCollisions[array_offset].x++;
						array_step++;
					}
				}
			} // Fim passo em paralelo


			/* Verificar os corpos rigidos que colidiram e montar a lista de partículas */
			uint array_rb_rb_col_offset = 0;

			for (int k = 0; k < numthreads; k++)
			{
				int array_offset = k * (num_elem_thread + 1);

				for (uint i = 0; i < mRBCollisions[array_offset].x; i++)
				{
					int rbidx = mRBCollisions[array_offset + i + 1].x;

					// Colidiu com corpo rígido
					if (mRBCollisions[array_offset + i + 1].y == 1)
					{
						// Copiar posicao relativas das particulas do corpo rigido que colidu
						memcpy(&mRBCollisionData.h_CollidedParticleRelativePos[array_rb_rb_col_offset],
							&mDynamicParticlePosRel[mDynamicParticleEndOffset[rbidx] - (mDynamicNumParticles[rbidx] - 1)],
							sizeof(float4) * mDynamicNumParticles[rbidx]);

						// Copiar flags utilizadas no scan segmentado
						memcpy(&mRBCollisionData.h_CollidedParticleFlag[array_rb_rb_col_offset],
							&mDynamicParticleStartFlag[mDynamicParticleEndOffset[rbidx] - (mDynamicNumParticles[rbidx] - 1)],
							sizeof(uint) * mDynamicNumParticles[rbidx]);

						array_rb_rb_col_offset += mDynamicNumParticles[rbidx];
						mNumParticlesRBRBCol += mDynamicNumParticles[rbidx];
						mRBCollisionData.h_CollidedParticleStartOffset[mNumCollisionRBRB] = array_rb_rb_col_offset - 1;
						mNumCollisionRBRB++;

					}
				}
			}


			Simulation::Application::CommonRigidBodyParams& commonRBParameters =
					Simulation::Application::ApplicationData::instance()->getCommonRigidBodyParams();
			//vl::Log::print("Thread: " + vl::String::fromInt(current_thread) +
				//	" From: " + vl::String::fromInt(start_element) + " To: " + vl::String::fromInt(end_element) + "\n");

			// Verificar se é necessário a etapa de narrow phase de corpos rígidos
			if (mNumCollisionRBRB > 0)
			{
				cutilSafeCall(cudaMemcpy(mRBCollisionData.d_CollidedParticlesPos,
					mRBCollisionData.h_CollidedParticleRelativePos,
					sizeof(float4) * mNumParticlesRBRBCol, cudaMemcpyHostToDevice));
				cutilSafeCall(cudaMemcpy(mRBCollisionData.d_CollidedParticleFlag,
					mRBCollisionData.h_CollidedParticleFlag,
					sizeof(uint) * mNumParticlesRBRBCol, cudaMemcpyHostToDevice));
				cutilSafeCall(cudaMemcpy(mRBCollisionData.d_CollidedParticleStartOffset,
					mRBCollisionData.h_CollidedParticleStartOffset,
					sizeof(int) * mNumCollisionRBRB, cudaMemcpyHostToDevice));


				//printGPUData(mRBCollisionData.d_CollidedParticlesPos,mNumParticlesRBRBCol );

				gCalcHashKeyWithOffset(mRBCollisionData.d_CollidedParticlesHashKey,
					mRBCollisionData.d_CollidedParticleId, m_dDynamicRBCenterMassPositions, m_dDynamicRBOrientation,
					mRBCollisionData.d_CollidedParticlesPos, m_dDynamicParticlesSizedTmp,
					mNumParticlesRBRBCol, commonRBParameters.getCellSize(), commonRBParameters.getDnNumBuckets());

				cudppSort(mSortHandle, mRBCollisionData.d_CollidedParticlesHashKey,
						mRBCollisionData.d_CollidedParticleId, commonRBParameters.st_keyBits, mNumParticlesRBRBCol);

				gFindCellOffsetOrderingPos(mRBCollisionData.d_CollidedBucketStart, mRBCollisionData.d_CollidedBucketEnd,
					mRBCollisionData.d_CollidedParticlesHashKey, mRBCollisionData.d_CollidedParticleId, m_dDynamicParticlesSizedTmp,
					mRBCollisionData.d_CollidedParticlesPos, mNumParticlesRBRBCol, commonRBParameters.getDnNumBuckets());
			}

		}



		void RigidBodySimulator::calculateRigidBodyForcesFromRB()
		{
			Simulation::Application::SimulationParams& simParams =
				Simulation::Application::ApplicationData::instance()->getSimulationParams();

			Simulation::Application::CommonRigidBodyParams& commonRBParameters =
				Simulation::Application::ApplicationData::instance()->getCommonRigidBodyParams();

			gCalculateRigidBodyForcesFromRB(mRBCollisionData.d_CollidedParticlesPos,
				m_dDynamicRBCenterMassVelocity, m_dDynamicRBCenterMassPositions,
				mRBCollisionData.d_CollidedParticleTorque, mRBCollisionData.d_CollidedParticlesForces,
				mRBCollisionData.d_CollidedBucketStart, mRBCollisionData.d_CollidedBucketEnd,
				mNumParticlesRBRBCol, commonRBParameters.getCellSize(), commonRBParameters.particleRadius,
				commonRBParameters.getDnNumBuckets(), commonRBParameters.getStNumBuckets(),
				m_dStaticRBRelativeParticlePositions, m_dStaticRBGridOffsetStart, m_dStaticRBGridOffsetEnd,
				-simParams.springForce, -simParams.dampingForce, mStaticEntities.size());

			//printGPUData(mRBCollisionData.d_CollidedParticlesForces, mNumParticlesRBRBCol);

			//printGPUData(mRBCollisionData.d_CollidedParticlesForces, mNumParticlesRBRBCol);
		}

		void RigidBodySimulator::calculateRigidBodyForcesFromFluid(float4* dfluidParticlePos, float4* dfluidParticleVel,
			uint* dfluidHashStart, uint* dfluidHashEnd)
		{
			Simulation::Application::CommonRigidBodyParams& commonRBParameters =
				Simulation::Application::ApplicationData::instance()->getCommonRigidBodyParams();

			Simulation::Application::CommonFluidParams& commonFluidParams =
					Simulation::Application::ApplicationData::instance()->getCommonFluidParams();


			/*gCalculateRigidBodyForcesFromFluids(mFluidCollisionData.d_CollidedParticlesPos,
				m_dDynamicRBCenterMassVelocity, m_dDynamicRBCenterMassPositions, mFluidCollisionData.d_CollidedParticlesForces,
				mFluidCollisionData.d_CollidedParticleTorque, mFluidCollisionData.d_CollidedBucketStart, mFluidCollisionData.d_CollidedBucketEnd, mNumParticlesRBFluidCol,
				commonRBParameters.getCellSize(), commonRBParameters.particleRadius, commonRBParameters.getDnNumBuckets(),
				dfluidParticlePos, dfluidParticleVel, dfluidHashStart, dfluidHashEnd, commonFluidParams.pointRadius * 2.0f,
				commonFluidParams.cell_size, commonFluidParams.getNumBuckets(), commonRBParameters.rb_to_fluid_relation,
				commonRBParameters.spring, commonRBParameters.damping);*/
		}

		void RigidBodySimulator::applyForces()
		{
			Simulation::Application::CommonRigidBodyParams& commonRBParameters =
					Simulation::Application::ApplicationData::instance()->getCommonRigidBodyParams();

			if (mNumCollisionRBRB > 0)
			{
				uint* dRBParticleForceRBIndex;
				uint* dRBParticleForceIndex;

				cutilSafeCall(cudaMalloc((void**)&dRBParticleForceRBIndex,
						sizeof(uint) * mDynamicParticlePosRel.size()));
				cutilSafeCall(cudaMalloc((void**)&dRBParticleForceIndex,
						sizeof(uint) * mDynamicParticlePosRel.size()));

				// classificar forças pelo indice de corpos rígidos a qual pertence
				gCopyRBDataToParticle(dRBParticleForceRBIndex,
					dRBParticleForceIndex, mRBCollisionData.d_CollidedParticlesForces,
					mNumParticlesRBRBCol);

				cudppSort(mSortHandle, dRBParticleForceRBIndex, dRBParticleForceIndex,
					commonRBParameters.st_keyBits, mNumParticlesRBRBCol);

				// Reordenar força e torque
				gReorderForceTorque(m_dDynamicParticlesSizedTmp, m_dDynamicParticlesSizedTmp2,
					mRBCollisionData.d_CollidedParticlesForces, mRBCollisionData.d_CollidedParticleTorque,
					dRBParticleForceIndex, mNumParticlesRBRBCol);

				//// Efetuar o scan
				cudppSegmentedScan(mScanHandle, mRBCollisionData.d_CollidedParticlesForces, m_dDynamicParticlesSizedTmp,
					mRBCollisionData.d_CollidedParticleFlag, mNumParticlesRBRBCol);
				cudppSegmentedScan(mScanHandle, mRBCollisionData.d_CollidedParticleTorque, m_dDynamicParticlesSizedTmp2,
					mRBCollisionData.d_CollidedParticleFlag, mNumParticlesRBRBCol);

				// Adicionar estas força ao total de forças do RB
				gAddForceTorque(mRBCollisionData.d_CollidedParticlesForces,
					mRBCollisionData.d_CollidedParticleTorque, m_dDynamicParticlesSizedTmp,
					m_dDynamicRBForces, m_dDynamicRBTorque,
					mRBCollisionData.d_CollidedParticleStartOffset, mNumCollisionRBRB);

				cutilSafeCall(cudaFree(dRBParticleForceRBIndex));
				cutilSafeCall(cudaFree(dRBParticleForceIndex));
			}

		}

		void RigidBodySimulator::integrate()
		{
			// Verificar se existem corpos rigidos para integrar
			if (mDynamicEntities.size() > 0)
			{
				Simulation::Application::CommonRigidBodyParams& commonRBParameters =
						Simulation::Application::ApplicationData::instance()->getCommonRigidBodyParams();

				gIntegrateRB(m_dDynamicRBCenterMassPositions, m_dDynamicRBMass,
					m_dDynamicRBCenterMassVelocity, m_dDynamicRBForces, mDynamicEntities.size(),
					commonRBParameters.deltaTimeIntegration, commonRBParameters.global_external_force,
					commonRBParameters.volume_min, commonRBParameters.volume_max, m_dDynamicAngularMomentum,
					m_dDynamicRBTorque, m_dDynamicInertiaTensor, m_dDynamicRBOrientation,
					m_dDynamicInvInertiaTensor, commonRBParameters.velocity_limit);


				cutilSafeCall(cudaMemcpy(mRBPositions, m_dDynamicRBCenterMassPositions, sizeof(float4)*
					mDynamicEntities.size(), cudaMemcpyDeviceToHost));

				cutilSafeCall(cudaMemcpy(mRBOrientations, m_dDynamicRBOrientation, sizeof(float4) *
						mDynamicEntities.size(), cudaMemcpyDeviceToHost));

				// Zerar forças
				cutilSafeCall(cudaMemset(m_dDynamicRBForces, 0, sizeof(float4) * mDynamicEntities.size()));
				cutilSafeCall(cudaMemset(m_dDynamicRBTorque, 0, sizeof(float4) * mDynamicEntities.size()));


				// Atualizar posicoes dos corpos rigidos
				 for (int i = 0; i < mDynamicEntities.size(); i++)
				 {
					 RigidBody* rb = mDynamicEntities[i];
					 //rb->setPosition(vl::fvec3(mRBPositions[i].x, mRBPositions[i].y, mRBPositions[i].z));

					 mc::math::Quatf quat(mRBOrientations[i].x, mRBOrientations[i].y,
							 mRBOrientations[i].z, mRBOrientations[i].w);
					 quat.normalizeQuaternion();

					 //rb->getTransform()->setWorldMatrix(
						//vl::fmat4::getTranslation(mRBPositions[i].x, mRBPositions[i].y, mRBPositions[i].z) *
						//quat.toMatrix4() * vl::fmat4::getScaling(3,3,3));
				 }

			}
		}

		void RigidBodySimulator::clearGPUData()
		{
			if (m_dStaticRBRelativeParticlePositions != NULL)
				cutilSafeCall(cudaFree(m_dStaticRBRelativeParticlePositions));

			if (m_dStaticRBGridOffsetStart != NULL)
				cutilSafeCall(cudaFree(m_dStaticRBGridOffsetStart));

			if (m_dStaticRBGridOffsetEnd != NULL)
				cutilSafeCall(cudaFree(m_dStaticRBGridOffsetEnd));

			if (m_dDynamicRBCenterMassPositions != NULL)
				cutilSafeCall(cudaFree(m_dDynamicRBCenterMassPositions));

			if (m_dDynamicParticleRelPosition != NULL)
				cutilSafeCall(cudaFree(m_dDynamicParticleRelPosition));

			if (m_dDynamicRBCenterMassVelocity != NULL)
				cutilSafeCall(cudaFree(m_dDynamicRBCenterMassVelocity));

			if (m_dDynamicRBMass != NULL)
				cutilSafeCall(cudaFree(m_dDynamicRBMass));

			if (m_dDynamicParticlesSizedTmp != NULL)
				cutilSafeCall(cudaFree(m_dDynamicParticlesSizedTmp));

			if (m_dDynamicParticlesSizedTmp2 != NULL)
				cutilSafeCall(cudaFree(m_dDynamicParticlesSizedTmp2));

			if (m_dDynamicRBForces != NULL)
				cutilSafeCall(cudaFree(m_dDynamicRBForces));

			if (m_dDynamicRBOrientation != NULL)
				cutilSafeCall(cudaFree(m_dDynamicRBOrientation));

			if (m_dDynamicAngularMomentum != NULL)
				cutilSafeCall(cudaFree(m_dDynamicAngularMomentum));

			if (m_dDynamicInertiaTensor != NULL)
				cutilSafeCall(cudaFree(m_dDynamicInertiaTensor));

			if (m_dDynamicInvInertiaTensor != NULL)
				cutilSafeCall(cudaFree(m_dDynamicInvInertiaTensor));

			if (mRBCollisions != NULL)
				free(mRBCollisions);

			if (mRBPositions != NULL)
				free(mRBPositions);

			if (mRBOrientations != NULL)
				free(mRBOrientations);


			mRBCollisionData.clear();
			mFluidCollisionData.clear();

		}


		void RigidBodySimulator::reset()
		{
			clearGPUData();
			createHashTable();
			updateStaticElements();
			updateDynamicRigidBody();
		}

		void RigidBodySimulator::release()
		{
			clearGPUData();
		}

	}
}
