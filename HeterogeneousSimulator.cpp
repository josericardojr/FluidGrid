/*
 * HeterogeneousSimulator.cu
 *
 *  Created on: 03/02/2011
 *      Author: josericardodasilvajunior
 */

#include "HeterogeneousSimulator.h"
#include "GPU_Kernels/K_FluidSim.cuh"
#include "GPU_Kernels/K_HashGrid.cuh"

Simulation::Simulators::HeterogeneousSimulator*
	Simulation::Simulators::HeterogeneousSimulator::mInstance = NULL;

namespace Simulation
{
	namespace Simulators
	{
		HeterogeneousSimulator* HeterogeneousSimulator::instance()
		{
			if (mInstance == NULL)
				mInstance = new HeterogeneousSimulator();

			return mInstance;
		}

		HeterogeneousSimulator::HeterogeneousSimulator()
		{
			// Setar parametros do kernel
			Simulation::Application::CommonFluidParams& commonParams =
					Application::ApplicationData::instance()->getCommonFluidParams();

			Simulation::Application::KernelsCache hKernelsCache;
			hKernelsCache.kPoly6 = (315.0f / (64.0f * Simulation::Application::PI *
					pow(commonParams.kernel_radius, 9)));
			hKernelsCache.kSpikyGradient = -(45.0f / (Simulation::Application::PI *
					pow(commonParams.kernel_radius, 6)));
			hKernelsCache.kViscosityLaplacian = 45.0f / (Simulation::Application::PI *
					pow(commonParams.kernel_radius, 6));

			gInitializeKernelParams(hKernelsCache);
			Simulation::Application::ApplicationData::instance()->setKernelsCache(hKernelsCache);

			// Setar parametros comuns da simulacao de fluidos
			gInitializeCommonFluidParams(commonParams);
			
		}

		void HeterogeneousSimulator::addSphFluid(Fluid::GPUSPHFluid* fluid)
		{
			mFluids.push_back(fluid);
		}

		void HeterogeneousSimulator::doTimeStep()
		{
			// Iniciar interacao
			//RigidBody::RigidBodySimulator::instance()->startInteraction();

			
			for (int i = 0; i < mFluids.size(); i++)
			{
				JRFXGL::Util::Clock clock;
				// Atualizar fluido
				mFluids[i] -> beginUpdate();

				clock.reset();
				findFluidNeighbourhood(mFluids[i]);
				

				// Sincronizar
				cutilSafeCall(cudaThreadSynchronize());
				printf("particle structuring: %f\n", clock.getMilliseconds());

				// Efetuar em paralelo a verificacao de colisao e o calculo da densidade
				// todo Na arquitetura fermi, colocamos essa etapa em um nucleo e a etapa seguinte em outro nucleo
				clock.reset();
				calculateFluidDensity(mFluids[i]);
				cutilSafeCall(cudaThreadSynchronize());
				printf("Fluid Density: %f\n", clock.getMilliseconds());

				// Verificar se houve interacao desse fluido com corpos rigidos (em paralelo)
			//	RigidBody::RigidBodySimulator::instance()->generateCollisionData(
				//		mFluids[i]->getAABB());

				// Verificar se e necessario interacao
	/*			if (RigidBody::RigidBodySimulator::instance()->getNumCollisionsWithFluid() > 0)
				{
					// Calculamos as forcas a serem aplicadas nos fluidos e corpos rigidos
					// todo como cada atualizacao de forcas nao existe dependencia, podemos colocar cada atualizacao em um nucleo da fermi

					RigidBody::RigidBodySimulator::instance()->calculateRigidBodyForcesFromFluid(
						mFluids[i]->getPosition(), mFluids[i]->getVelocity(),
						mFluids[i]->getHashStart(), mFluids[i]->getHashEnd());

					calculateFluidForces(mFluids[i],
						RigidBody::RigidBodySimulator::instance()->getNumStaticRigidBodies() > 0,
						true);

					// Sincronizar
					cutilSafeCall(cudaThreadSynchronize());
				}
				else*/
				{
					clock.reset();
					calculateFluidForces(mFluids[i], false,
//						RigidBody::RigidBodySimulator::instance()->getNumStaticRigidBodies() > 0,
						false);
					cudaThreadSynchronize();
					printf("forces: %f\n", clock.getMilliseconds());
				}

				// Aplicar forcas
				//RigidBody::RigidBodySimulator::instance()->applyForces();

				// Integrar
				//RigidBody::RigidBodySimulator::instance()->integrate();

				clock.reset();
				integrateFluidForces(mFluids[i]);
				cudaThreadSynchronize();
				printf("integration: %f\n", clock.getMilliseconds());

				mFluids[i]-> endUpdate();
			}

			//RigidBody::RigidBodySimulator::instance()->endInteraction();
		}



		void HeterogeneousSimulator::findFluidNeighbourhood(Fluid::GPUSPHFluid *fluid)
		{
			Simulation::Application::CommonFluidParams& commonParams =
					Application::ApplicationData::instance()->getCommonFluidParams();

			gCalcHashKey(fluid->getHashKey(), fluid->getParticleID(), fluid->getPosition(),
				fluid->getNumParticles(), commonParams.getNumBuckets(), commonParams.cell_size);



			// Efetuar o sort
			gSort(fluid->getHashKey(), fluid->getParticleID(), fluid->getNumParticles());
			

			// Calcular offset de cada celula
			gFindCellOffsetOrderingPosVel(fluid->getHashStart(), fluid->getHashEnd(), 
				fluid->getHashKey(), fluid->getParticleID(), fluid->getPosition(), 
				fluid->getOrderedPosition(), fluid->getVelocity(), fluid->getOrderedVelocities(),
				fluid->getNumParticles(), commonParams.getNumBuckets());
		}


		void HeterogeneousSimulator::calculateFluidDensity(Fluid::GPUSPHFluid *fluid)
		{
			gCalcMassDensity(fluid->getOrderedPosition(), fluid->getDensityField(),
				fluid -> getHashKey(), fluid -> getHashStart(), 
				fluid -> getHashEnd(), fluid -> getNumParticles(), fluid->getGasConstant(), 
				fluid->getRestDensity(), fluid->getMass());
		}

		void HeterogeneousSimulator::calculateFluidForces(Fluid::GPUSPHFluid *fluid,
				bool useStaticRB, bool useDynamicRB)
		{
			Simulation::Application::CommonFluidParams& commonFluidParams =
				Application::ApplicationData::instance()->getCommonFluidParams();

			Simulation::Application::CommonRigidBodyParams& commonRigidBodyParams =
				Application::ApplicationData::instance()->getCommonRigidBodyParams();

			//RigidBody::RigidBodySimulator* rbInstance = Simulation::RigidBody::RigidBodySimulator::instance();

			// Verificar se existem corpos rigidos estaticos na cena
			//if (useStaticRB && useDynamicRB)
		/*	{
				RigidBody::CollisionData* fluidCollisionData = rbInstance->getFluidCollisionData();
				gSPHCalculateForcesWithSRB(fluid->getOrderedPosition(), fluid->getOrderedVelocities(),
					fluid->getForces(), fluid->getVelocity(), fluid->getDensityField(), fluid->getParticleID(),
					fluid->getHashStart(), fluid->getHashEnd(), fluid->getNumParticles(), fluid->getMass(),
					fluid->getViscosity(), commonFluidParams.pointRadius * 2.0f, commonFluidParams.fluid_to_srb_relation, 0.5f * 3.0f,
					commonRigidBodyParams.getCellSize(), -10, -5, Simulation::RigidBody::RigidBodySimulator::instance()->getStaticRBParticlePositions(),
					Simulation::RigidBody::RigidBodySimulator::instance()->getStaticRBGridStartOffset(),
					Simulation::RigidBody::RigidBodySimulator::instance()->getStaticRBGridEndOffset(),
					commonRigidBodyParams.getStNumBuckets(), fluidCollisionData->d_CollidedParticlesPos,
					rbInstance->getDynamicVelocities(), fluidCollisionData->d_CollidedBucketStart,
					fluidCollisionData->d_CollidedBucketEnd, commonRigidBodyParams.getDnNumBuckets(), useStaticRB, useDynamicRB);
			}*/


			gSPHCalculateForces(fluid->getOrderedPosition(), fluid->getOrderedVelocities(), 
				fluid->getForces(), fluid->getDensityField(), fluid->getHashStart(), 
				fluid->getHashEnd(), fluid->getVelocity(), fluid->getParticleID(),
				fluid->getNumParticles(), fluid->getMass(), fluid->getViscosity());
				
				


		/*	gSPHCalculateForcesWithSRB(fluid->getOrderedPosition(), fluid->getOrderedVelocities(),
					fluid->getForces(), fluid->getVelocity(), fluid->getDensityField(), fluid->getParticleID(),
					fluid->getHashStart(), fluid->getHashEnd(), fluid->getNumParticles(), fluid->getMass(),
					fluid->getViscosity(), commonFluidParams.pointRadius, commonFluidParams.fluid_to_srb_relation, 5.5f,
					commonRigidBodyParams.st_cell_size, 20.5f, 20.5f, Simulation::RigidBody::RigidBodySimulator::instance()->getStaticRBParticlePositions(),
					Simulation::RigidBody::RigidBodySimulator::instance()->getStaticRBGridStartOffset(),
					Simulation::RigidBody::RigidBodySimulator::instance()->getStaticRBGridEndOffset(),
					commonRigidBodyParams.getStNumBuckets());*/
		}

		void HeterogeneousSimulator::integrateFluidForces(Fluid::GPUSPHFluid *fluid)
		{
			gSPHIntegrate(fluid->getPosition(), fluid->getVelocity(), fluid->getForces(),
				fluid->getNumParticles(), fluid->getDensityField(), fluid->getHashStart(), 
				fluid->getHashEnd());
		}
	}
}
