/*
 * HeterogeneousSimulator.cuh
 *
 *  Created on: 02/02/2011
 *      Author: josericardodasilvajunior
 */

#ifndef __HeterogeneousSimulator_CUH__
#define __HeterogeneousSimulator_CUH__

#include <vector>
#include "GPUSPHFluid.h"
#include "ApplicationData.cuh"
//#include "RigidBodySimulator.h"

namespace Simulation
{
	namespace Simulators
	{
		/*
		 * Classe responsavel pelo gerenciamento das tarefas da simulacao de corpos rigidos e
		 * fluidos
		 */
		class HeterogeneousSimulator
		{
		private:
			HeterogeneousSimulator();

		private:
			/* Metodos relacioinados com a simulacao de fluidos */
			void findFluidNeighbourhood(Fluid::GPUSPHFluid *fluid);
			void calculateFluidDensity(Fluid::GPUSPHFluid *fluid);
			void calculateFluidForces(Fluid::GPUSPHFluid *fluid,
				bool useStaticRB, bool useDynamicRB);
			void integrateFluidForces(Fluid::GPUSPHFluid *fluid);




		public:
			static HeterogeneousSimulator* instance();

			void addSphFluid(Fluid::GPUSPHFluid* fluid);

			void doTimeStep();


		private:
			std::vector<Simulation::Fluid::GPUSPHFluid*>         mFluids;
			static HeterogeneousSimulator*    mInstance;


		};
	}
}

#endif /* __HeterogeneousSimulator_CUH__ */
