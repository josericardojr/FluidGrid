/*
 * ApplicationData.cpp
 *
 *  Created on: 31/01/2011
 *      Author: josericardodasilvajunior
 */
#include "ApplicationData.cuh"

Simulation::Application::ApplicationData*
	Simulation::Application::ApplicationData::mInstance = NULL;

namespace Simulation
{
	namespace Application
	{

		ApplicationData* ApplicationData::instance()
		{
			if (mInstance == NULL)
				mInstance = new ApplicationData();

			return mInstance;
		}
	}
}
