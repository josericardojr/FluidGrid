

//#include "SDL.h"
#include "AppMain.h"
//#include "../../Inc/Util/Error.h"
//#include "../../Inc/Lua/LuaState.h"
//#include "../../Inc/GUI/GUIManager.h"



	//---------------
	// CLASSES
	//---------------

	/**
	 *
	 */
	AppMain::AppMain ()
	{
		node.attachEntity(JRFXGL::Graphics::ModelLoader::loadObjModel("resources/teapot.obj"));
		m.reset(new JRFXGL::Graphics::Material());
		node.getEntity()->setMaterial(m);
		paused = true;

		//plane = JRFXGL::Graphics::ModelLoader::loadObjModel("resources/plane.obj");
		//plane->setMaterial(&m);
	}
		

	/**
	 *
	 */
	AppMain::~AppMain ()
	{ 
		//delete plane;
		delete mFluid;
		delete mBlockEmitter;
	}

	bool AppMain::OnInit()
	{
		cutilSafeCall(cudaDeviceReset());
		cutilSafeCall(cudaGLSetGLDevice(0));

		camera.getCamera().getProjection().setAspectRatio(800, 600);
		camera.mouseSpeed = 60;
		camera.speed = 500;
		initializeFluid();

		return true;
	}

	void AppMain::OnUpdate(float dt)
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_DEPTH_TEST);
		dl.enable();
		camera.update(dt);

		if (!paused)
			Simulation::Simulators::HeterogeneousSimulator::instance()->doTimeStep();
	}

    void AppMain::OnGUI()
    {
       JRFXGL::Graphics::GUIManager& mgr = JRFXGL::Graphics::GUIManager::getInstance();
        mgr.label(std::string("teste"), 10, 20);
    }

	void AppMain::OnRender       ()
	{
		glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
		camera.getCamera().applyTransform();
		dl.apply();
		mFluid->render();
	}

	void AppMain::initializeFluid()
	{
		// Reinicializar GPU

		// Parametros da simulacao
		Simulation::Application::CommonFluidParams commonFluidParams;
		commonFluidParams.volume_min = make_float3(-20.0f, -20.0f, -20.0f);
		commonFluidParams.volume_max = make_float3(20.0f, 85.0f, 40.0f);
		commonFluidParams.global_external_force = make_float3(0, -9.8f, 0);
		commonFluidParams.simulation_scale = 0.009f;
		commonFluidParams.deltaTimeIntegration = 0.004f;
		commonFluidParams.bucket_size = make_int3(64, 64, 64);
		commonFluidParams.keyBits = 24;
		commonFluidParams.velocity_limit = 600;
		commonFluidParams.pointScale = (float)(800/ tanf(90.0*0.5f*3.1415f/180.0f ));
		commonFluidParams.fluid_to_srb_relation = make_int3(1, 1, 1);

		// Parametros do fluido
		Simulation::Application::FluidParams fluidParams;
		fluidParams.max_particles = 24000;
		fluidParams.mass = ((128*1024.0f)/fluidParams.max_particles) * 0.0002f;
		fluidParams.restDensity = 1000.0f;
		fluidParams.stiffness = 0.2f;
		fluidParams.viscosity = 1.0f;
		fluidParams.gasConstant = 1.5;
		fluidParams.particle_radius = 0.3f;

		// Distancia entre as particulas
		float dist = pow(fluidParams.mass / fluidParams.restDensity, 1.0f/3.0f);
		float spacing = dist *0.87f / commonFluidParams.simulation_scale;
		commonFluidParams.kernel_radius = 2.0f * dist * 0.87f;       // m
		commonFluidParams.kernel_radius_squared = commonFluidParams.kernel_radius *
				commonFluidParams.kernel_radius;
		commonFluidParams.pointRadius = spacing;

		float pv = commonFluidParams.kernel_radius / commonFluidParams.simulation_scale;
		mc::math::Vec3f p( pv * 1.0f, pv * 1.0f, pv * 1.0f);
		commonFluidParams.cell_size = make_float3( p.x, p.y, p.z);

		Simulation::Application::ApplicationData::instance()->setCommonFluidParams(
				commonFluidParams);

		// Emissor de partículas
		mBlockEmitter = new Simulation::Emitter::BlockEmitter(
			mc::math::Vec3f(-20.0f, -20.0f, -20.0f),
			mc::math::Vec3f(20.0f, 85.0f, 0.0f), spacing);

		mFluid = new Simulation::Fluid::GPUSPHFluid(fluidParams, mBlockEmitter);

		Simulation::Fluid::FluidParticleMaterialPtr fluidMaterial(
				new Simulation::Fluid::FluidParticleMaterial());
		fluidMaterial->setPointScale((float)(800.0f / tanf(90.0f*0.5f*3.1415f/180.0f)));
		fluidMaterial->setPointRadius(commonFluidParams.pointRadius);
		mFluid->setMaterial(fluidMaterial);
		Simulation::Simulators::HeterogeneousSimulator::instance()->addSphFluid(mFluid);
	}



