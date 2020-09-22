//---------------------------------------------------------------------------
#ifndef __APP_MAIN_H__
#define __APP_MAIN_H__
//---------------------------------------------------------------------------

#include "jrfxgl.h"
#include "GPUSPHFluid.h"
#include "Emitter.h"
#include "HeterogeneousSimulator.h"
#include "FluidParticleRenderMaterial.h"

	class AppMain : public JRFXGL::App::AppRunner
	{
	  public:
		virtual ~AppMain ();
		
	  public:
		AppMain ();

	  protected:
		virtual void OnRender       ();
		virtual bool OnInit         ();
        virtual void OnGUI          ();
		virtual void OnUpdate       (float dt);
		virtual void OnKeyDown(const SDL_keysym & keysym)
		{
			if (keysym.sym == SDLK_ESCAPE)
				OnQuit();

			if (keysym.sym == SDLK_p)
				paused = !paused;
		}

	  private:
		void initializeFluid();

	  private:
		JRFXGL::Graphics::FPSCameraController camera;
		JRFXGL::Graphics::Node node;
		JRFXGL::Graphics::Mesh *plane;
		boost::shared_ptr<JRFXGL::Graphics::Material> m;
		JRFXGL::Graphics::DirectionalLight dl;
		Simulation::Fluid::GPUSPHFluid*      mFluid;
		Simulation::Emitter::IEmitter*  mBlockEmitter;
		bool paused;
	};


//---------------------------------------------------------------------------
#endif

