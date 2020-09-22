/*
 * RigidBody.h
 *
 *  Created on: 03/02/2011
 *      Author: josericardodasilvajunior
 */

#ifndef RIGIDBODY_H_
#define RIGIDBODY_H_

#include <vector_types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector_functions.h>
#include "jrfxgl.h"

namespace JRFXGL
{
	namespace RigidBodySim
	{
		struct sInertiaTensor
		{
			float4 Ix;
			float4 Iy;
			float4 Iz;
		};

	/*	class RigidBodyParticleRenderer : public vl::Geometry
		{
		public:
			RigidBodyParticleRenderer(vl::ArrayFloat3* array, float radius);

		protected:
			virtual void render_Implementation(const vl::Actor* actor,
				const vl::Shader* shader, const vl::Camera* camera, vl::OpenGLContext* gl_context) const;

		private:
			vl::ref<vl::ArrayFloat3> mArray;
			float           mRadius;
			GLUquadricObj* sphere;
		};*/

		/*
		 * Elemento rigid body que pode ser simulado
		 */
		class RigidBody : public JRFXGL::Graphics::Node
		{
		private:
			RigidBody();

		public:
			~RigidBody();
			static RigidBody* loadFromFile(std::string filename);


			float4*       getRelativeParticlePositions(){ return mRelativeParticlePositions; }
			int           getNumParticles(){ return mNumParticles; }
			float         getParticleRadius(){ return mParticleRadius; }
			float         getMass(){ return mMass; }
			bool          isStatic(){ return mIsStatic; }
			//bool          BBCollide(vl::AABB aabb);
			//bool          BBCollide(RigidBody* rb);
			sInertiaTensor getInertiaTensor(){ return inertiaTensor; }
			sInertiaTensor getInvInertiaTensor(){ return invInertiaTensor; }
			//void           calcBounds();
			//vl::AABB      getBB(){ return mBoundingBox.transformed(vl::mat4::getTranslation(getPosition())); }


		private:
			float4*          mRelativeParticlePositions;
			int              mNumParticles;
			bool             mIsStatic;
			float            mParticleRadius;
			float            mMass;
			std::string      mOriginalMeshFilename;
			sInertiaTensor   inertiaTensor;
			sInertiaTensor   invInertiaTensor;
			//vl::AABB         mBoundingBox;
		};
	}
}

#endif /* RIGIDBODY_H_ */
