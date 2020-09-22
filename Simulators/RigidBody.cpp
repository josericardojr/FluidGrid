/*
 * RigidBody.cpp
 *
 *  Created on: 03/02/2011
 *      Author: josericardodasilvajunior
 */

#include "RigidBody.h"


namespace JRFXGL
{
	namespace RigidBodySim
	{
		/*RigidBodyParticleRenderer::RigidBodyParticleRenderer(vl::ArrayFloat3* array, float radius)
		{
			mArray = array;
			mRadius = radius;

			sphere = gluNewQuadric();
			gluQuadricNormals(sphere, GLU_SMOOTH);
			gluQuadricTexture(sphere, GL_TRUE);
		}

		void RigidBodyParticleRenderer::render_Implementation(const vl::Actor* actor,
			const vl::Shader* shader, const vl::Camera* camera, vl::OpenGLContext* gl_context) const
		{
			Geometry::render_Implementation(actor, shader, camera, gl_context);
			/*VL_CHECK_OGL()

			//glPushMatrix();
			//const vl::mat4& m = actor->transform()->localMatrix();
			//glMultMatrixf(m.ptr());
			for (int i = 0; i < mArray->size(); i++)
			{
				glPushMatrix();
				glTranslatef(mArray->at(i).x(), mArray->at(i).y(), mArray->at(i).z());
				gluSphere(sphere, mRadius, 16, 8);
				glPopMatrix();
			}
			//glPopMatrix();
			VL_CHECK_OGL()
		}*/


		RigidBody::RigidBody()
		{
		}

		RigidBody::~RigidBody()
		{
			delete mRelativeParticlePositions;
		}

		RigidBody* RigidBody::loadFromFile(std::string filename)
		{
			RigidBody *rb = new RigidBody();
			const char* cc = filename.c_str();

			// Abrir arquivo
			std::ifstream file(cc);

			std::string buffer;

			// Recuperar tipo: estatico ou dinamico
			getline(file, buffer);
			//buffer = vl::String::trimStdString(buffer);
			if (buffer.compare("static") == 0)
			{
				rb->mIsStatic = true;
				rb->mMass = 0;
			}
			else
			{
				getline(file, buffer);
				//buffer = vl::String::trimStdString(buffer);
				rb->mIsStatic = false;
				rb->mMass = atof(buffer.c_str());
			}

			// Recuperar raio
			getline(file, buffer);
			//buffer = vl::String::trimStdString(buffer);
			rb->mParticleRadius = atof(buffer.c_str());

			// recuperar mesh original
			getline(file, buffer);
			//buffer = vl::String::trimStdString(buffer);
			rb->mOriginalMeshFilename = buffer;

			// recuperar o numero de pontos
			getline(file, buffer);
			//buffer = vl::String::trimStdString(buffer);
			rb->mNumParticles = atof(buffer.c_str());
			rb->mRelativeParticlePositions = new float4[rb->mNumParticles];

			// recuperar posicoes
			for (int i = 0; i < rb->mNumParticles; i++)
			{
				getline(file, buffer);
				std::istringstream line(buffer);
				std::string s1, s2, s3;
				line >> s1>>s2>>s3;
				rb->mRelativeParticlePositions[i].x = atof(s1.c_str());
				rb->mRelativeParticlePositions[i].y = atof(s2.c_str());
				rb->mRelativeParticlePositions[i].z = atof(s3.c_str());
				rb->mRelativeParticlePositions[i].w = rb->mParticleRadius;
			}

			file.close();

			// Criar scene node e carregar modelo
			//rb->attachEntity(JRFXGL::Graphics::ModelLoader::loadObjModel(rb->mOriginalMeshFilename.c_str()));
			//rb->getEntity()->setMaterial(new JRFXGL::Graphics::Material());

			rb->inertiaTensor.Ix = make_float4(200 * (rb->mMass / 12), 0.0f, 0.0f, 0.0f);
			rb->inertiaTensor.Iy = make_float4(0.0f, 200 * (rb->mMass / 12), 0.0f, 0.0f);
			rb->inertiaTensor.Iz = make_float4(0.0f, 0.0f, 200 * (rb->mMass / 12), 0.0f);

			rb->invInertiaTensor.Ix = make_float4(0.00125f, 0.0f, 0.0f, 0.0f);
			rb->invInertiaTensor.Iy = make_float4(0.0f, 0.00125f, 0.0f, 0.0f);
			rb->invInertiaTensor.Iz = make_float4(0.0f, 0.0f, 0.00125f, 0.0f);

			//rb->calcBounds();
			return rb;
		}

	/*	void RigidBody::calcBounds()
		{
			vl::fvec3 minBB;
			vl::fvec3 maxBB;


			for (int i = 0; i < mNumParticles; i++)
			{

				if (i == 0)
				{
					minBB = maxBB = vl::fvec3(mRelativeParticlePositions[i].x,
							mRelativeParticlePositions[i].y,mRelativeParticlePositions[i].z);
				}
				else
				{
					if (mRelativeParticlePositions[i].x > maxBB.x()) maxBB.x() = mRelativeParticlePositions[i].x;
					if (mRelativeParticlePositions[i].y > maxBB.y()) maxBB.y() = mRelativeParticlePositions[i].y;
					if (mRelativeParticlePositions[i].z > maxBB.z()) maxBB.z() = mRelativeParticlePositions[i].z;

					if (mRelativeParticlePositions[i].x < minBB.x()) minBB.x() = mRelativeParticlePositions[i].x;
					if (mRelativeParticlePositions[i].y < minBB.y()) minBB.y() = mRelativeParticlePositions[i].y;
					if (mRelativeParticlePositions[i].z < minBB.z()) minBB.z() = mRelativeParticlePositions[i].z;
				}
			}

			mBoundingBox.setMinCorner(minBB);
			mBoundingBox.setMaxCorner(maxBB);
		}

		bool RigidBody::BBCollide(vl::AABB aabb)
		{
			//return boundingBoxSafe().intersects(aabb);
			//return boundingBoxSafe().intersects(aabb);
			return true;
		}

		bool RigidBody::BBCollide(RigidBody* rb)
		{
			//return mBoundingBox.intersects(rb->getBB());
			//return BBCollide(rb->getBB());
			return true;
		}*/
	}
}
