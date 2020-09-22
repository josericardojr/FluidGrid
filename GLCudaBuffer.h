#ifndef __VLSimBuffer_h__
#define __VLSimBuffer_h__

#include <JRFXGL.h>
#include <cuda_gl_interop.h>
#include <cutil_inline.h>

namespace Simulation
{
	template <class T>
	class GLCudaBuffer
	{
	public:
		GLCudaBuffer() :
			mBuffer(NULL),
			mMapped(false)
		{
		}

		~GLCudaBuffer()
		{
			if (mBuffer)
			{
				unmapBuffer();

				// Remover registro
				cutilSafeCall(cudaGraphicsUnregisterResource(this->mBufferCudaMap));
				mBuffer = NULL;
			}
		}

		void setBuffer(JRFXGL::Graphics::GLVertexBuffer*  vlBuffer)
		{
			if (mBuffer)
			{
				unmapBuffer();

				// Remover registro
				cutilSafeCall(cudaGraphicsUnregisterResource(mBufferCudaMap));
				mBuffer = NULL;
			}

			mBuffer = vlBuffer;

			// Registrar buffer
			cutilSafeCall(cudaGraphicsGLRegisterBuffer(&this->mBufferCudaMap,
					mBuffer->id(), cudaGraphicsMapFlagsNone));
		}

		void mapBuffer()
		{
			assert(mBuffer && "JRFXGL::Graphics::GLVertexBuffer is NULL");

			if (!mMapped)
			{

				cutilSafeCall(cudaGraphicsMapResources(1, &this->mBufferCudaMap, 0));

				size_t num_bytes;
				cutilSafeCall(cudaGraphicsResourceGetMappedPointer((void**)&mDevice_data,
						&num_bytes, this->mBufferCudaMap));

				mMapped = true;
			}

		}

		void unmapBuffer()
		{
			if (mMapped)
			{
				assert(mBuffer && "JRFXGL::Graphics::GLVertexBuffer is NULL");

				cutilSafeCall(cudaGraphicsUnmapResources(1, &this->mBufferCudaMap, 0));

				mMapped = false;
			}
		}

		T getData(){ return mDevice_data; }

	private:
		JRFXGL::Graphics::GLVertexBuffer* mBuffer;
		struct cudaGraphicsResource*      mBufferCudaMap;
		T                                 mDevice_data;
		bool                              mMapped;
	};
	
	typedef GLCudaBuffer<mc::math::Vec4f*> GLCudaBufferVec4;
}
#endif
