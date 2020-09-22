#ifndef __K_MATHUTILS_CU__

#define __K_MATHUTILS_CU__

#include "K_DeviceFunctions.cu"
#include "K_SimulationLib.cuh"
#include <vector_types.h>
#include <vector_functions.h>
#include <cutil_math.h>
#include <cutil_inline.h>



namespace SimulationLib {
	namespace MathUtils {

		/*
		 * Efetuar multiplicacao de matrizes
		 */
		static __device__ void MatrixMultD(float4 *mat1, float4 *mat2, float4 *result)
		{
			result[0].x = mat1[0].x * mat2[0].x + mat1[0].y * mat2[1].x + mat1[0].z * mat2[2].x;
			result[0].y = mat1[0].x * mat2[0].y + mat1[0].y * mat2[1].y + mat1[0].z * mat2[2].y;
			result[0].z = mat1[0].x * mat2[0].z + mat1[0].y * mat2[1].z + mat1[0].z * mat2[2].z;

			result[1].x = mat1[1].x * mat2[0].x + mat1[1].y * mat2[1].x + mat1[1].z * mat2[2].x;
			result[1].y = mat1[1].x * mat2[0].y + mat1[1].y * mat2[1].y + mat1[1].z * mat2[2].y;
			result[1].z = mat1[1].x * mat2[0].z + mat1[1].y * mat2[1].z + mat1[1].z * mat2[2].z;

			result[2].x = mat1[2].x * mat2[0].x + mat1[2].y * mat2[1].x + mat1[2].z * mat2[2].x;
			result[2].y = mat1[2].x * mat2[0].y + mat1[2].y * mat2[1].y + mat1[2].z * mat2[2].y;
			result[2].z = mat1[2].x * mat2[0].z + mat1[2].y * mat2[1].z + mat1[2].z * mat2[2].z;
		}

		/*
		 * Efetuar adicao de matrizes
		 */
		static __device__ void MatrixAddD(float4* mat1, float4* mat2, float4* result)
		{
			result[0].x = mat1[0].x + mat2[0].x;
			result[0].y = mat1[0].y + mat2[0].y;
			result[0].z = mat1[0].z + mat2[0].z;

			result[1].x = mat1[1].x + mat2[1].x;
			result[1].y = mat1[1].y + mat2[1].y;
			result[1].z = mat1[1].z + mat2[1].z;

			result[2].x = mat1[2].x + mat2[2].x;
			result[2].y = mat1[2].y + mat2[2].y;
			result[2].z = mat1[2].z + mat2[2].z;
		}


		/*
		 * Retorna o produto vetorial de um vetor
		 */
		static __device__ float4 Vec3CrossD(float4 vec1, float4 vec2)
		{
			return make_float4((vec1.y * vec2.z) - (vec1.z * vec2.y),
							   (vec1.z * vec2.x) - (vec1.x * vec2.z),
							   (vec1.x * vec2.y) - (vec1.y * vec2.x), 0);
		}

		/*
		 * Retorna o produto escalar entre dois vetores
		 */
		static __device__ float Vec3DotD(float4 vec1, float4 vec2)
		{
			return (vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z);
		}

		static __device__ float4 Quaternion_NormalizeD(float4 q)
		{
			float norm = sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);

			q.x /= norm;
			q.y /= norm;
			q.z /= norm;
			q.w /= norm;

			return q;
		}

		/*
		 * Retorna o conjugado de um quaternion
		 */
		static __device__ float4 Quaternion_conjugateD(float4 quat)
		{
			return make_float4(-quat.x, -quat.y, -quat.z, quat.w);
		}

		/*
		 * Transforma uma matriz em um quaternion
		 */
		static __device__ float4 MatrixToQuaternionD(float4* matrix)
		{
			float4 q;
			float tr, s;

			tr = matrix[0].x + matrix[1].y + matrix[2].z;

			if (tr >= 0)
			{
				s = sqrtf(tr + 1);
				q.w = 0.5 * s;
				s = 0.5 / s;
				q.x = (matrix[2].y - matrix[1].z) * s;
				q.y = (matrix[0].z - matrix[2].x) * s;
				q.z = (matrix[1].x - matrix[0].y) * s;
			}
			else
			{
				int i = 0;

				if (matrix[1].y > matrix[0].x)
					i = 1;

				if (i == 0)
					if (matrix[2].z > matrix[i].x) i = 2;
				if (i == 1)
					if (matrix[2].z > matrix[i].y) i = 2;

				switch (i)
				{
				case 0:
					s = sqrtf((matrix[0].x - (matrix[1].y + matrix[2].z)) + 1);
					q.x = 0.5 * s;
					s = 0.5 / s;
					q.y = (matrix[0].y + matrix[1].x) * s;
					q.z = (matrix[2].x + matrix[0].z) * s;
					q.w = (matrix[2].y - matrix[1].z) * s;
					break;

				case 1:
					s = sqrtf((matrix[1].y - (matrix[2].z + matrix[0].x)) + 1);
					q.y = 0.5 * s;
					s = 0.5 / s;
					q.z = (matrix[1].z + matrix[2].y) * s;
					q.x = (matrix[0].y + matrix[1].x) * s;
					q.w = (matrix[0].z - matrix[2].x) * s;
					break;

				case 2:
					s = sqrtf((matrix[2].z - (matrix[0].x + matrix[1].y)) + 1);
					q.z = 0.5 * s;
					s = 0.5 / s;
					q.x = (matrix[2].x + matrix[0].z) * s;
					q.y = (matrix[1].z + matrix[2].y) * s;
					q.w = (matrix[1].x - matrix[0].y) * s;
				}
			}

			return q;

		}

		/*
		 * Metodo utilizado para calcular o produto entre dois quaternions
		 */
		static __device__ float4 Quaternion_ProductD(float4 quat1, float4 quat2)
		{
			float w = quat1.w * quat2.w - Vec3DotD(quat1, quat2);
			float4 cross = Vec3CrossD(quat1, quat2);

			return make_float4( (quat1.w * quat2.x) + (quat2.w * quat1.x) + cross.x,
				(quat1.w * quat2.y) + (quat2.w * quat1.y) + cross.y,
				(quat1.w * quat2.z) + (quat2.w * quat1.z) + cross.z,
				w);
		}

		/*
		 * Metodo utilizado para transformar uma posicacao de acordo com a
		 * orientacao dada
		 */
		static __device__ float4 Quaternion_TransformPointD(float4 pos, float4 orient)
		{
			pos.w = 0;

			return Quaternion_ProductD(
				Quaternion_ProductD(Quaternion_conjugateD(orient), pos), orient);
		}

		static __device__ void Quaternion_ToMatrixD(float4 q, float4* result)
		{
			result[0].x = 1 - (2 * q.y * q.y) - (2 * q.z * q.z);
			result[0].y = (2 * q.x * q.y) + (2 * q.w * q.z);
			result[0].z = (2 * q.x * q.z) - (2 * q.w * q.y);

			result[1].x = (2 * q.x * q.y) - (2 * q.w * q.z);
			result[1].y = 1 - (2 * q.x * q.x) - (2 * q.z * q.z);
			result[1].z = (2 * q.y * q.z) + (2 * q.w * q.x);

			result[2].x = (2 * q.x * q.z) + (2 * q.w * q.y);
			result[2].y = (2 * q.y * q.z) - (2 * q.w * q.x);
			result[2].z = 1 - (2 * q.x * q.x) - (2 * q.y * q.y);
		}

		/*
		 * Efetuar transposta de uma matriz
		 */
		static __device__ void Matrix_TransposeD(float4 *mat, float4 *result)
		{
			result[0].x = mat[0].x; result[0].y = mat[1].x; result[0].z = mat[2].x;
			result[1].x = mat[0].y; result[1].y = mat[1].y; result[1].z = mat[2].y;
			result[2].x = mat[0].z; result[2].y = mat[1].z; result[2].z = mat[2].z;
		}

		/*
		 * Calcular tensor de inercia
		 */
		static __device__ float4 AngularVelocityD(float4* invTensorInertia, float4 angularMomentum, float4 orientation)
		{
			float4 orientationMatrix[3];
			float4 orientationMatrixTranspose[3];
			float4 temp[3];
			float4 globalTensorInertia[3];


			float4 result;

			// Calcular matriz do quaternion
			Quaternion_ToMatrixD(orientation, orientationMatrix);

			// Calcular matriz transposta
			Matrix_TransposeD(orientationMatrix, orientationMatrixTranspose);

			// Calcular tensor de inércia no espaço global
			MatrixMultD(orientationMatrix, invTensorInertia, temp);
			MatrixMultD(temp, orientationMatrixTranspose, globalTensorInertia);

			result.x = globalTensorInertia[0].x * angularMomentum.x + globalTensorInertia[0].y * angularMomentum.y + globalTensorInertia[0].z * angularMomentum.z;
			result.y = globalTensorInertia[1].x * angularMomentum.x + globalTensorInertia[1].y * angularMomentum.y + globalTensorInertia[1].z * angularMomentum.z;
			result.z = globalTensorInertia[2].x * angularMomentum.x + globalTensorInertia[2].y * angularMomentum.y + globalTensorInertia[2].z * angularMomentum.z;

			return result;
		}
	}
}

#endif
