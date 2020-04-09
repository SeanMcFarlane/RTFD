#include "project.h"

namespace CUDA { // CUDA implementation
	
	__global__
	void gpu_lin_solve(float* x, float* x0, float a, float c, int N, int pad, int bnd) {
		int const dim = N + bnd;
		int const i = blockIdx.x*blockDim.x+threadIdx.x;
		int const j = blockIdx.y*blockDim.y+threadIdx.y;
		if (i < dim-pad && j < dim-pad && i > pad && j > pad)
		{
			const float left = x[IX(i - 1, j)];
			const float right = x[IX(i + 1, j)];
			const float above = x[IX(i, j + 1)];
			const float below = x[IX(i, j - 1)];
			const float sum = left + right + above + below;
			const float prod = sum * a;
			const float addOrig = prod + x0[IX(i, j)];
			const float result = addOrig / c;
			x[IX(i, j)] = result;
		}
	}

	__global__
	void gpu_set_bnd(int N, int pad, int bnd, uint32_t b, float *x)
	{
		//int const dim = N + bnd;
		int const i = blockIdx.x*blockDim.x+threadIdx.x;
		int const j = blockIdx.y*blockDim.y+threadIdx.y;

		if (i == pad-1) {
			x[IX(i, j)] = (b == 1) ? -x[IX(i+1, j)] : x[IX(i+1, j)];
		}		
		
		else if (i == N+pad) {
			x[IX(i, j)] = (b == 1) ? -x[IX(i-1, j)] : x[IX(i-1, j)];
		}

		else if (j == pad-1) {
			x[IX(i, j)] = (b == 2) ? -x[IX(i, j+1)] : x[IX(i, j+1)];
		}		

		else if (j == N+pad) {
			x[IX(i, j)] = (b == 2) ? -x[IX(i, j-1)] : x[IX(i, j-1)];
		}

		//x[IX(pad-1, pad-1)] = 0.5f * (x[IX(pad, pad-1)] + x[IX(pad-1, pad)]);
		//x[IX(pad-1, N + pad)] = 0.5f * (x[IX(pad, N+pad)] + x[IX(pad-1, N+pad-1)]);
		//x[IX(N+pad, pad-1)] = 0.5f * (x[IX(N+pad-1,pad-1)] + x[IX(N+pad, pad)]);
		//x[IX(N+pad, N+pad)] = 0.5f * (x[IX(N+pad-1, N+pad)] + x[IX(N+pad, N+pad-1)]);
	}

	__host__
	void add_source(uint32_t N, float *x, float *s, float dt)
	{
		uint32_t i, size = (N + bnd) * (N + bnd);
		for (i = 0; i < size; i++)
			x[i] += dt * s[i];
	}


	void lin_solve(uint32_t N, uint32_t b, float *x, float *x0, float a, float c)
	{
		//const int blocksize = 64;
		auto const size = sizeof(float)*array_size;
		//printf("Size = %i\n", (int)size);
		//printf("array_size = %i\n", (int)array_size);

		//auto const num_blocks = ceil(array_size / static_cast<int>(blocksize));

		float *dev_x, *dev_x0;
		cudaMalloc((void**)&dev_x, size);
		cudaMalloc((void**)&dev_x0, size);
		
		dim3 blockSize = dim3(8, 8);

		int dim = N + bnd;
		int bx = (dim + blockSize.x - 1) / blockSize.x;
		int by = (dim + blockSize.y - 1) / blockSize.y;
		dim3 gridSize = dim3(bx, by);

		uint32_t k;

		cudaMemcpy(dev_x, x, size, cudaMemcpyHostToDevice);
		cudaMemcpy(dev_x0, x0, size, cudaMemcpyHostToDevice);

		for (k = 0; k < 20; k++)
		{
			CUDA::gpu_lin_solve<<<gridSize, blockSize>>>(dev_x, dev_x0, a, c, N, pad, bnd);
			CUDA::gpu_set_bnd<<<gridSize, blockSize>>>(N, pad, bnd, b, dev_x);
		}

		cudaMemcpy(x, dev_x, size, cudaMemcpyDeviceToHost);
		cudaMemcpy(x0, dev_x0, size, cudaMemcpyDeviceToHost);

		cudaFree(dev_x);
		cudaFree(dev_x0);
	}


	void set_bnd(uint32_t N, uint32_t b, float *x)
	{
		uint32_t i;
		uint32_t pad = bnd/2;
		for (i = pad; i < N+pad; i++)
		{
			x[IX(pad-1, i)] = b == 1 ? -x[IX(pad, i)] : x[IX(pad, i)];
			x[IX(N+pad, i)] = b == 1 ? -x[IX(N+pad-1, i)] : x[IX(N+pad-1, i)];
			x[IX(i, pad-1)] = b == 2 ? -x[IX(i, pad)] : x[IX(i, pad)];
			x[IX(i, N+pad)] = b == 2 ? -x[IX(i, N+pad-1)] : x[IX(i, N+pad-1)];
		}
		x[IX(pad-1, pad-1)] = 0.5f * (x[IX(pad, pad-1)] + x[IX(pad-1, pad)]);
		x[IX(pad-1, N + pad)] = 0.5f * (x[IX(pad, N+pad)] + x[IX(pad-1, N+pad-1)]);
		x[IX(N+pad, pad-1)] = 0.5f * (x[IX(N+pad-1,pad-1)] + x[IX(N+pad, pad)]);
		x[IX(N+pad, N+pad)] = 0.5f * (x[IX(N+pad-1, N+pad)] + x[IX(N+pad, N+pad-1)]);
	}


	//void lin_solve(uint32_t N, uint32_t b, float* x, float* x0, float a, float c)
	//{
	//	uint32_t i, j, k;
	//	for (k = 0; k < 20; k++){
	//		for ( j=4; j<=N+4; j++ ) { 
	//			for ( i=4; i<=N+4; i+=4 ) {
	//				x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
	//			}
	//		}
	//		CUDA::set_bnd(N, b, x);
	//	}
	//}

	void diffuse(uint32_t N, uint32_t b, float *x, float *x0, float diff, float dt)
	{
		float a = dt * diff * N * N;
		CUDA::lin_solve(N, b, x, x0, a, 1 + 4 * a);
	}

	void advect(uint32_t N, uint32_t b, float *d, float *d0, float *u, float *v, float dt)
	{
		uint32_t i, j, i0, j0, i1, j1;
		float x, y, s0, t0, s1, t1, dt0;

		dt0 = dt * N;
		FOR_EACH_CELL
			x = i - dt0 * u[IX(i, j)];
			y = j - dt0 * v[IX(i, j)];
			if (x < 0.5f)
				x = 0.5f;
			if (x > N + 0.5f)
				x = N + 0.5f;
			i0 = (int)x;
			i1 = i0 + 1;
			if (y < 0.5f)
				y = 0.5f;
			if (y > N + 0.5f)
				y = N + 0.5f;
			j0 = (int)y;
			j1 = j0 + 1;
			s1 = x - i0;
			s0 = 1 - s1;
			t1 = y - j0;
			t0 = 1 - t1;
			d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
						s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
		END_FOR
		CUDA::set_bnd(N, b, d);
	}

	void project(uint32_t N, float *u, float *v, float *p, float *div)
	{
		uint32_t i, j;
		const float mul = -0.5f/(float)N;
		FOR_EACH_CELL
			const float sub1 = u[IX(i + 1, j)] - u[IX(i - 1, j)];
			const float sub2 = v[IX(i, j + 1)] - v[IX(i, j - 1)];
			const float sum = sub1 + sub2;
			const float result = sum * mul;
		 	div[IX(i, j)] = result;
			p[IX(i, j)] = 0;
		END_FOR
		CUDA::set_bnd(N, 0, div);
		CUDA::set_bnd(N, 0, p);

		CUDA::lin_solve(N, 0, p, div, 1, 4);

		FOR_EACH_CELL
			u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
			v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
		END_FOR
		CUDA::set_bnd(N, 1, u);
		CUDA::set_bnd(N, 2, v);
	}

	void dens_step(uint32_t N, float *x, float *x0, float *u, float *v, float diff, float dt)
	{
		CUDA::add_source(N, x, x0, dt);
		SWAP(x0, x);
		CUDA::diffuse(N, 0, x, x0, diff, dt);
		SWAP(x0, x);
		CUDA::advect(N, 0, x, x0, u, v, dt);
	}
	
	void vel_step(uint32_t N, float *u, float *v, float *u0, float *v0, float visc, float dt)
	{
		CUDA::add_source(N, u, u0, dt);
		SWAP(u0, u);
		CUDA::diffuse(N, 1, u, u0, visc, dt);
		CUDA::add_source(N, v, v0, dt);				
		SWAP(v0, v);
		CUDA::diffuse(N, 2, v, v0, visc, dt);
		CUDA::project(N, u, v, u0, v0);
		SWAP(u0, u);
		SWAP(v0, v);
		CUDA::advect(N, 1, u, u0, u0, v0, dt);
		CUDA::advect(N, 2, v, v0, u0, v0, dt);
		CUDA::project(N, u, v, u0, v0);
	}
}