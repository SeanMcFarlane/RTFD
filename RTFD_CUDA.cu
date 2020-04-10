#include "project.h"

#define BSIZE 16
#define JOBSPERTHREAD 4
#define DEVICE_CHECK_OOB if(i < d_pad || i >= d_dim-d_pad || j < d_pad || j >= d_dim-d_pad) {printf("DEVICE BOUNDS ERROR: Cell (%i,%i) OOB at line %d\n", i,j, __LINE__);}
#define INBOUNDS(i,j) (i>=d_pad && i<d_dim-d_pad && j>=d_pad && j<d_dim-d_pad)
#define FOR_JOBS_IN_BOUNDS for (uint32_t jobnum = 0; jobnum < JOBSPERTHREAD; jobnum++) { if INBOUNDS(i,j) 
#define END_JOBS i++; }	
#define DIX(i,j) ((i+d_pad)+(d_N+d_bnd)*(j))


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
	if (code != cudaSuccess) 
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

namespace CUDA { // CUDA implementation

	//device globals
	__device__ uint32_t d_dim, d_pad, d_bnd, d_N;

	//host globals
	dim3 blockSize;
	uint32_t bx;
	uint32_t by;
	dim3 gridSize;
	dim3 threadCount;


	//
	// CUDA Initialization
	//
	__host__
	int allocate_data_cuda_pinned(void** ptr, size_t size)
	{
		if (cudaHostAlloc(ptr, size, cudaHostAllocDefault) == cudaSuccess) {
			return 1;
		}
		else return 0;
	}

	__global__
	void gpu_init_cuda_globals(uint32_t N, uint32_t dim, uint32_t bnd, uint32_t pad)
	{
		d_dim = dim;
		d_pad = pad;
		d_bnd = bnd;
		d_N = N;
	}

	__host__
	void init_cuda_globals(uint32_t N, uint32_t dim, uint32_t bnd, uint32_t pad)
	{
		gpu_init_cuda_globals<<<1,1>>>(N, dim, bnd, pad);

		blockSize = dim3(BSIZE, BSIZE);
		bx = (int)ceil((dim) / blockSize.x);
		by = (int)ceil((dim) / blockSize.y);
		gridSize = dim3(bx, by);
		threadCount = dim3(BSIZE/JOBSPERTHREAD, BSIZE);
	}

	//
	// Device Kernels
	//

	__global__
	void gpu_project_1(float mul, float* u, float* v, float* p, float* div) {
		uint32_t i = blockIdx.x*(JOBSPERTHREAD*blockDim.x)+threadIdx.x*JOBSPERTHREAD;
		uint32_t const j = blockIdx.y*blockDim.y+threadIdx.y;

		FOR_JOBS_IN_BOUNDS { 
			DEVICE_CHECK_OOB
			const float sub1 = u[DIX(i + 1, j)] - u[DIX(i - 1, j)];
			const float sub2 = v[DIX(i, j + 1)] - v[DIX(i, j - 1)];
			const float sum = sub1 + sub2;
			const float result = sum * mul;
			div[DIX(i, j)] = result;
			p[DIX(i, j)] = 0;
		} END_JOBS
	}

	__global__
	void gpu_project_2(float* p, float* u, float* v) {
		uint32_t i = blockIdx.x*(JOBSPERTHREAD*blockDim.x)+threadIdx.x*JOBSPERTHREAD;
		const uint32_t j = blockIdx.y*blockDim.y+threadIdx.y;
		const uint32_t mul = 0.5f * d_N;
		FOR_JOBS_IN_BOUNDS { 
			DEVICE_CHECK_OOB
			const float a = p[DIX(i + 1, j)] - p[DIX(i - 1, j)];
			u[DIX(i, j)] = u[DIX(i, j)] - mul*a;
			const float b = p[DIX(i, j + 1)] - p[DIX(i, j - 1)];
			v[DIX(i, j)] = v[DIX(i, j)] - mul*b;
		} END_JOBS
	}

	__global__
	void gpu_lin_solve(float* x, float* x0, const float a, const float c) {
		uint32_t i = blockIdx.x*(JOBSPERTHREAD*blockDim.x)+threadIdx.x*JOBSPERTHREAD;
		uint32_t const j = blockIdx.y*blockDim.y+threadIdx.y;

		FOR_JOBS_IN_BOUNDS { 
			DEVICE_CHECK_OOB
			const float left = x[DIX(i - 1, j)];
			const float right = x[DIX(i + 1, j)];
			const float above = x[DIX(i, j + 1)];
			const float below = x[DIX(i, j - 1)];
			const float sum = left + right + above + below;
			const float prod = sum * a;
			const float addOrig = prod + x0[DIX(i, j)];
			const float result = addOrig / c;
			x[DIX(i, j)] = result;
		} END_JOBS
	}

	__global__
	void gpu_advect(float *d_d, float *d_d0, float *d_u, float *d_v, const float dt){
		uint32_t i = blockIdx.x*(JOBSPERTHREAD*blockDim.x)+threadIdx.x*JOBSPERTHREAD;
		uint32_t const j = blockIdx.y*blockDim.y+threadIdx.y;

		uint32_t i0, j0, i1, j1;
		float x, y, s0, t0, s1, t1, dt0;

		FOR_JOBS_IN_BOUNDS { 
			DEVICE_CHECK_OOB
			dt0 = dt * d_N;
			x = i - dt0 * d_u[DIX(i, j)];
			y = j - dt0 * d_v[DIX(i, j)];
			if (x < 0.5f)
				x = 0.5f;
			if (x > d_N + 0.5f)
				x = d_N + 0.5f;
			i0 = (int)x;
			i1 = i0 + 1;
			if (y < 0.5f)
				y = 0.5f;
			if (y > d_N + 0.5f)
				y = d_N + 0.5f;
			j0 = (int)y;
			j1 = j0 + 1;
			s1 = x - i0;
			s0 = 1 - s1;
			t1 = y - j0;
			t0 = 1 - t1;
			d_d[DIX(i, j)] = s0 * (t0 * d_d0[DIX(i0, j0)] + t1 * d_d0[DIX(i0, j1)]) +
				s1 * (t0 * d_d0[DIX(i1, j0)] + t1 * d_d0[DIX(i1, j1)]);
		} END_JOBS
	}

	__global__
	void gpu_set_bnd(uint32_t b, float *x)
	{
		uint32_t const n = blockIdx.x*blockDim.x+threadIdx.x;
				
		const uint32_t ub = d_N+d_pad;
		const uint32_t lb = d_pad-1;
		if(n < ub && n > lb){
			x[DIX(lb, n)] = (b == 1) ? -x[DIX(lb+1, n)] : x[DIX(lb+1, n)];	
			x[DIX(ub, n)] = (b == 1) ? -x[DIX(ub-1, n)] : x[DIX(ub-1, n)];
			x[DIX(n, lb)] = (b == 2) ? -x[DIX(n, lb+1)] : x[DIX(n, lb+1)];
			x[DIX(n, ub)] = (b == 2) ? -x[DIX(n, ub-1)] : x[DIX(n, ub-1)];

		}

		//x[DIX(d_pad-1, d_pad-1)] = 0.5f * (x[DIX(d_pad, d_pad-1)] + x[DIX(d_pad-1, d_pad)]);
		//x[DIX(d_pad-1, d_N + d_pad)] = 0.5f * (x[DIX(d_pad, d_N+d_pad)] + x[DIX(d_pad-1, d_N+d_pad-1)]);
		//x[DIX(N+d_pad, d_pad-1)] = 0.5f * (x[DIX(N+d_pad-1,d_pad-1)] + x[DIX(N+d_pad, d_pad)]);
		//x[DIX(N+d_pad, d_N+d_pad)] = 0.5f * (x[DIX(N+d_pad-1, d_N+d_pad)] + x[DIX(N+d_pad, d_N+d_pad-1)]);
	}

	__global__
	void gpu_add_source(float *x, float *s, const float dt){
		uint32_t i = blockIdx.x*(JOBSPERTHREAD*blockDim.x)+threadIdx.x*JOBSPERTHREAD;
		uint32_t const j = blockIdx.y*blockDim.y+threadIdx.y;

		FOR_JOBS_IN_BOUNDS { 
			DEVICE_CHECK_OOB
			x[DIX(i,j)] += dt * s[DIX(i,j)];
		} END_JOBS
	}

	__host__
	void add_source(float *x, float *s, float dt)
	{
		uint32_t i, size = array_size;
		for (i = 0; i < size; i++)
			x[i] += dt * s[i];
	}

	void lin_solve(uint32_t b, float *x, float *x0, float a, float c)
	{
		auto const size = sizeof(float)*array_size;
		float *dev_x, *dev_x0;
		gpuErrchk( cudaMalloc((void**)&dev_x, size));
		gpuErrchk( cudaMalloc((void**)&dev_x0, size));
		gpuErrchk( cudaMemcpy(dev_x, x, size, cudaMemcpyHostToDevice));
		gpuErrchk( cudaMemcpy(dev_x0, x0, size, cudaMemcpyHostToDevice));

		uint32_t k;
		for (k = 0; k < 20; k++)
		{
			CUDA::gpu_lin_solve<<<gridSize, threadCount>>>(dev_x, dev_x0, a, c);
			CUDA::gpu_set_bnd<<<bx, BSIZE>>>(b, dev_x);
		}

		gpuErrchk(cudaMemcpy(x, dev_x, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(x0, dev_x0, size, cudaMemcpyDeviceToHost));

		gpuErrchk(cudaFree(dev_x));
		gpuErrchk(cudaFree(dev_x0));
	}

	void lin_solve_2(uint32_t b, float *d_x, float *d_x0, float a, float c)
	{
		uint32_t k;
		for (k = 0; k < 20; k++)
		{
			CUDA::gpu_lin_solve<<<gridSize, threadCount>>>(d_x, d_x0, a, c);
			CUDA::gpu_set_bnd<<<bx, BSIZE>>>(b, d_x);
		}
	}

	void diffuse(uint32_t b, float *x, float *x0, float diff, float dt)
	{
		float a = dt * diff * N * N;
		CUDA::lin_solve(b, x, x0, a, 1 + 4 * a);
	}

	void diffuse_2(uint32_t b, float *d_x, float *d_x0, float diff, float dt)
	{
		float a = dt * diff * N * N;
		CUDA::lin_solve_2(b, d_x, d_x0, a, 1 + 4 * a);
	}

	void advect(uint32_t b, float *d, float *d0, float *u, float *v, float dt)
	{
		auto const size = sizeof(float)*array_size;

		float *d_d, *d_d0, *d_u, *d_v;			 
		gpuErrchk(cudaMalloc((void**)&d_d, size	));
		gpuErrchk(cudaMalloc((void**)&d_d0, size));
		gpuErrchk(cudaMalloc((void**)&d_u, size	));
		gpuErrchk(cudaMalloc((void**)&d_v, size	));

		gpuErrchk(cudaMemcpy(d_d, d, size, cudaMemcpyHostToDevice	));
		gpuErrchk(cudaMemcpy(d_d0, d0, size, cudaMemcpyHostToDevice	));
		gpuErrchk(cudaMemcpy(d_u, u, size, cudaMemcpyHostToDevice	));
		gpuErrchk(cudaMemcpy(d_v, v, size, cudaMemcpyHostToDevice	));

		CUDA::gpu_advect<<<gridSize, threadCount>>>(d_d, d_d0, d_u, d_v, dt);
		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(b, d_d);

		gpuErrchk(cudaMemcpy(d, d_d, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaFree(d_d));
		gpuErrchk(cudaMemcpy(d0, d_d0, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaFree(d_d0));
		gpuErrchk(cudaMemcpy(u, d_u, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaFree(d_u));
		gpuErrchk(cudaMemcpy(v, d_v, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaFree(d_v));

	}

	void project(float *u, float *v, float *p, float *div){

		auto const size = sizeof(float)*array_size;

		float *d_u, *d_v, *d_p, *d_div;

		gpuErrchk(cudaMalloc((void**)&d_v, size		));
		gpuErrchk(cudaMalloc((void**)&d_p, size		));
		gpuErrchk(cudaMalloc((void**)&d_u, size		));
		gpuErrchk(cudaMalloc((void**)&d_div, size	));

		gpuErrchk(cudaMemcpy(d_u, u, size, cudaMemcpyHostToDevice		));
		gpuErrchk(cudaMemcpy(d_v, v, size, cudaMemcpyHostToDevice		));
		gpuErrchk(cudaMemcpy(d_p, p, size, cudaMemcpyHostToDevice		));
		gpuErrchk(cudaMemcpy(d_div, div, size, cudaMemcpyHostToDevice	));

		const float mul = -0.5f/(float)N;
		CUDA::gpu_project_1<<<gridSize, threadCount>>>(mul, d_u, d_v, d_p, d_div);

		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(0, d_div);
		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(0, d_p);

		uint32_t k;
		for (k = 0; k < 20; k++)
		{
			CUDA::gpu_lin_solve<<<gridSize, threadCount>>>(d_p, d_div, 1, 4);
			CUDA::gpu_set_bnd<<<bx, BSIZE>>>(0, d_p);
		}

		CUDA::gpu_project_2<<<gridSize, threadCount>>>(d_p, d_u, d_v);

		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(1, d_u);
		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(2, d_v);
																		 
		gpuErrchk(cudaMemcpy(p, d_p, size, cudaMemcpyDeviceToHost		));
		gpuErrchk(cudaMemcpy(div, d_div, size, cudaMemcpyDeviceToHost	));
		gpuErrchk(cudaMemcpy(u, d_u, size, cudaMemcpyDeviceToHost		));
		gpuErrchk(cudaMemcpy(v, d_v, size, cudaMemcpyDeviceToHost		));

		gpuErrchk(cudaFree(d_p	));
		gpuErrchk(cudaFree(d_div));
		gpuErrchk(cudaFree(d_u	));
		gpuErrchk(cudaFree(d_v	));
	}

	void project_2(float *d_u, float *d_v, float *d_p, float *d_div)
	{
		const float mul = -0.5f/(float)N;
		CUDA::gpu_project_1<<<gridSize, threadCount>>>(mul, d_u, d_v, d_p, d_div);

		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(0, d_div);
		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(0, d_p);

		uint32_t k;
		for (k = 0; k < 20; k++)
		{
			CUDA::gpu_lin_solve<<<gridSize, threadCount>>>(d_p, d_div, 1, 4);
			CUDA::gpu_set_bnd<<<bx, BSIZE>>>(0, d_p);
		}

		CUDA::gpu_project_2<<<gridSize, threadCount>>>(d_p, d_u, d_v);

		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(1, d_u);
		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(2, d_v);
	}

	void dens_step(float *x, float *x0, float *u, float *v, float diff, float dt)
	{
	
		/*
			Prepare device memory
		*/

		auto const size = sizeof(float)*array_size;
		float *d_x, *d_x0, *d_u, *d_v;	
		gpuErrchk(cudaMalloc((void**)&d_x, size));
		gpuErrchk(cudaMalloc((void**)&d_x0, size));
		gpuErrchk(cudaMalloc((void**)&d_u, size	));
		gpuErrchk(cudaMalloc((void**)&d_v, size	));
		gpuErrchk(cudaMemcpy(d_x, x, size, cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_x0, x0, size, cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_u, u, size, cudaMemcpyHostToDevice	));
		gpuErrchk(cudaMemcpy(d_v, v, size, cudaMemcpyHostToDevice	));

		/*
			Compute solution
		*/
		
		CUDA::gpu_add_source<<<gridSize, threadCount>>>(d_x, d_x0, dt);
		SWAP(d_x0, d_x);
		CUDA::diffuse_2(0, d_x, d_x0, diff, dt);
		SWAP(d_x0, d_x);
		CUDA::gpu_advect<<<gridSize, threadCount>>>(d_x, d_x0, d_u, d_v, dt);
		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(0, d_u);

		/*
			Return result to CPU memory
		*/

		gpuErrchk(cudaMemcpy(x, d_x, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(x0, d_x0, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(u, d_u, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(v, d_v, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaFree(d_x));
		gpuErrchk(cudaFree(d_x0));
		gpuErrchk(cudaFree(d_u));
		gpuErrchk(cudaFree(d_v));


	}
	
	void vel_step(float *u, float *v, float *u0, float *v0, float visc, float dt)
	{

		/*
			Prepare device memory
		*/
		auto const size = sizeof(float)*array_size;
		float *d_u, *d_v, *d_u0, *d_v0; 	
		float *d_v0_temp, *d_u0_temp; 	
		gpuErrchk(cudaMalloc((void**)&d_u, size	));
		gpuErrchk(cudaMalloc((void**)&d_v, size	));
		gpuErrchk(cudaMalloc((void**)&d_u0, size));
		gpuErrchk(cudaMalloc((void**)&d_v0, size));

		gpuErrchk(cudaMalloc((void**)&d_v0_temp, size));
		gpuErrchk(cudaMalloc((void**)&d_u0_temp, size));

		gpuErrchk(cudaMemcpy(d_u, u, size, cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_v, v, size, cudaMemcpyHostToDevice));
		gpuErrchk(cudaMemcpy(d_u0, u0, size, cudaMemcpyHostToDevice	));
		gpuErrchk(cudaMemcpy(d_v0, v0, size, cudaMemcpyHostToDevice	));


		/*
			Compute solution
		*/
		
		CUDA::gpu_add_source<<<gridSize, threadCount>>>(d_u, d_u0, dt);
		SWAP(d_u0, d_u);
		CUDA::diffuse_2(1, d_u, d_u0, visc, dt);
		
		CUDA::gpu_add_source<<<gridSize, threadCount>>>(d_v, d_v0, dt);
		SWAP(d_v0, d_v);
		CUDA::diffuse_2(2, d_v, d_v0, visc, dt);
		
		CUDA::project_2(d_u, d_v, d_u0, d_v0);
		
		SWAP(d_u0, d_u);
		SWAP(d_v0, d_v);

		gpuErrchk(cudaMemcpy(d_u0_temp, d_u0, size, cudaMemcpyDeviceToDevice ));
		gpuErrchk(cudaMemcpy(d_v0_temp, d_v0, size, cudaMemcpyDeviceToDevice ));
		
		CUDA::gpu_advect<<<gridSize, threadCount>>>(d_u, d_u0_temp, d_u0, d_v0, dt);
		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(1, d_u);

		CUDA::gpu_advect<<<gridSize, threadCount>>>(d_v, d_v0_temp, d_u0, d_v0, dt);
		CUDA::gpu_set_bnd<<<bx, BSIZE>>>(2, d_v);
		
		CUDA::project_2(d_u, d_v, d_u0, d_v0);

		/*
			Return result to CPU memory
		*/

		gpuErrchk(cudaMemcpy(u, d_u, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(v, d_v, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(u0, d_u0, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaMemcpy(v0, d_v0, size, cudaMemcpyDeviceToHost));
		gpuErrchk(cudaFree(d_u));
		gpuErrchk(cudaFree(d_v));
		gpuErrchk(cudaFree(d_u0));
		gpuErrchk(cudaFree(d_v0));
		gpuErrchk(cudaFree(d_u0_temp));
		gpuErrchk(cudaFree(d_v0_temp));

		
	}
}

