#include "project.h"

namespace SIMD { // SIMD implementation

	#ifdef DEBUG_LOG
	void printm128(const __m128 val){
		float __attribute__ ((aligned(16))) *output = (float*)_mm_malloc(4*(sizeof(float)), 16);
		_mm_store_ps( output, val);
		for(int n = 0; n < 4; n++)
		{
			DPRINT(output[n]<<", ");
		}
		DPRINT("\n");
	}
	#else
	void printm128(const __m128 val){}
	#endif

	/*
		SIMD 128 index shift functions
	*/

	static const __m128 simd_rshift2(const __m128 _t1, const __m128 _t2)
	{ // 0123 4567
		return _mm_shuffle_ps( _t1, _t2, _MM_SHUFFLE(1,0,3,2)); // 2345		
	}

	static const __m128 simd_rshift1(const __m128 _t1, const __m128 _t2)
	{ // 0123 2345
		__m128 const shift2 = simd_rshift2(_t1, _t2); // 2345
		//DPRINT("RSHIFT2 ");
		//printm128(shift2);
		return _mm_shuffle_ps( _t1, shift2, _MM_SHUFFLE(2,1,2,1)); // 1234 (X, 1, 2, X) (X, 4, 3, X)
	}

	static const __m128 simd_lshift2(const __m128 _t0, const __m128 _t1)
	{ //-4-3-2-1 0123
		return _mm_shuffle_ps( _t0, _t1, _MM_SHUFFLE(1,0,3,2)); // -2,-1,0,1		
	}

	static const __m128 simd_lshift1(const __m128 _t0, const __m128 _t1)
	{ //-4-3-2-1 0123
		__m128 const shift2 = simd_lshift2(_t0, _t1); // -2,-1,0,1
		//DPRINT("LSHIFT2 ");
		//printm128(shift2);
		return _mm_shuffle_ps( shift2, _t1, _MM_SHUFFLE(2,1,2,1)); // -1,0,1,2
	}

	void add_source(uint32_t N, float *x, float *s, float dt)
	{
		uint32_t i, size = (N + 2) * (N + 2);
		for (i = 0; i < size; i++)
			x[i] += dt * s[i];
	}

	void set_bnd(uint32_t N, uint32_t b, float *x)
	{
		uint32_t i;

		for (i = 1; i <= N; i++)
		{
			x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
			x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
			x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
			x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
		}
		x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
		x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
		x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
		x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
	}

	void m128_test(){
		DPRINT("pre-test\n");
		__m128 _t0 = _mm_set_ps(-1, -2, -3, -4);
		__m128 _t1 = _mm_set_ps(3, 2, 1, 0);
		__m128 _t2 = _mm_set_ps(7, 6, 5, 4);
		DPRINT("_t0 ");
		printm128(_t0);
		DPRINT("_t1 ");
		printm128(_t1);
		DPRINT("_t2 ");
		printm128(_t2);
		__m128 const rshift1 = simd_rshift1(_t1, _t2);
		__m128 const lshift1 = simd_lshift1(_t0, _t1);	
		__m128 const rshift2 = simd_rshift2(_t1, _t2); // 2345	
		DPRINT("rshift1 ");
		printm128(rshift1);
		DPRINT("lshift1 ");
		printm128(lshift1);
	}

	/*
		Potential improvements for SIMD solution:
		- SIMD'ify advect and project
	*/

	void lin_solve(uint32_t N, uint32_t b, float *x, float *x0, float a, float c)
	{
		float cInv = 1.0f/c;
		__m128 _a = _mm_set_ps(a,a,a,a);
		__m128 _c = _mm_set_ps(cInv,cInv,cInv,cInv);
		//m128_test();
		uint32_t i, j, k;


		for (k = 0; k < 20; k++)
		{
			//#pragma omp parallel for
			for ( j=1; j<=N; j++ ) { 
				__m128 mid0;
				__m128 mid1 = _mm_load_ps( &( x[IX(0, j)] ) );
				__m128 mid2 = _mm_load_ps( &( x[IX(4, j)] ) );

				for ( i=4; i<=N-4; i+=4 ) {
					DPRINT("i,j=("<<i<<","<<j<<")\n");
					DPRINT("IX(i,j)="<<IX(i,j)<<"\n");

					mid0 = mid1;
					mid1 = mid2;
					mid2 =_mm_load_ps( &( x[IX(i+4, j)] ) );

					__m128 const center = _mm_load_ps( &( x0[IX(i, j)] ) );

					__m128 const left = simd_lshift1(mid0, mid1);
					__m128 const right = simd_rshift1(mid1, mid2);

					__m128 const up = _mm_load_ps( &( x[IX(i, j+1)] ) );
					__m128 const down = _mm_load_ps( &( x[IX(i, j-1)] ) );

					DPRINT("All loads work!\n")

					__m128 const adj1 = _mm_add_ps(up, down);
					__m128 const adj2 = _mm_add_ps(left, right);
					__m128 const adjacents = _mm_add_ps(adj1, adj2);

					__m128 const scaledAdj = _mm_mul_ps(_a, adjacents);
					__m128 const numerator = _mm_add_ps(center, scaledAdj);
					__m128 const result =  _mm_mul_ps(numerator, _c);

					_mm_store_ps( &(x[IX(i, j)]), result);
					
					DPRINT("SIMD lin_solve success.\n");

				}
			}
			SIMD::set_bnd(N, b, x);
		}
	}

	void diffuse(uint32_t N, uint32_t b, float *x, float *x0, float diff, float dt)
	{
		float a = dt * diff * N * N;
		SIMD::lin_solve(N, b, x, x0, a, 1 + 4 * a);
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
		SIMD::set_bnd(N, b, d);
	}

	void project(uint32_t N, float *u, float *v, float *p, float *div)
	{
		uint32_t i, j;
		float mul = -0.5f;
		const __m128 multiplier = _mm_set_ps(mul,mul,mul,mul);
		const __m128 N128 = _mm_set_ps(N,N,N,N);
		const __m128 zeroes = _mm_setzero_ps();
		for ( j=1; j<=N; j++ ) { 

			__m128 u0;
			__m128 u1 = _mm_load_ps( &( u[IX(0, j)] ) );
			__m128 u2 = _mm_load_ps( &( u[IX(4, j)] ) );

			for ( i=4; i<=N-4; i+=4 ) {
				u0=u1;
				u1=u2;
				u2 = _mm_load_ps( &( u[IX(i+4, j)] ) );

				const __m128 vAbove = _mm_load_ps( &( u[IX(i, j+1)] ) );
				const __m128 vBelow = _mm_load_ps( &( u[IX(i, j-1)] ) );
				
				const __m128 ul = simd_lshift1(u0, u1);
				const __m128 ur = simd_rshift1(u1, u2);

				const __m128 uDiff = _mm_sub_ps(ur, ul);
				const __m128 vDiff = _mm_sub_ps(vAbove, vBelow);

				const __m128 sum = _mm_add_ps(uDiff, vDiff);
				const __m128 result = _mm_mul_ps(sum, multiplier);
				const __m128 result2 = _mm_div_ps(result, N128);
				//const __m128 resultF = _mm_floor_ps(result);
				_mm_store_ps( &(div[IX(i, j)]), result2);
				_mm_store_ps( &(p[IX(i, j)]), zeroes);
			}
		}
		
		FOR_EACH_CELL
			div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] + v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
			p[IX(i, j)] = 0;
		END_FOR
		SIMD::set_bnd(N, 0, div);
		SIMD::set_bnd(N, 0, p);

		SIMD::lin_solve(N, 0, p, div, 1, 4);

		FOR_EACH_CELL
			u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
			v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
		END_FOR
		SIMD::set_bnd(N, 1, u);
		SIMD::set_bnd(N, 2, v);
	}

	void dens_step(uint32_t N, float *x, float *x0, float *u, float *v, float diff, float dt)
	{
		SIMD::add_source(N, x, x0, dt);
		SWAP(x0, x);
		SIMD::diffuse(N, 0, x, x0, diff, dt);
		SWAP(x0, x);
		SIMD::advect(N, 0, x, x0, u, v, dt);
	}

	void vel_step(uint32_t N, float *u, float *v, float *u0, float *v0, float visc, float dt)
	{
		SIMD::add_source(N, u, u0, dt);
		SIMD::add_source(N, v, v0, dt);
		SWAP(u0, u);
		SIMD::diffuse(N, 1, u, u0, visc, dt);
		SWAP(v0, v);
		SIMD::diffuse(N, 2, v, v0, visc, dt);
		base::project(N, u, v, u0, v0);
		SWAP(u0, u);
		SWAP(v0, v);
		SIMD::advect(N, 1, u, u0, u0, v0, dt);
		SIMD::advect(N, 2, v, v0, u0, v0, dt);
		SIMD::project(N, u, v, u0, v0);
	}
}