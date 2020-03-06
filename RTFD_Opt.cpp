#include "project.h"

// Single core optimized functions for Interim Report 1

namespace opt { 
	void add_source(uint32_t N, float *x, float *s, float dt)
	{
		uint32_t i, size=(N+2)*(N+2);
		for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
	}

	void set_bnd ( uint32_t N, uint32_t b, float * x )
	{
		uint32_t i;
		//DPRINT("set_bnd begin\n");
		for ( i=1 ; i<=N ; i++ ) {
			x[ZIX(0  ,i)] = b==1 ? -x[ZIX(1,i)] : x[ZIX(1,i)];
			x[ZIX(N+1,i)] = b==1 ? -x[ZIX(N,i)] : x[ZIX(N,i)];
			x[ZIX(i,0  )] = b==2 ? -x[ZIX(i,1)] : x[ZIX(i,1)];
			x[ZIX(i,N+1)] = b==2 ? -x[ZIX(i,N)] : x[ZIX(i,N)];
		}
		x[ZIX(0  ,0  )] = 0.5f*(x[ZIX(1,0  )]+x[ZIX(0  ,1)]);
		x[ZIX(0  ,N+1)] = 0.5f*(x[ZIX(1,N+1)]+x[ZIX(0  ,N)]);
		x[ZIX(N+1,0  )] = 0.5f*(x[ZIX(N,0  )]+x[ZIX(N+1,1)]);
		x[ZIX(N+1,N+1)] = 0.5f*(x[ZIX(N,N+1)]+x[ZIX(N+1,N)]);
		//DPRINT("set_bnd end\n");
	}

	void lin_solve ( uint32_t N, uint32_t b, float * x, float * x0, float a, float c )
	{
		//DPRINT("lin_solve begin\n");
		uint32_t i, j, k, zX, zY;
		for ( k=0 ; k<20 ; k++ ) { 
			for( zY=0; zY<zonesInRow; zY++){
				for( zX=0; zX<zonesInRow; zX++){
					for ( j=zY*zoneLen; j<(zY+1)*zoneLen; j++ ) {
						for ( i=zX*zoneLen; i<(zX+1)*zoneLen; i++ ) { 
							if(i <= 0 || j <= 0 || i > N || j > N){continue;}
							x[ZIX(i, j)] = (x0[ZIX(i, j)] + a * (x[ZIX(i - 1, j)] + x[ZIX(i + 1, j)] + x[ZIX(i, j - 1)] + x[ZIX(i, j + 1)])) / c;
						}
					}
				}
			}
			// FOR_EACH_CELL
			// 	x[ZIX(i, j)] = (x0[ZIX(i, j)] + a * (x[ZIX(i - 1, j)] + x[ZIX(i + 1, j)] + x[ZIX(i, j - 1)] + x[ZIX(i, j + 1)])) / c;
			// END_FOR
			opt::set_bnd ( N, b, x );
		}
		//DPRINT("lin_solve end\n");
	}

	void diffuse ( uint32_t N, uint32_t b, float * x, float * x0, float diff, float dt )
	{
		float a=dt*diff*N*N;
		opt::lin_solve ( N, b, x, x0, a, 1+4*a );
	}

	void advect ( uint32_t N, uint32_t b, float * d, float * d0, float * u, float * v, float dt )
	{
		uint32_t i, j, i0, j0, i1, j1;
		float x, y, s0, t0, s1, t1, dt0;

		dt0 = dt*N;
		FOR_EACH_CELL
			x = i-dt0*u[ZIX(i,j)]; y = j-dt0*v[ZIX(i,j)];
			if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
			if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
			s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
			d[ZIX(i,j)] = s0*(t0*d0[ZIX(i0,j0)]+t1*d0[ZIX(i0,j1)])+
						s1*(t0*d0[ZIX(i1,j0)]+t1*d0[ZIX(i1,j1)]);
		END_FOR
		opt::set_bnd ( N, b, d );
	}

	void project ( uint32_t N, float * u, float * v, float * p, float * div )
	{
		uint32_t i, j;

		FOR_EACH_CELL
			div[ZIX(i,j)] = -0.5f*(u[ZIX(i+1,j)]-u[ZIX(i-1,j)]+v[ZIX(i,j+1)]-v[ZIX(i,j-1)])/N;
			p[ZIX(i,j)] = 0;
		END_FOR	
		opt::set_bnd ( N, 0, div ); opt::set_bnd ( N, 0, p );

		opt::lin_solve ( N, 0, p, div, 1, 4 );

		FOR_EACH_CELL
			u[ZIX(i,j)] -= 0.5f*N*(p[ZIX(i+1,j)]-p[ZIX(i-1,j)]);
			v[ZIX(i,j)] -= 0.5f*N*(p[ZIX(i,j+1)]-p[ZIX(i,j-1)]);
		END_FOR
		opt::set_bnd ( N, 1, u ); opt::set_bnd ( N, 2, v );
	}

	void dens_step(uint32_t N, float *x, float *x0, float *u, float *v, float diff, float dt)
	{
		opt::add_source(N, x, x0, dt);
		SWAP(x0, x);
		opt::diffuse(N, 0, x, x0, diff, dt);
		SWAP(x0, x);
		opt::advect(N, 0, x, x0, u, v, dt);
	}

	void vel_step ( uint32_t N, float * u, float * v, float * u0, float * v0, float visc, float dt )
	{
		opt::add_source ( N, u, u0, dt ); opt::add_source ( N, v, v0, dt );
		SWAP ( u0, u ); opt::diffuse ( N, 1, u, u0, visc, dt );
		SWAP ( v0, v ); opt::diffuse ( N, 2, v, v0, visc, dt );
		opt::project ( N, u, v, u0, v0 );
		SWAP ( u0, u ); SWAP ( v0, v );
		opt::advect ( N, 1, u, u0, u0, v0, dt ); opt::advect ( N, 2, v, v0, u0, v0, dt );
		opt::project ( N, u, v, u0, v0 );
	}
}
