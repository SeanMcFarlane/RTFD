/*
  ======================================================================
   demo_prof.cpp - SFML free version of the demo that requires no extra libraries
  ----------------------------------------------------------------------
   Author : Sean Mcfarlane
   Creation Date : Jan 2019

   Description:

    This project was created based on the paper "Real-Time Fluid Dynamics for Games" by Jos Stam.
    Code is a modified version of the sample code provided alongside the paper, which was adapted to use C++ and SFML rather than C and OpenGL.
    This was created for the purpose of testing optimization methods as part of my CSC485C coursework.

  =======================================================================
*/

#include <stdlib.h>
#include <stdio.h>
#include "project.h"

/* global variables */

static float dt, diff, visc;
static float force, source;

//static float * u, * v, * u_prev, * v_prev;
//static float * dens, * dens_prev;

float __attribute__ ((aligned(16))) *u, *v, *u_prev, *v_prev;
float __attribute__ ((aligned(16))) *dens, *dens_prev;
float __attribute__ ((aligned(16))) *test;

static uint32_t win_id;
static int win_x, win_y;
static uint32_t mouse_down[3];
static int omx, omy, mx, my;

float timeSpeed;
const uint32_t frameRate = 60;
const float visLen = 30.0f;
const uint32_t resolution = 512;
uint32_t iterations;
uint32_t cur_iter = 0;
uint32_t N;
uint32_t bnd;

bool profiling;
int optim_mode;

/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/

static void clear_data ( void )
{
	uint32_t i, size=(N+bnd)*(N+bnd);

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
	}
	DPRINT("Data cleared\n");
}

static uint32_t allocate_data_simd ( void )
{
	DPRINT("ALLOCATING DATA SIMD\n");
	uint32_t size = (N+bnd)*(N+bnd);

	u 			= (float*)_mm_malloc(size*(sizeof(float)), 16);
	v 			= (float*)_mm_malloc(size*(sizeof(float)), 16);
	u_prev 		= (float*)_mm_malloc(size*(sizeof(float)), 16);
	v_prev 		= (float*)_mm_malloc(size*(sizeof(float)), 16);
	dens 		= (float*)_mm_malloc(size*(sizeof(float)), 16);
	dens_prev 	= (float*)_mm_malloc(size*(sizeof(float)), 16);
	test 		= (float*)_mm_malloc(size*(sizeof(float)), 16);

	if ( !u || !v || !u_prev || !v_prev || !dens || !dens_prev ) {
		fprintf ( stderr, "cannot allocate data\n" );
		return ( 0 );
	}

	DPRINT("ALLOCATED DATA SIMD SUCCESSFULLY\n");

	return ( 1 );
}

static void clear_data_simd ( void )
{
	uint32_t i, size=(N+bnd)*(N+bnd);

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
	}
	DPRINT("Data cleared\n");
}

static void ProfileLoop()
{
	while(cur_iter < iterations){
		cur_iter++;
		switch(optim_mode){
			case 0: // SIMD (IX)
			{
                // Constant force is added each iteration to demonstrate simulation without needing input.
				SIMD::add_force(N / 4, N / 4, 100.0f, 0); 
				SIMD::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				SIMD::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				SIMD::dens_step(N, dens, dens_prev, u, v, diff, dt);
				break;
			}
			case 1:
			{
				// Constant force is added each iteration to demonstrate simulation without needing input.
				base::add_force(N / 4, N / 4, 100.0f, 0); 
				base::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				base::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				base::dens_step(N, dens, dens_prev, u, v, diff, dt);
				break;
			}
			case 2:
			{
				// Constant force is added each iteration to demonstrate simulation without needing input.
				opt::add_force(N / 4, N / 4, 100.0f, 0); 
				opt::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				opt::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				opt::dens_step(N, dens, dens_prev, u, v, diff, dt);
				break;
			}
			case 3:
			{
				// Constant force is added each iteration to demonstrate simulation without needing input.
				parallel::add_force(N / 4, N / 4, 100.0f, 0); 
				parallel::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				parallel::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				parallel::dens_step(N, dens, dens_prev, u, v, diff, dt);
				break;
			}
			case 4: // SIMD again
			{
                // Constant force is added each iteration to demonstrate simulation without needing input.
				SIMD::add_force(N / 4, N / 4, 100.0f, 0); 
				SIMD::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
				SIMD::vel_step(N, u, v, u_prev, v_prev, visc, dt);
				SIMD::dens_step(N, dens, dens_prev, u, v, diff, dt);
				break;
			}
		}
	}

}

namespace opt
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		if (i < 1 || i > N || j < 1 || j > N)
			return;
		u[ZIX(i, j)] += dt * xForce;
		v[ZIX(i, j)] += dt * yForce;
	}
} // namespace opt

namespace parallel
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		if (i < 1 || i > N || j < 1 || j > N)
			return;
		u[ZIX(i, j)] += dt * xForce;
		v[ZIX(i, j)] += dt * yForce;
	}
} // namespace opt


namespace base
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		if (i < 1 || i > N || j < 1 || j > N)
			return;
		u[IX(i, j)] += dt * xForce;
		v[IX(i, j)] += dt * yForce;
	}
} // namespace base

namespace SIMD
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		if (i < 1 || i > N || j < 1 || j > N)
			return;
		u[IX(i, j)] += dt * xForce;
		v[IX(i, j)] += dt * yForce;
	}
} // namespace SIMD

/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main ( int argc, char ** argv )
{

	if ( argc != 1 && argc != 5 ) {
		fprintf ( stderr, "Usage: demo.exe <profiling mode[0-1]> <select implementation[0-4]> <resolution[int]> <iterations[int]>\n");
		return 1;
	}

	if ( argc == 1 ) {
		profiling = true;
		optim_mode = 0;
		N = 256;
		iterations = 1000;
		fprintf ( stderr, "Using defaults: profiler mode, SIMD, 256x256, 1000 iterations\n");
	} else {
		profiling = true;
		optim_mode = atoi(argv[2]);
		N = atoi(argv[3]);
		iterations = atoi(argv[4]);
	}

	if(optim_mode==0||optim_mode==4){bnd = 8;}
	else{bnd = 2;}
	N-=bnd;
	
	timeSpeed = 1.0f;
	diff = 0.0f;
	visc = 0.0f;
	force = 1000.0f;
	source = 100.0f;

    printf( "Beginning test...\n" );

	if (!allocate_data_simd()) return 1;
	clear_data();

	win_x = resolution;
	win_y = resolution;

	if (profiling)
	{
		dt = timeSpeed / 60.0f;
	}
	else{
		dt = timeSpeed / (float)frameRate;
	}

    ProfileLoop();
    
    printf( "Test complete.\n" );
	return 0;
}