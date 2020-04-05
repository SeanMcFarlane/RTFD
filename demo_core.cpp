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

/*
  ----------------------------------------------------------------------
	Global Variables
  ----------------------------------------------------------------------
*/

float dt, diff, visc;
float force, source;
float __attribute__ ((aligned(16))) *u, *v, *u_prev, *v_prev;
float __attribute__ ((aligned(16))) *dens, *dens_prev;
float __attribute__ ((aligned(16))) *test;
float timeSpeed;
const uint32_t frameRate = 60;
const uint32_t resolution = 1024;
uint32_t iterations;
uint32_t cur_iter = 0;
uint32_t N;
uint32_t bnd;
uint32_t pad;
bool profiling;
int optim_mode;
const uint32_t zoneLen = 4;
const uint32_t zoneSize = 16;
const uint32_t divShift = 2; //Bit shift amount to perform division
const uint32_t zonesInRow = 64; //Should be equal to N/zoneLen.

/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/

void clear_data(){
	uint32_t i, size=(N+bnd)*(N+bnd);

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
	}
	DPRINT("Data cleared\n");
}

uint32_t allocate_data_simd()
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

// static void clear_data_simd ( void )
// {
// 	uint32_t i, size=(N+bnd)*(N+bnd);

// 	for ( i=0 ; i<size ; i++ ) {
// 		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
// 	}
// 	DPRINT("Data cleared\n");
// }

void Simulate(uint32_t optim_mode)
{
	switch(optim_mode){
		case 0: // Baseline
		{
			// Constant force is added each iteration to demonstrate simulation without needing input.
			base::add_force(N / 4, N / 4, 100.0f, 0); 
			base::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
			base::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			base::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
		case 1: // Singlethreaded Optimized
		{
			// Constant force is added each iteration to demonstrate simulation without needing input.
			opt::add_force(N / 4, N / 4, 100.0f, 0); 
			opt::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
			opt::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			opt::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
		case 2: // Parallelized
		{
			// Constant force is added each iteration to demonstrate simulation without needing input.
			parallel::add_force(N / 4, N / 4, 100.0f, 0); 
			parallel::add_force(3 * N / 4, 3 * N / 4, -100.0f, 0);
			parallel::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			parallel::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
		case 3: // SIMD (IX)
		{
			// Constant force is added each iteration to demonstrate simulation without needing input.
			SIMD::add_force(N / 4 + pad, N / 4 + pad, 100.0f, 0); 
			SIMD::add_force(3 * N / 4 + pad, 3 * N / 4 + pad, -100.0f, 0);
			SIMD::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			SIMD::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
		case 4: // SIMD & Parallelized
		{
			// Constant force is added each iteration to demonstrate simulation without needing input.
			SIMD_PARA::add_force(N / 4 + pad, N / 4 + pad, 100.0f, 0); 
			SIMD_PARA::add_force(3 * N / 4 + pad, 3 * N / 4 + pad, -100.0f, 0);
			SIMD_PARA::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			SIMD_PARA::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
	}
}

namespace opt
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		if (i < pad || i > N || j < pad || j > N)
			return;
		u[ZIX(i, j)] += dt * xForce;
		v[ZIX(i, j)] += dt * yForce;
	}
} // namespace opt

namespace parallel
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		opt::add_force(i,j,xForce,yForce);
	}
} // namespace parallel


namespace base
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
        if (i < pad || i > N || j < pad || j > N)
			return;
		u[IX(i, j)] += dt * xForce;
		v[IX(i, j)] += dt * yForce;
	}
} // namespace base

namespace SIMD
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		base::add_force(i, j, xForce, yForce);
	}
} // namespace SIMD

namespace SIMD_PARA
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		base::add_force(i, j, xForce, yForce);
	}
} // namespace SIMD_PARA