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

//SIMD ARRAYS
//float __attribute__ ((aligned(16))) *u, *v, *u_prev, *v_prev;
//float __attribute__ ((aligned(16))) *dens, *dens_prev;
//float __attribute__ ((aligned(16))) *test;

//CUDA ARRAYS
float *u, *v, *u_prev, *v_prev;
float *dens, *dens_prev;
float *test;


float timeSpeed;
const uint32_t frameRate = 120;
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
uint32_t array_size;

/*
  ----------------------------------------------------------------------
   free/clear/allocate simulation data
  ----------------------------------------------------------------------
*/

void clear_data(){
	uint32_t i, size=array_size;

	for ( i=0 ; i<size ; i++ ) {
		u[i] = v[i] = u_prev[i] = v_prev[i] = dens[i] = dens_prev[i] = 0.0f;
	}
	DPRINT("Data cleared\n");
}

int allocate_data_simd()
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

int allocate_data(void)
{
	int size = (N + bnd) * (N + bnd);

	u = (float*)malloc(size * sizeof(float));
	v = (float*)malloc(size * sizeof(float));
	u_prev = (float*)malloc(size * sizeof(float));
	v_prev = (float*)malloc(size * sizeof(float));
	dens = (float*)malloc(size * sizeof(float));
	dens_prev = (float*)malloc(size * sizeof(float));

	if (!u || !v || !u_prev || !v_prev || !dens || !dens_prev) {
		fprintf(stderr, "cannot allocate data\n");
		return (0);
	}

	return (1);
}

void Simulate(uint32_t optim_mode)
{
	switch(optim_mode){
		case 0: // Baseline
		{
			// Constant smoke sources added at two points in the field.
			base::add_density(N / 4, N / 4, 1024, 4);
			base::add_density((3 * N / 4), (3 * N / 4), 1024, 4);
			// Constant force is added each iteration to demonstrate turbulence without needing input.
			base::add_force(N / 4, N / 4, 75.0f, 75.0f); 
			base::add_force(3 * N / 4, 3 * N / 4, -75.0f, -75.0f);
			base::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			base::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
		case 1: // Singlethreaded Optimized
		{
			// Constant smoke sources are added at two points in the field.
			opt::add_density(N / 4, N / 4, 1024, 4);
			opt::add_density((3 * N / 4), (3 * N / 4), 1024, 4);
			// Constant force is added each iteration to demonstrate turbulence without needing input.
			opt::add_force(N / 4, N / 4, 75.0f, 75.0f);
			opt::add_force(3 * N / 4, 3 * N / 4, -75.0f, -75.0f);
			opt::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			opt::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
		case 2: // Parallelized
		{
			// Constant smoke sources are added at two points in the field.
			parallel::add_density(N / 4, N / 4, 1024, 4);
			parallel::add_density((3 * N / 4), (3 * N / 4), 1024, 4);
			// Constant force is added each iteration to demonstrate turbulence without needing input.
			parallel::add_force(N / 4, N / 4, 75.0f, 75.0f);
			parallel::add_force(3 * N / 4, 3 * N / 4, -75.0f, -75.0f);
			parallel::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			parallel::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
		case 3: // SIMD (IX)
		{
			// Constant smoke sources are added at two points in the field.
			SIMD::add_density(N / 4, N / 4, 1024, 4);
			SIMD::add_density((3 * N / 4), (3 * N / 4), 1024, 4);
			// Constant force is added each iteration to demonstrate turbulence without needing input.
			SIMD::add_force(N / 4 + pad, N / 4 + pad, 75.0f, 75.0f);
			SIMD::add_force(3 * N / 4 + pad, 3 * N / 4 + pad, -75.0f, -75.0f);
			SIMD::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			SIMD::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
		case 4: // SIMD & Parallelized
		{
			// Constant smoke sources are added at two points in the field.
			SIMD_PARA::add_density(N / 4, N / 4, 1024, 4);
			SIMD_PARA::add_density((3 * N / 4), (3 * N / 4), 1024, 4);
			// Constant force is added each iteration to demonstrate turbulence without needing input.
			SIMD_PARA::add_force(N / 4 + pad, N / 4 + pad, 75.0f, 75.0f);
			SIMD_PARA::add_force(3 * N / 4 + pad, 3 * N / 4 + pad, -75.0f, -75.0f);
			SIMD_PARA::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			SIMD_PARA::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
		case 5: // CUDA
		{
			// Constant smoke sources are added at two points in the field.
			CUDA::add_density(N / 4, N / 4, 1024, 4);
			CUDA::add_density((3 * N / 4), (3 * N / 4), 1024, 4);
			// Constant force is added each iteration to demonstrate turbulence without needing input.
			CUDA::add_force(N / 4 + pad, N / 4 + pad, 150.0f, 150.0f);
			CUDA::add_force(3 * N / 4 + pad, 3 * N / 4 + pad, -150.0f, -150.0f);
			CUDA::vel_step(N, u, v, u_prev, v_prev, visc, dt);
			CUDA::dens_step(N, dens, dens_prev, u, v, diff, dt);
			break;
		}
	}
}

//
// Row-by-row indexed functions
// 

namespace base
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
        if (i < pad || i > N || j < pad || j > N)
			return;
		u[IX(i, j)] += dt * xForce;
		v[IX(i, j)] += dt * yForce;
	}
	void add_density(int i, int j, float density, int diameter) {
		if (i < pad || i > N || j < pad || j > N) {
			return;
		}
		if (diameter > 1) {
			int iter = 0;
			for (int i_offset = -diameter / 2; i_offset <= diameter / 2; i_offset++) {
				for (int j_offset = -diameter / 2; j_offset <= diameter / 2; j_offset++) {
					if (i_offset * i_offset + j_offset * j_offset > (diameter / 2) * (diameter / 2)) { continue; }
					uint32_t index = IX(i + i_offset, j + j_offset);
					if (index < array_size && index >= 0) {
						dens_prev[IX(i+i_offset, j+j_offset)] += density;
					}
				}
			}
		}
	}
} // namespace base

namespace SIMD
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		base::add_force(i, j, xForce, yForce);
	}
	void add_density(int i, int j, float density, int diameter) {
		base::add_density(i, j, density, diameter);
	}
} // namespace SIMD

namespace SIMD_PARA
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce){
		base::add_force(i, j, xForce, yForce);
	}
	void add_density(int i, int j, float density, int diameter) {
		base::add_density(i, j, density, diameter);
	}
} // namespace SIMD_PARA

namespace CUDA
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce) {
		base::add_force(i, j, xForce, yForce);
	}
	void add_density(int i, int j, float density, int diameter) {
		base::add_density(i, j, density, diameter);
	}
} // namespace CUDA

//
// Zone-indexed functions
//

namespace opt
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce) {
		if (i < pad || i > N || j < pad || j > N)
			return;
		u[ZIX(i, j)] += dt * xForce;
		v[ZIX(i, j)] += dt * yForce;
	}

	void add_density(int i, int j, float density, int diameter) {
		if (i < pad || i > N || j < pad || j > N)
			return;

		if (diameter > 1) {
			for (int i_offset = -diameter / 2; i_offset <= diameter / 2; i_offset++) {
				for (int j_offset = -diameter / 2; j_offset <= diameter / 2; j_offset++) {
					if (i_offset * i_offset + j_offset * j_offset > (diameter / 2) * (diameter / 2)) { continue; }
					uint32_t index = ZIX(i + i_offset, j + j_offset);
					if (index < array_size && index >= 0) {
						dens_prev[ZIX(i + i_offset, j + j_offset)] += density;
					}
				}
			}
		}
	}
} // namespace opt

namespace parallel
{
	void add_force(uint32_t i, uint32_t j, float xForce, float yForce) {
		opt::add_force(i, j, xForce, yForce);
	}
	void add_density(int i, int j, float density, int diameter) {
		opt::add_density(i, j, density, diameter);
	}
} // namespace parallel
