#pragma once

#define _USE_MATH_DEFINES
#include <iostream>
#include <assert.h>
#include <cmath>
#include <omp.h>    // for multi-core parallelism
#include "nmmintrin.h" // for SSE4.2
#include "immintrin.h" // for AVX

#define IX(i,j) ((i)+(N+bnd)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( j=bnd/2; j<N+(bnd/2); j++ ) { for ( i=bnd/2 ; i<N+(bnd/2); i++ ) { 
#define END_FOR }}
#define FOR_EACH_CELL_FULL for ( j=0; j<N+bnd; j++ ) { for ( i=0 ; i<N+bnd; i++ ) { 
//#define DEBUG_LOG
//#define DISPLAY_COORDS
#ifdef DEBUG_LOG 
#define DPRINT(str){std::cout << str;}
#else
#define DPRINT(str)
#endif
#define UNUSED(x) (void)(x)

extern float dt, diff, visc;
extern float force, source;
extern float __attribute__ ((aligned(16))) *u, *v, *u_prev, *v_prev;
extern float __attribute__ ((aligned(16))) *dens, *dens_prev;
extern float __attribute__ ((aligned(16))) *test;
extern float timeSpeed;
extern const uint32_t frameRate;
extern const uint32_t resolution;
extern uint32_t iterations;
extern uint32_t cur_iter;
extern uint32_t N;
extern uint32_t bnd;
extern uint32_t pad;
extern bool profiling;
extern int optim_mode;
extern const uint32_t zoneLen;
extern const uint32_t zoneSize;
extern const uint32_t divShift; //Bit shift amount to perform division
extern const uint32_t zonesInRow; //Should be equal to N/zoneLen.
extern uint32_t array_size;

static const uint32_t ZIX(const uint32_t x, const uint32_t y)
{
    const uint32_t zoneXc = (x >> divShift);
    const uint32_t zoneYc = (y >> divShift);
    const uint32_t zoneIndexc = zoneXc + zoneYc * zonesInRow;
    const uint32_t localXc = (x - (zoneXc << divShift));
    const uint32_t localYc = (y - (zoneYc << divShift));
    #ifdef DEBUG_LOG
    //DPRINT("Result (" << x << ", " << y << ")=[" << result << "]\n");
    const uint32_t result = (zoneIndexc * zoneSize) + localXc + (localYc * zoneLen);
    if(result < 0 || result > (N+bnd)*(N+bnd))
    {
        DPRINT("Result [" << result << "] is outside of range ["<<(N+bnd)*(N+bnd)<<"]!\n");
        DPRINT("\n(x,y)=(" << x << ", " << y << ")\n");
        DPRINT("Zone (" << zoneXc << ", " << zoneYc << ")\n");
        DPRINT("LocalXY (" << localXc << ", " << localYc << ")\n");

        throw "Result index out of range!";
    }
    #endif
    return (zoneIndexc * zoneSize) + localXc + (localYc * zoneLen);
}

namespace parallel{
    void dens_step(uint32_t N, float *x, float *x0, float *u, float *v, float diff, float dt);
    void vel_step(uint32_t N, float *u, float *v, float *u0, float *v0, float visc, float dt);
    void add_source(int i, int j, float density, int diameter);
    void add_force(uint32_t i, uint32_t j, float xForce, float yForce);
    void project(uint32_t N, float *u, float *v, float *p, float *div);
	void render_velocity();
    void set_bnd(uint32_t N, uint32_t b, float *x);
}

namespace opt{
    void dens_step(uint32_t N, float *x, float *x0, float *u, float *v, float diff, float dt);
    void vel_step(uint32_t N, float *u, float *v, float *u0, float *v0, float visc, float dt);
    void add_source(int i, int j, float density, int diameter);
    void add_force(uint32_t i, uint32_t j, float xForce, float yForce);
    void project(uint32_t N, float *u, float *v, float *p, float *div);
	void render_velocity();
    void set_bnd(uint32_t N, uint32_t b, float *x);
}

namespace base{
    void dens_step(uint32_t N, float *x, float *x0, float *u, float *v, float diff, float dt);
    void vel_step(uint32_t N, float *u, float *v, float *u0, float *v0, float visc, float dt);
    void add_source(int i, int j, float density, int diameter);
    void add_force(uint32_t i, uint32_t j, float xForce, float yForce);
    void project(uint32_t N, float *u, float *v, float *p, float *div);
	void render_velocity();
    void set_bnd(uint32_t N, uint32_t b, float *x);
} // namespace base

namespace SIMD{
    void dens_step(uint32_t N, float *x, float *x0, float *u, float *v, float diff, float dt);
    void vel_step(uint32_t N, float *u, float *v, float *u0, float *v0, float visc, float dt);
    void add_source(int i, int j, float density, int diameter);
    void add_force(uint32_t i, uint32_t j, float xForce, float yForce);
    void project(uint32_t N, float *u, float *v, float *p, float *div);
	void render_velocity();
    void set_bnd(uint32_t N, uint32_t b, float *x);
    void m128_test();
    void m128_test2(float *test);
} // namespace SIMD

namespace SIMD_PARA{
    void dens_step(uint32_t N, float *x, float *x0, float *u, float *v, float diff, float dt);
    void vel_step(uint32_t N, float *u, float *v, float *u0, float *v0, float visc, float dt);
    void add_source(int i, int j, float density, int diameter);
    void add_force(uint32_t i, uint32_t j, float xForce, float yForce);
    void project(uint32_t N, float *u, float *v, float *p, float *div);
	void render_velocity();
    void set_bnd(uint32_t N, uint32_t b, float *x);
    void m128_test();
    void m128_test2(float *test);
} // namespace SIMD_PARA

void Simulate(uint32_t optim_mode);
void clear_data();
uint32_t allocate_data_simd();
