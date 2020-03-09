#pragma once

#define _USE_MATH_DEFINES
#include <SFML/Graphics.hpp>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <omp.h>    // for multi-core parallelism

#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N; i++ ) { for ( j=1 ; j<=N; j++ ) {
#define END_FOR }}
//#define DEBUG_LOG
//#define DISPLAY_COORDS
#ifdef DEBUG_LOG 
#define DPRINT(str){std::cout << str;}
#else
#define DPRINT(str)
#endif

extern sf::RenderWindow *window;
extern sf::Clock *g_clock;
extern float timeSpeed;
extern const uint32_t frameRate;
extern const float visLen;
extern const uint32_t resolution;
extern uint32_t iterations;
extern uint32_t N; //Note to self: Making N const does not improve performance.

extern bool profiling;
extern int optim_mode;
const uint32_t zoneLen = 4;
const uint32_t zoneSize = 16;
const uint32_t divShift = 2; //Bit shift amount to perform division
const uint32_t zonesInRow = 16; //Should be equal to N/zoneLen.

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
    if(result < 0 || result > (N+2)*(N+2))
    {
        DPRINT("Result [" << result << "] is outside of range ["<<(N+2)*(N+2)<<"]!\n");
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
    void add_force(uint32_t i, uint32_t j, float xForce, float yForce);
	void render_velocity();
}

namespace opt{
    void dens_step(uint32_t N, float *x, float *x0, float *u, float *v, float diff, float dt);
    void vel_step(uint32_t N, float *u, float *v, float *u0, float *v0, float visc, float dt);
    void add_force(uint32_t i, uint32_t j, float xForce, float yForce);
	void render_velocity();
}

namespace base{
    void dens_step(uint32_t N, float *x, float *x0, float *u, float *v, float diff, float dt);
    void vel_step(uint32_t N, float *u, float *v, float *u0, float *v0, float visc, float dt);
    void add_force(uint32_t i, uint32_t j, float xForce, float yForce);
	void render_velocity();
} // namespace opt