#include <stdlib.h>
#include <stdio.h>
#include "project.h"

/*
  ----------------------------------------------------------------------
   GUI-Free Main --- main routine
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

	if(optim_mode==3||optim_mode==4){bnd = 8;}
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

	dt = timeSpeed / 60.0f;

	while(cur_iter < iterations){
		cur_iter++;
    	Simulate(optim_mode);
	}
    
    printf( "Test complete.\n" );
	return 0;
}