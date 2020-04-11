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

	if ( argc != 1 && argc != 4 && argc !=6) {
		fprintf ( stderr, "Usage: \tdemo.exe <select implementation[0-5]> <resolution[int]> <iterations[int]>\n or:\tdemo.exe 5 <resolution[int]> <iterations[int]> <CUDA block size[int]> <CUDA cells per thread[int]>");
		return 1;
	}
	profiling = true;
	if ( argc == 1 ) {
		optim_mode = 4;
		N = 256;
		iterations = 1000;
		printf ("Using defaults: profiling mode, SIMD+Multithreading, 256x256, 1000 iterations\n");
	} else if(argc == 6){
		optim_mode = atoi(argv[1]);
		N = atoi(argv[2]);
		iterations = atoi(argv[3]);

		if(optim_mode != 5){
			fprintf ( stderr, "Using too many arguments for a non-CUDA mode.\n");
			return 1;
		}

		CUDA::BSIZE = atoi(argv[4]);
		if(atoi(argv[5]) <= CUDA::BSIZE){
			CUDA::CELLSPERTHREAD = atoi(argv[5]);
		}
		else{
			fprintf ( stderr, "Cells per thread cannot exceed block size.\n");
			return 1;
		}
	} else {
		optim_mode = atoi(argv[1]);
		N = atoi(argv[2]);
		iterations = atoi(argv[3]);
		CUDA::BSIZE = 32;
		CUDA::CELLSPERTHREAD = 4;
	}

	if(optim_mode==3||optim_mode==4){bnd = 8;}
	else{bnd = 2;}
	N-=bnd;
	pad = bnd/2;
	dim = N + bnd;
	array_size = dim * dim;
	timeSpeed = 1.0f;
	diff = 0.00001f;
	visc = 0.00001f;
	force = 50.0f;
	source = 512.0f;

    printf( "Beginning test...\n" );

	if (optim_mode == 5) {
		if (!CUDA::allocate_data())
		{
			fprintf(stderr, "Cuda allocation failed.\n");
			return 1;
		}
		CUDA::init_cuda_globals(N, dim, bnd, pad);
	}
	else {
		if (!allocate_data_simd()) return 1;
	}

	clear_data();


	printf( "Data allocated...\n" );

	dt = timeSpeed / 60.0f;

	while(cur_iter < iterations){
		cur_iter++;
    	Simulate(optim_mode);
	}
    
	if(optim_mode==5){
		CUDA::dealloc_cuda_globals();
	}

    printf( "Test complete.\n" );
	return 0;
}