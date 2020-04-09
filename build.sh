#!/bin/bash

g++ -c -fopenmp -mavx -O3 -g demo_core.cpp
g++ -c -fopenmp -mavx -O3 -g demo_perf.cpp
g++ -c -fopenmp -mavx -O3 -g demo_render.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_Base.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_Opt.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_Parallel.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_SIMD.cpp
nvcc -c RTFD_CUDA.cu -o RTFD_CUDA.obj
g++ -fopenmp -mavx -O3 -g demo_core.o demo_render.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o RTFD_CUDA.obj -o ./Linux/demo_render.exe -lsfml-graphics -lsfml-window -lsfml-system -L/usr/local/cuda/lib64 -lcuda -lcudart
g++ -fopenmp -mavx -O3 -g demo_core.o demo_perf.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o RTFD_CUDA.obj -o ./Linux/demo_perf.exe -L/usr/local/cuda/lib64 -lcuda -lcudart
rm *.o
rm *.obj
echo "Build complete ✓"