#!/bin/bash

g++ -c -fopenmp -mavx -O3 -g demo_core.cpp
g++ -c -fopenmp -mavx -O3 -g demo_perf.cpp
g++ -c -fopenmp -mavx -O3 -g demo_render.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_Base.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_Opt.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_Parallel.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_SIMD.cpp
g++ -fopenmp -mavx -O3 -g demo_core.o demo_render.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o -o demo_render.exe -lsfml-graphics -lsfml-window -lsfml-system
g++ -fopenmp -mavx -O3 -g demo_core.o demo_perf.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o -o demo_perf.exe
rm *.o
echo "Build complete ✓"