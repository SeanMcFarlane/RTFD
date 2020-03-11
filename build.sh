#!/bin/bash

g++ -c -fopenmp -g -O2 demo.cpp
g++ -c -fopenmp -g -O2 RTFD_Base.cpp
g++ -c -fopenmp -g -O2 RTFD_Opt.cpp
g++ -c -fopenmp -g -O2 RTFD_Parallel.cpp
g++ -c -fopenmp -g -O2 RTFD_SIMD.cpp
g++ -fopenmp -g -O2 demo.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o -o demo.exe -lsfml-graphics -lsfml-window -lsfml-system
