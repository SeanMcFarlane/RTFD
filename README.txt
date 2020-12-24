Demonstration: https://www.youtube.com/watch?v=6nD7S7YMc8E
This project was created based on the paper "Real-Time Fluid Dynamics for Games" by Jos Stam.[1]
Code is a modified version of the sample code provided alongside the paper, which was adapted to use C++ and SFML rather than C and OpenGL.
This was created for the purpose of testing optimization methods as part of my CSC485C coursework.

Core algorithms for the simulation can be found in RTFD_[VERSION].cpp files 
Experiment setup code can be found in demo.cpp
Global variables and indexing functions can be found in project.h

LAUNCH OPTIONS: ./demo_render.exe [profiling mode (0,1)] [implementation (0,1,2,3)] [resolution] [iterations]
Implementation options are as follows:
0 - Base
1 - Single-Core Optimized
2 - Multithreaded Optimized
3 - SIMD
4 - SIMD+Multithreaded
5 - CUDA

REQUIRED PACKAGES: 
libsfml-dev

This project is compatible with both 64-bit Windows and Linux systems.

An SFML-free executable is included in the project, useful when you don't need any GUI and just want to profile 
performance. To use this version, execute demo_prof.exe instead.

[1] Stam, Jos. (2003). Real-Time Fluid Dynamics for Games. 
https://www.researchgate.net/publication/2560062_Real-Time_Fluid_Dynamics_for_Games

Sean David Mcfarlane - seanmcfarlane115@gmail.com
