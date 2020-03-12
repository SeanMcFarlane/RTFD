#!/bin/bash

g++ -c -fopenmp -mavx -O3 -g demo.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_Base.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_Opt.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_Parallel.cpp
g++ -c -fopenmp -mavx -O3 -g RTFD_SIMD.cpp
g++ -fopenmp -mavx -O3 -g demo.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o -o demo.exe -lsfml-graphics -lsfml-window -lsfml-system
rm *.o

#x86_64-w64-mingw32-g++ -c -static -static-libgcc -static-libstdc++ -DSFML_STATIC -fopenmp -mavx -g demo.cpp
#x86_64-w64-mingw32-g++ -c -static -static-libgcc -static-libstdc++ -DSFML_STATIC -fopenmp -mavx -g RTFD_Base.cpp
#x86_64-w64-mingw32-g++ -c -static -static-libgcc -static-libstdc++ -DSFML_STATIC -fopenmp -mavx -g RTFD_Opt.cpp
#x86_64-w64-mingw32-g++ -c -static -static-libgcc -static-libstdc++ -DSFML_STATIC -fopenmp -mavx -g RTFD_Parallel.cpp
#x86_64-w64-mingw32-g++ -c -static -static-libgcc -static-libstdc++ -DSFML_STATIC -fopenmp -mavx -g RTFD_SIMD.cpp
#x86_64-w64-mingw32-g++ -static -static-libgcc -static-libstdc++ -DSFML_STATIC -fopenmp -mavx -g demo.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o -o Windows/demo_windows.exe -lsfml-graphics -lsfml-window -lsfml-system
#rm *.o