ECHO BEGIN COMPILE
g++ -c -fopenmp -mavx -O3 -ggdb demo_core.cpp
g++ -c -fopenmp -mavx -O3 -ggdb demo_perf.cpp
g++ -c -fopenmp -mavx -O3 -ggdb demo_render.cpp
g++ -c -fopenmp -mavx -O3 -ggdb RTFD_Base.cpp
g++ -c -fopenmp -mavx -O3 -ggdb RTFD_Opt.cpp
g++ -c -fopenmp -mavx -O3 -ggdb RTFD_Parallel.cpp
g++ -c -fopenmp -mavx -O3 -ggdb RTFD_SIMD.cpp
g++ -static -Wl,--stack,8388608 -fopenmp -mavx -O3 -ggdb demo_core.o demo_render.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o -o demo_render_windows.exe -lsfml-graphics -lsfml-window -lsfml-system
g++ -static -Wl,--stack,8388608 -fopenmp -mavx -O3 -ggdb demo_core.o demo_perf.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o -o demo_perf_windows.exe
REM rm –path ./*.o
ECHO Build complete
