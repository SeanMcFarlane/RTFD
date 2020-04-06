ECHO BEGIN COMPILE
g++ -static -static-libgcc -static-libstdc++ -c -fopenmp -mavx -O3 -ggdb demo_core.cpp
g++ -static -static-libgcc -static-libstdc++ -c -fopenmp -mavx -O3 -ggdb demo_perf.cpp
g++ -static -static-libgcc -static-libstdc++ -c -fopenmp -mavx -O3 -ggdb demo_render.cpp
g++ -static -static-libgcc -static-libstdc++ -c -fopenmp -mavx -O3 -ggdb RTFD_Base.cpp
g++ -static -static-libgcc -static-libstdc++ -c -fopenmp -mavx -O3 -ggdb RTFD_Opt.cpp
g++ -static -static-libgcc -static-libstdc++ -c -fopenmp -mavx -O3 -ggdb RTFD_Parallel.cpp
g++ -static -static-libgcc -static-libstdc++ -c -fopenmp -mavx -O3 -ggdb RTFD_SIMD.cpp
g++ -static -static-libgcc -static-libstdc++ -Wl,--stack,8388608 -fopenmp -mavx -O3 -ggdb demo_core.o demo_render.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o -o demo_render_windows.exe -lsfml-graphics -lsfml-window -lsfml-system
g++ -static -static-libgcc -static-libstdc++ -Wl,--stack,8388608 -fopenmp -mavx -O3 -ggdb demo_core.o demo_perf.o RTFD_Base.o RTFD_Opt.o RTFD_Parallel.o RTFD_SIMD.o -o demo_perf_windows.exe
del "*.o" /f /q
ECHO Build complete
