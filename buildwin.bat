
ECHO BEGIN COMPILE
ECHO OFF
CL /c /openmp /arch:AVX /O2  /DEBUG:FULL /Zi  demo_core.cpp
CL /c /EHsc /openmp /arch:AVX /O2  /DEBUG:FULL /Zi  /IWin64\include demo_perf.cpp
CL /c /EHsc /openmp /arch:AVX /O2  /DEBUG:FULL /Zi  /IWin64\include demo_render.cpp
CL /c /O2  /DEBUG:FULL /Zi  RTFD_Base.cpp
CL /c /O2  /DEBUG:FULL /Zi  RTFD_Opt.cpp
CL /c /openmp /O2  /DEBUG:FULL /Zi  RTFD_Parallel.cpp
CL /c /EHsc /openmp /arch:AVX /O2  /DEBUG:FULL /Zi  RTFD_SIMD.cpp
nvcc -c -g -G RTFD_CUDA.cu -arch=sm_50 ^
 -gencode=arch=compute_50,code=sm_50 ^
 -gencode=arch=compute_52,code=sm_52 ^
 -gencode=arch=compute_60,code=sm_60 ^
 -gencode=arch=compute_61,code=sm_61 ^
 -gencode=arch=compute_70,code=sm_70 ^
 -gencode=arch=compute_70,code=compute_70

CL /openmp /arch:AVX /O2  /DEBUG:FULL /Zi  /F 8388608 demo_core.obj demo_render.obj RTFD_Base.obj RTFD_Opt.obj RTFD_Parallel.obj RTFD_SIMD.obj RTFD_CUDA.obj /Fe:Win64/demo_render_windows.exe^
 "Win64\lib\sfml-graphics.lib" "Win64\lib\sfml-system.lib" "Win64\lib\sfml-window.lib" "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.2\lib\Win32\*.lib" "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.2\lib\x64\*.lib" 

CL /openmp /arch:AVX /O2  /DEBUG:FULL /Zi  /F 8388608 demo_core.obj demo_perf.obj RTFD_Base.obj RTFD_Opt.obj RTFD_Parallel.obj RTFD_SIMD.obj RTFD_CUDA.obj /Fe:Win64/demo_perf_windows.exe^
 "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.2\lib\Win32\*.lib" "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.2\lib\x64\*.lib" 
 
del "*.obj" /f /q
ECHO Build complete
