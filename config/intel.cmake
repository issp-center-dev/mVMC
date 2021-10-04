# for Intel Compiler

set(CMAKE_C_COMPILER "icc" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "icpc" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wformat -Werror=format-security")
set(CMAKE_C_FLAGS_RELEASE "-Wno-unknown-pragmas -O3 -DNDEBUG -DHAVE_SSE2" CACHE STRING "" FORCE)
#set(CMAKE_CXX_FLAGS="-std=c++11" FORCE)

set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -DHAVE_SSE2" CACHE STRING "" FORCE)

# for Intel MKL
set(BLA_VENDOR "Intel10_64lp" CACHE STRING "" FORCE)

# for BLIS & block pfaffian
set(BLIS_ARTIFACT_CONFIG "intel64")
set(PFAFFIAN_BLOCKED OFF)
