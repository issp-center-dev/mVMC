# for AOCC

set(CMAKE_C_COMPILER "clang" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "clang++" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall -Wformat -Werror=format-security")
set(CMAKE_C_FLAGS_RELEASE "-Wno-unknown-pragmas -O3 -DNDEBUG -DHAVE_SSE2" CACHE STRING "" FORCE)
#set(CMAKE_CXX_FLAGS="-std=c++11")

set(CMAKE_Fortran_COMPILER "flang" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -DHAVE_SSE2" CACHE STRING "" FORCE)

# OpenMP, libatomic
set(CMAKE_EXE_LINKER_FLAGS "-fopenmp -latomic") 

# for AOCL
set(BLA_VENDOR "FLAME" CACHE STRING "" FORCE)

# for BLIS & block pfaffian
set(BLIS_ARTIFACT_CONFIG "amd64")
#set(PFAFFIAN_BLOCKED OFF)
