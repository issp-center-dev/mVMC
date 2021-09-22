# for Apple clang compiler
# additional libomp and gfortran installation required
# mac computers are suggested to use this configuration for better performance
set(CMAKE_C_COMPILER "clang" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "clang++" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall  -Wformat -Werror=format-security")
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wno-unknown-pragmas -Wno-logical-not-parentheses")
set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "" FORCE)

# OpenMP with libomp
set(CMAKE_EXE_LINKER_FLAGS "-L/usr/local/lib -lomp") 
set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "" FORCE) 
set(OpenMP_C_LIB_NAMES "")
set(CMAKE_OSX_SYSROOT "")

