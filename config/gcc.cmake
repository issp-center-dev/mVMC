# for GCC Compiler
set(CMAKE_C_COMPILER "gcc" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall  -Wformat -Werror=format-security")
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wno-unknown-pragmas ")
Set(CMAKE_Fortran_COMPILER "gfortran" CACHE STRING "" FORCE)
