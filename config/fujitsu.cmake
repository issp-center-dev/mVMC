# for Fujitsu Compiler
set(CMAKE_C_COMPILER "mpifccpx" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "mpiFCCpx" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG -Kfast,parallel -Kmemalias,alias_const" CACHE STRING "" FORCE)
set(OpenMP_C_FLAGS "-Kopenmp" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "mpifrtpx" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-DNDEBUG  -DFUJITSU -Kfast,parallel" CACHE STRING "" FORCE)
set(OpenMP_Fortran_FLAGS "-Kopenmp" CACHE STRING "" FORCE)
set(MPI_EXTRA_LIBRARY "-lmpi_f90 -lmpi_f77" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-Kopenmp" CACHE STRING "" FORCE)

# for SSL2
set(BLAS_LIBRARIES "-SSL2" CACHE STRING "" FORCE)
set(LAPACK_LIBRARIES ${BLAS_LIBRARIES} CACHE STRING "" FORCE)
