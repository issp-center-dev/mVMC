# for Fujitsu Compiler in Clang-mode
set(CMAKE_C_COMPILER "mpifccpx" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "mpifCCpx" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-Nclang -Ofast -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-Nclang -Ofast -DNDEBUG" CACHE STRING "" FORCE)

set(CMAKE_Fortran_COMPILER "mpifrtpx" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-Kfast,parallel -DNDEBUG -DFUJITSU" CACHE STRING "" FORCE)

set(OpenMP_C_FLAGS "-fopenmp" CACHE STRING "" FORCE)
set(OpenMP_Fortran_FLAGS "-Kopenmp" CACHE STRING "" FORCE)

# for SSL2
set(BLAS_LIBRARIES "-SSL2" CACHE STRING "" FORCE)
set(LAPACK_LIBRARIES ${BLAS_LIBRARIES} CACHE STRING "" FORCE)

# for BLIS & block pfaffian
set(BLIS_ARTIFACT_CONFIG "a64fx")
set(PFAFFIAN_BLOCKED ON)
