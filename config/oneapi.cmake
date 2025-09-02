# for Intel Compiler

set(CMAKE_C_COMPILER "icx" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "icpx" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall -Wformat -Werror=format-security")
set(CMAKE_C_FLAGS_RELEASE "-Wno-unknown-pragmas -O2 -DNDEBUG -xHost" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE " -Wno-unknown-pragmas -O2 -DNDEBUG -xHost" CACHE STRING "" FORCE)

set(CMAKE_Fortran_COMPILER "ifx" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -xHost" CACHE STRING "" FORCE)

# for Intel MKL
set(BLA_VENDOR "Intel10_64lp" CACHE STRING "" FORCE)

# for BLIS & block pfaffian
set(BLIS_ARTIFACT_CONFIG "intel64")

if(USE_SCALAPACK)
  if(SCALAPACK_LIBRARIES MATCHES "")
     set(SCALAPACK_LIBRARIES "-qmkl=cluster")
     #set(SCALAPACK_LIBRARIES "-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm")
  endif(SCALAPACK_LIBRARIES MATCHES "")

  message(STATUS "SCALAPACK_LIBRARIES is ${SCALAPACK_LIBRARIES}")
endif(USE_SCALAPACK)

