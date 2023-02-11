cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(mVMC NONE)

option(USE_SCALAPACK "Use Scalapack" OFF)
option(PFAFFIAN_BLOCKED "Use blocked-update Pfaffian to speed up." OFF)
option(USE_GEMMT "Use GEMMT. Recommended regardless blocked-Pfaffian-update." ON)

add_definitions(-D_mVMC)
if(CONFIG)
  message(STATUS "Loading configration: " ${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
  include(${PROJECT_SOURCE_DIR}/config/${CONFIG}.cmake)
endif(CONFIG)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

option(GIT_SUBMODULE_UPDATE "execute `git submodule update` automatically" ON)

# First, enables C language only.
# External packages only use their C API.
enable_language(C)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

# TODO: Is this really needed?
if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran} ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  # TODO: Really needs separation?
  if("${CMAKE_C_COMPILER_VERSION}" VERSION_LESS "15.0.0.20140528")
     set(OMP_FLAG_Intel "-openmp")
  else()
    set(OMP_FLAG_Intel "-qopenmp")
  endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OMP_FLAG_Intel}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OMP_FLAG_Intel}")
else()
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif(OPENMP_FOUND)	
endif()

find_package(MPI)
if(MPI_C_FOUND)
  include_directories(${MPI_C_INCLUDE_PATH})
  add_definitions(-D_mpi_use)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
endif(MPI_C_FOUND)

# BLAS/LAPACK/Pfapack77 needs Fortran interface.
enable_language(Fortran)

find_package(LAPACK)
if(USE_SCALAPACK MATCHES OFF)
  if(LAPACK_FOUND)
    add_definitions(-D_lapack)
  endif(LAPACK_FOUND)
endif()

if(PFAFFIAN_BLOCKED
    OR USE_GEMMT)
  set(Require_BLIS ON)
endif()

if(BLA_VENDOR MATCHES "Intel10*"
    OR BLAS_LIBRARIES MATCHES "/*mkl*")
  # Don't require BLIS when MKL is used.
  add_definitions(-DMKL)
  add_definitions(-DBLAS_EXTERNAL)
  add_definitions(-DF77_COMPLEX_RET_INTEL)
  set(Require_BLIS OFF)
elseif(BLA_VENDOR MATCHES "FLA*"
    OR BLAS_LIBRARIES MATCHES "/*blis*")
  # Skip extra BLIS if it's already the BLAS vendor.
  list(GET BLAS_LIBRARIES 0 BLIS_FIRST_LIB)
  get_filename_component(BLIS_LIB_DIR ${BLIS_FIRST_LIB} DIRECTORY)
  include_directories(${BLIS_LIB_DIR}/../include)
  include_directories(${BLIS_LIB_DIR}/../include/blis)
  set(Require_BLIS OFF)
else()
  # BLAS vendor preference:
  # External > BLIS > Reference
  if(DEFINED BLA_VENDOR)
    add_definitions(-DBLAS_EXTERNAL)
  endif()
endif()

# Build and enable tests
# testing setup
# enable_testing() must be called in the top-level CMakeLists.txt before any add_subdirectory() is called.
option(Testing "Enable testing" ON)
if (Testing)
  enable_testing(test)
  add_subdirectory(test)
endif()


# git submodule update
# ref: https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
set(MVMC ON) # Option for upstream StdFace.
find_package(Git QUIET)
set(STDFACE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/StdFace")
if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  # Update submodules as needed
  if(GIT_SUBMODULE_UPDATE)
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()

  if(NOT EXISTS "${STDFACE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR "The submodule StdFace were not downloaded! GIT_SUBMODULE_UPDATE was turned off or failed. Please update submodules and try again.")
  endif()
endif()
add_subdirectory("${STDFACE_DIR}")

# C++ support for pfupdates & ltl2inv.
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

# Require Eigen
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

if(Require_BLIS)
  add_definitions(-D_BLIS)
  include("download_blis_artifact.cmake")
else(Require_BLIS)
  # Use bundled blis.h
  include_directories(src/common/deps)
endif(Require_BLIS)

if(PFAFFIAN_BLOCKED)
  add_definitions(-D_pf_block_update)
endif(PFAFFIAN_BLOCKED)
add_subdirectory(src/pfupdates)
if(Require_BLIS)
  add_dependencies(pfupdates blis_include)
endif(Require_BLIS)

if (Document)
  add_subdirectory(doc)
endif(Document)

add_subdirectory(src/ComplexUHF)
add_subdirectory(src/pfapack/fortran)
add_subdirectory(src/ltl2inv)
  if(Require_BLIS)
    add_dependencies(ltl2inv blis_include)
  endif(Require_BLIS)
add_subdirectory(src/mVMC)
add_subdirectory(tool)
