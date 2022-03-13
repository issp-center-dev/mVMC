cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(mVMC NONE)

option(USE_SCALAPACK "Use Scalapack" OFF)
option(PFAFFIAN_BLOCKED "Use blocked-update Pfaffian to speed up." OFF)

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
enable_language(C CXX Fortran)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_MACOSX_RPATH 1)

if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran} ${CMAKE_Fortran_FLAGS_RELEASE}")
endif()

if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
  if("${CMAKE_C_COMPILER_VERSION}" VERSION_LESS "15.0.0.20140528")
     set(OMP_FLAG_Intel "-openmp")
  else()
    set(OMP_FLAG_Intel "-qopenmp")
  endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OMP_FLAG_Intel}")
else()
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  endif(OPENMP_FOUND)	
endif()

find_package(MPI)
if(MPI_C_FOUND)
  include_directories(${MPI_C_INCLUDE_PATH})
  add_definitions(-D_mpi_use)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
endif(MPI_C_FOUND)

find_package(LAPACK)
if(USE_SCALAPACK MATCHES OFF)
  if(LAPACK_FOUND)
    add_definitions(-D_lapack)
  endif(LAPACK_FOUND)
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

if(PFAFFIAN_BLOCKED)
  add_definitions(-D_pf_block_update)
  include("download_blis_artifact.cmake")
  # Must set BLIS artifact BEFORE adding pfupdates target.
  add_subdirectory(src/pfupdates)
  add_dependencies(pfupdates blis_include)
else(PFAFFIAN_BLOCKED)
  # Use bundled blis.h
  include_directories(src/common/deps)
endif(PFAFFIAN_BLOCKED)
add_subdirectory(src/ltl2inv)

if (Document)
  add_subdirectory(doc)
endif(Document)

add_subdirectory(src/ComplexUHF)
add_subdirectory(src/pfapack/fortran)
add_subdirectory(src/mVMC)
add_subdirectory(tool)
