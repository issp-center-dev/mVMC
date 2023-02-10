# include guard
cmake_minimum_required(VERSION 3.0)

if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

if("${ARCHITECTURE}" STREQUAL "x86_64")
  add_definitions(-DF77_COMPLEX_RET_INTEL)
endif()
include_directories(../pfapack/c_interface)

# Require Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(ltl2inv STATIC ilaenv_lauum.cc gemmt2blas.cc ltl2inv.cc)
target_compile_definitions(ltl2inv PRIVATE -D_CC_IMPL)

