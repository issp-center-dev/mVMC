# include guard
cmake_minimum_required(VERSION 2.8.0 )

if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_definitions(-DBLAS_EXTERNAL)
add_definitions(-DF77_COMPLEX_RET_INTEL)
include_directories(../common)

add_library(ltl2inv STATIC ltl2inv.cc)
target_compile_definitions(ltl2inv PRIVATE -D_CC_IMPL)

