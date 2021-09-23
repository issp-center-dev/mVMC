# include guard
cmake_minimum_required(VERSION 2.8.0 )

if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_definitions(-DBLAS_EXTERNAL)
add_definitions(-DF77_COMPLEX_RET_INTEL)
include_directories(../pfaffine/src)

add_library(pfupdates STATIC
  pf_interface.cc
  ../pfaffine/src/skpfa.cc
  ../pfaffine/src/sktdf.cc
  ../pfaffine/src/sktdi.cc)
target_compile_definitions(pfupdates PRIVATE -D_CC_IMPL)

