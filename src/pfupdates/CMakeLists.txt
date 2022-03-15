# include guard
cmake_minimum_required(VERSION 2.8.0 )

if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(../common)
include_directories(../ltl2inv)

add_library(pfupdates STATIC pf_interface.cc)
target_compile_definitions(pfupdates PRIVATE -D_CC_IMPL)

