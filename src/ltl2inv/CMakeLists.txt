# include guard
cmake_minimum_required(VERSION 2.8.0 )

if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(../common)

# TODO: Move blalink_gemmt.c to other subprojects?
add_library(ltl2inv STATIC ltl2inv.cc blalink_gemmt.c)
target_compile_definitions(ltl2inv PRIVATE -D_CC_IMPL)

