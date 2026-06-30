# include guard
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

# F77_COMPLEX_RET_INTEL is a BLAS vendor ABI choice, not a CPU architecture
# choice. The top-level CMakeLists.txt enables it for MKL/Intel BLAS only.
include_directories(../common)

# blalink_gemmt.c requires BLIS (bli_*gemmt); only needed when USE_GEMMT is ON.
set(SOURCES_ltl2inv ltl2inv.cc ilaenv_lauum.cc ilaenv_wrap.f90)
if(USE_GEMMT)
  list(APPEND SOURCES_ltl2inv blalink_gemmt.c)
endif()
add_library(ltl2inv STATIC ${SOURCES_ltl2inv})
target_compile_definitions(ltl2inv PRIVATE -D_CC_IMPL)
