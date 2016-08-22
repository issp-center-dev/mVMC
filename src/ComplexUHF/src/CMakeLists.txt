# include guard
cmake_minimum_required(VERSION 2.8.0 )

if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of UHF.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(include)

set(SOURCES_UHF
        UHFmain.c output.c cal_energy.c cal_cisajs.c green.c makeham.c diag.c initial.c check.c matrixlapack.c readdef.c
 )


add_executable(UHF ${SOURCES_UHF} )
target_link_libraries(UHF ${LAPACK_LIBRARIES})
if(MPI_FOUND)
  target_link_libraries(UHF ${MPI_C_LIBRARIES})
endif(MPI_FOUND)

