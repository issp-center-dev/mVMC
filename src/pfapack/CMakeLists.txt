# include guard
cmake_minimum_required(VERSION 2.8)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

set(SOURCES_PFAPACK dlasktrd.f dskmv.f dskr2.f dsktd2.f dsktrd.f zlasktrd.f zskmv.f zskr2.f zsktd2.f zsktrd.f dlasktrf.f dskpfa.f dskr2k.f dsktf2.f dsktrf.f zlasktrf.f zskpfa.f zskr2k.f zsktf2.f zsktrf.f)
include_directories(${MPI_Fortran_INCLUDE_PATH})
add_library(pfapack ${SOURCES_PFAPACK})
target_link_libraries(pfapack ${MPI_Fortran_LIBRARIES} ${LAPACK_LIBRARIES})
install(TARGETS pfapack
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin)
