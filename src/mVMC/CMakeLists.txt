# include guard
cmake_minimum_required(VERSION 2.8.0 )

add_definitions(-D_mVMC)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(include)
include_directories(../common)
include_directories("${STDFACE_DIR}/src")

add_definitions(-DMEXP=19937)

set(SOURCES_vmcmain
        vmcmain.c physcal_lanczos.c splitloop.c 
 )

set(SOURCES_sfmt
        ../sfmt/SFMT.c   
 )

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../pfapack)

add_executable(vmcdry.out vmcdry.c)
target_link_libraries(vmcdry.out StdFace m)
add_executable(vmc.out ${SOURCES_vmcmain} ${SOURCES_sfmt})
target_link_libraries(vmc.out StdFace)
target_link_libraries(vmc.out pfapack)
if(PFAFFIAN_BLOCKED)
  target_link_libraries(vmc.out pfupdates blis pthread)
endif(PFAFFIAN_BLOCKED)
target_link_libraries(vmc.out ${LAPACK_LIBRARIES} m)

if(USE_SCALAPACK)
  string(REGEX REPLACE "-L[ ]+" "-L" sc_libs ${SCALAPACK_LIBRARIES})
  string(REGEX REPLACE "[ ]+" ";" sc_libs ${sc_libs})
  foreach(sc_lib IN LISTS sc_libs)
    target_link_libraries(vmc.out ${sc_lib})
  endforeach(sc_lib)
  message(STATUS "SCALAPACK_LIBRARIES: ${SCALAPACK_LIBRARIES}")
endif(USE_SCALAPACK)

if(MPI_FOUND)
  target_link_libraries(vmc.out ${MPI_C_LIBRARIES})
endif(MPI_FOUND)
install(TARGETS vmcdry.out RUNTIME DESTINATION bin)
install(TARGETS vmc.out RUNTIME DESTINATION bin)
add_definitions(-D_mVMC)
