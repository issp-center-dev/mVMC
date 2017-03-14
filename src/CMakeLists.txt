# include guard
cmake_minimum_required(VERSION 2.8.0 )

add_definitions(-D_mVMC)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_subdirectory(pfapack)

include_directories(include)
include_directories(sfmt)
add_definitions(-DMEXP=19937)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/pfapack)

set(SOURCES_vmcmain
        vmcmain.c   
 )

set(SOURCES_sfmt
        sfmt/SFMT.c   
 )

add_library(StdFace STATIC StdFace/ChainLattice.c StdFace/HoneycombLattice.c StdFace/SquareLattice.c StdFace/StdFace_main.c StdFace/StdFace_ModelUtil.c StdFace/TriangularLattice.c StdFace/Ladder.c StdFace/Kagome.c StdFace/Orthorhombic.c StdFace/Pyrochlore.c StdFace/Wannier90.c StdFace/FCOrtho.c)


add_executable(vmcdry.out vmcdry.c)
target_link_libraries(vmcdry.out StdFace)
add_executable(vmc.out ${SOURCES_vmcmain} ${SOURCES_sfmt})

target_link_libraries(vmc.out pfapack StdFace ${LAPACK_LIBRARIES})
if(MPI_FOUND)
  target_link_libraries(vmc.out ${MPI_C_LIBRARIES})
endif(MPI_FOUND)
install(TARGETS vmc.out RUNTIME DESTINATION bin)

