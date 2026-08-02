add_definitions(-D_mVMC)
# include guard
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

include_directories(include)
include_directories(../common)
include_directories("${STDFACE_DIR}/src")

add_definitions(-DMEXP=19937)

set(SOURCES_vmcmain
        vmcmain.c absolute_pfaffian.c classic_pfaffian_matrix.c
        classic_pfaffian_state.c
        classic_pfaffian_collective.c
        physcal_lanczos.c physcal_lanczos2.c
        lanczos2_contract.c splitloop.c
 )

set(SOURCES_sfmt
        ../sfmt/SFMT.c   
 )

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../pfapack)

add_executable(vmcdry.out vmcdry.c)
target_link_libraries(vmcdry.out StdFace_mvmc m)
add_executable(vmc.out ${SOURCES_vmcmain} ${SOURCES_sfmt})
if(Testing)
  # Keep fault injection, reference samplers, and audits out of production.
  target_sources(vmc.out PRIVATE classic_pfaffian_sampler.c
      classic_pfaffian_audit.c krylov_fock_reference.c
      classic_krylov_amplitude.c classic_krylov_model.c)
  target_compile_definitions(vmc.out PRIVATE MVMC_ENABLE_FAULT_INJECTION
      MVMC_ENABLE_PFAFFIAN_SAMPLER_REFERENCE
      MVMC_ENABLE_PFAFFIAN_STATE_AUDIT
      MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
      MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
endif(Testing)
target_link_libraries(vmc.out StdFace_mvmc m)
target_link_libraries(vmc.out pfapack)
target_link_libraries(vmc.out ltl2inv)
if(PFAFFIAN_BLOCKED)
  target_link_libraries(vmc.out pfupdates)
endif(PFAFFIAN_BLOCKED)
if(Require_BLIS)
  target_link_libraries(vmc.out blis pthread)
endif(Require_BLIS)
target_link_libraries(vmc.out ${LAPACK_LIBRARIES} m)

if(USE_SCALAPACK)
  string(REGEX REPLACE "-L[ ]+" "-L" sc_libs "${SCALAPACK_LIBRARIES}")
  string(REGEX REPLACE "[ ]+" ";" sc_libs "${sc_libs}")
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
