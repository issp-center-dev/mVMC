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
        classic_pfaffian_state.c classic_krylov_configuration.c
        classic_pfaffian_collective.c
        physcal_lanczos.c physcal_lanczos2.c
        power_lanczos_contract.c power_lanczos_fp_metadata.c
        lanczos2_contract.c splitloop.c
        krylov_fock_reference.c krylov_fock_proposal.c bounded_krylov_engine.c
        bounded_krylov_collective.c krylov_positive_guide.c
        krylov_positive_sampler.c krylov_streaming_statistics.c
        krylov_matrix_measurement.c krylov_gevp_solver.c
        krylov_final_state.c krylov_final_state_chain.c
 )

set(SOURCES_sfmt
        ../sfmt/SFMT.c   
 )

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../pfapack)

add_executable(vmcdry.out vmcdry.c)
target_link_libraries(vmcdry.out StdFace_mvmc m)
add_executable(vmc.out ${SOURCES_vmcmain} ${SOURCES_sfmt})
target_compile_definitions(vmc.out PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
    MVMC_ENABLE_POWER_LANCZOS_P5_CORE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(vmc.out PRIVATE -fno-fast-math -ffp-contract=off)
  target_compile_definitions(vmc.out PRIVATE
      MVMC_POWER_LANCZOS_FP_CONTRACT_OFF_REQUESTED=1)
endif()
if(Testing)
  # Keep fault injection, reference samplers, and audits out of production.
  target_sources(vmc.out PRIVATE classic_pfaffian_sampler.c
      classic_pfaffian_audit.c
      classic_krylov_amplitude.c classic_krylov_model.c)
  target_compile_definitions(vmc.out PRIVATE MVMC_ENABLE_FAULT_INJECTION
      MVMC_ENABLE_PFAFFIAN_SAMPLER_REFERENCE
      MVMC_ENABLE_PFAFFIAN_STATE_AUDIT
      MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
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
