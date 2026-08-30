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
        classic_krylov_projection.c classic_krylov_amplitude.c
        classic_krylov_model.c power_lanczos_classic_bridge.c
        classic_pfaffian_collective.c
        physcal_lanczos.c physcal_lanczos2.c
        power_lanczos_contract.c power_lanczos_fp_metadata.c
        power_lanczos_observable_census.c
        power_lanczos_observable_registry.c
        power_lanczos_observable_action.c
        power_lanczos_observable_evaluator.c
        power_lanczos_observable_blocks.c
        power_lanczos_observable_preflight.c
        power_lanczos_observable_confirmation.c
        power_lanczos_sha256.c
        power_lanczos_rng.c
        power_lanczos_chain_rng.c
        power_lanczos_snapshot.c
        power_lanczos_production_session.c
        power_lanczos_production_session_prepare.c
        power_lanczos_production_session_stages.c
        power_lanczos_primitive_trace.c
        power_lanczos_numeric_evidence.c
        power_lanczos_stabilization_statistics.c
        power_lanczos_stabilization_output.c
        power_lanczos_corrected_dispatch.c
        power_lanczos_json_writer.c
        power_lanczos_output_transaction.c
        power_lanczos_block_collective.c
        power_lanczos_selector.c
        lanczos2_contract.c splitloop.c
        krylov_fock_reference.c krylov_fock_proposal.c bounded_krylov_engine.c
        bounded_krylov_collective.c krylov_positive_guide.c
        krylov_positive_sampler.c krylov_streaming_statistics.c
        krylov_matrix_measurement.c krylov_gevp_solver.c
        power_lanczos_gevp.c
        krylov_final_state.c krylov_final_state_chain.c
 )

set(SOURCES_sfmt
        ../sfmt/SFMT.c   
 )

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../pfapack)

add_executable(vmcdry.out vmcdry.c)
target_link_libraries(vmcdry.out StdFace_mvmc m)

function(configure_vmc_executable target include_testing_support)
  add_executable(${target} ${SOURCES_vmcmain} ${SOURCES_sfmt})
  target_compile_definitions(${target} PRIVATE
      MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
      MVMC_ENABLE_POWER_LANCZOS_P5_CORE)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE
        -fno-fast-math -ffp-contract=off)
    target_compile_definitions(${target} PRIVATE
        MVMC_POWER_LANCZOS_FP_CONTRACT_OFF_REQUESTED=1)
  endif()
  if(include_testing_support)
    # Keep fault injection, reference samplers, and audits out of production.
    target_sources(${target} PRIVATE classic_pfaffian_sampler.c
        classic_pfaffian_audit.c)
    target_compile_definitions(${target} PRIVATE MVMC_ENABLE_FAULT_INJECTION
        MVMC_ENABLE_PFAFFIAN_SAMPLER_REFERENCE
        MVMC_ENABLE_PFAFFIAN_STATE_AUDIT
        MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
  endif()
  target_link_libraries(${target} StdFace_mvmc m pfapack ltl2inv)
  if(PFAFFIAN_BLOCKED)
    target_link_libraries(${target} pfupdates)
  endif()
  if(Require_BLIS)
    target_link_libraries(${target} blis pthread)
  endif()
  target_link_libraries(${target} ${LAPACK_LIBRARIES} m)

  if(USE_SCALAPACK)
    string(REGEX REPLACE "-L[ ]+" "-L" sc_libs
        "${SCALAPACK_LIBRARIES}")
    string(REGEX REPLACE "[ ]+" ";" sc_libs "${sc_libs}")
    foreach(sc_lib IN LISTS sc_libs)
      target_link_libraries(${target} ${sc_lib})
    endforeach(sc_lib)
  endif()

  if(MPI_FOUND)
    target_link_libraries(${target} ${MPI_C_LIBRARIES})
  endif()
endfunction()

configure_vmc_executable(vmc.out ${Testing})
if(Testing)
  # CTest builds this target on demand and applies every cumulative production
  # symbol policy to a binary compiled without Testing-only hooks.
  configure_vmc_executable(vmc.production-policy.out FALSE)
  set_target_properties(vmc.production-policy.out PROPERTIES
      EXCLUDE_FROM_ALL TRUE)
endif()

if(USE_SCALAPACK)
  message(STATUS "SCALAPACK_LIBRARIES: ${SCALAPACK_LIBRARIES}")
endif()
install(TARGETS vmcdry.out RUNTIME DESTINATION bin)
install(TARGETS vmc.out RUNTIME DESTINATION bin)
add_definitions(-D_mVMC)
