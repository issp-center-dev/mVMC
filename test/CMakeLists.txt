if(NOT PYTHON_EXECUTABLE)
    find_package(PythonInterp 3.6 REQUIRED)
endif()

foreach(datafile ${DATA_FILES})
    configure_file(${datafile} ${CMAKE_BINARY_DIR}/test/ COPYONLY)
endforeach()

if(LAPACK_FOUND)
    add_executable(blas_dgemm_canary blas_dgemm_canary.c)
    target_link_libraries(blas_dgemm_canary ${LAPACK_LIBRARIES} m)
    add_test(NAME BLAS_DGEMM_Canary COMMAND blas_dgemm_canary)
endif()

add_executable(nbody_operator_unit
    nbody_operator_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/nbody_operator.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/backflow_nbody_scratch.c)
target_include_directories(nbody_operator_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
add_test(NAME NBodyOperator_Unit COMMAND nbody_operator_unit)

add_executable(lanczos2_contract_unit
    lanczos2_contract_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/lanczos2_contract.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_contract.c)
target_include_directories(lanczos2_contract_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
add_test(NAME Lanczos2Contract_Unit COMMAND lanczos2_contract_unit)

add_executable(power_lanczos_contract_unit
    power_lanczos_contract_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_contract.c)
target_include_directories(power_lanczos_contract_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_contract_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
add_test(NAME PowerLanczosContract_Unit
    COMMAND power_lanczos_contract_unit)

add_executable(power_lanczos_production_symbol_probe
    power_lanczos_production_symbol_probe.c)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_production_symbol_probe PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
add_test(NAME PowerLanczosProductionSymbols_OracleForbidden
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/runtest_power_lanczos_production_symbols.py
            --python ${PYTHON_EXECUTABLE}
            --verifier
            ${CMAKE_CURRENT_SOURCE_DIR}/python/verify_power_lanczos_p6b_production_symbols.py
            --binary $<TARGET_FILE:power_lanczos_production_symbol_probe>
            --policy
            ${CMAKE_CURRENT_SOURCE_DIR}/python/data/power_lanczos_symbol_policy_oracle_probe.json)

add_executable(power_lanczos_observable_census_unit
    power_lanczos_observable_census_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_census.c)
target_include_directories(power_lanczos_observable_census_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_observable_census_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
add_test(NAME PowerLanczosObservableCensus_Unit
    COMMAND power_lanczos_observable_census_unit)

add_executable(power_lanczos_observable_registry_unit
    power_lanczos_observable_registry_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_census.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_registry.c)
target_include_directories(power_lanczos_observable_registry_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_observable_registry_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
add_test(NAME PowerLanczosObservableRegistry_Unit
    COMMAND power_lanczos_observable_registry_unit)

add_executable(power_lanczos_observable_action_unit
    power_lanczos_observable_action_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_action.c)
target_include_directories(power_lanczos_observable_action_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_observable_action_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
add_test(NAME PowerLanczosObservableAction_Unit
    COMMAND power_lanczos_observable_action_unit)

add_executable(power_lanczos_observable_evaluator_unit
    power_lanczos_observable_evaluator_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_evaluator.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_action.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_include_directories(power_lanczos_observable_evaluator_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(power_lanczos_observable_evaluator_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_observable_evaluator_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_observable_evaluator_unit
    pfapack ltl2inv ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(power_lanczos_observable_evaluator_unit
        blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(power_lanczos_observable_evaluator_unit
        ${MPI_C_LIBRARIES})
endif()
add_test(NAME PowerLanczosObservableEvaluator_Unit
    COMMAND power_lanczos_observable_evaluator_unit)

add_executable(power_lanczos_observable_blocks_unit
    power_lanczos_observable_blocks_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_blocks.c)
target_include_directories(power_lanczos_observable_blocks_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(power_lanczos_observable_blocks_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_observable_blocks_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_observable_blocks_unit m)
add_test(NAME PowerLanczosObservableBlocks_Unit
    COMMAND power_lanczos_observable_blocks_unit)

add_executable(power_lanczos_observable_preflight_unit
    power_lanczos_observable_preflight_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_preflight.c)
target_include_directories(power_lanczos_observable_preflight_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_observable_preflight_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_observable_preflight_unit m)
add_test(NAME PowerLanczosObservablePreflight_Unit
    COMMAND power_lanczos_observable_preflight_unit)

add_executable(power_lanczos_production_session_unit
    power_lanczos_production_session_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_sha256.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_rng.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_snapshot.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_production_session.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_production_session_prepare.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_production_session_stages.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_block_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_chain_rng.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_census.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_action.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_evaluator.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_blocks.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_confirmation.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_primitive_trace.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_stabilization_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_stabilization_output.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_json_writer.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_matrix_measurement.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_streaming_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_gevp.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_gevp_solver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state_chain.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_include_directories(power_lanczos_production_session_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(power_lanczos_production_session_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
    MVMC_ENABLE_POWER_LANCZOS_P5_CORE
    MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_production_session_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_production_session_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(power_lanczos_production_session_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(power_lanczos_production_session_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(power_lanczos_production_session_unit
        ${MPI_C_LIBRARIES})
endif()
add_test(NAME PowerLanczosProductionSession_Unit
    COMMAND power_lanczos_production_session_unit)
if(MPI_FOUND)
    add_test(NAME PowerLanczosProductionSession_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_production_session_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME PowerLanczosProductionSession_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_production_session_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
        PowerLanczosProductionSession_Unit_Rank2_mpi PROPERTIES PROCESSORS 2)
    set_tests_properties(
        PowerLanczosProductionSession_Unit_Rank4_mpi PROPERTIES PROCESSORS 4)
endif()

add_executable(power_lanczos_chain_rng_unit
    power_lanczos_chain_rng_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_rng.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_chain_rng.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_include_directories(power_lanczos_chain_rng_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(power_lanczos_chain_rng_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_chain_rng_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_chain_rng_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(power_lanczos_chain_rng_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(power_lanczos_chain_rng_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(power_lanczos_chain_rng_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME PowerLanczosChainRng_Unit
    COMMAND power_lanczos_chain_rng_unit)
if(MPI_FOUND)
    add_test(NAME PowerLanczosChainRng_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_chain_rng_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME PowerLanczosChainRng_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_chain_rng_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
        PowerLanczosChainRng_Unit_Rank2_mpi PROPERTIES PROCESSORS 2)
    set_tests_properties(
        PowerLanczosChainRng_Unit_Rank4_mpi PROPERTIES PROCESSORS 4)
endif()

add_executable(power_lanczos_block_collective_unit
    power_lanczos_block_collective_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_block_collective.c)
target_include_directories(power_lanczos_block_collective_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(power_lanczos_block_collective_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_block_collective_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_block_collective_unit m)
if(MPI_FOUND)
    target_link_libraries(power_lanczos_block_collective_unit
        ${MPI_C_LIBRARIES})
endif()
add_test(NAME PowerLanczosBlockCollective_Unit
    COMMAND power_lanczos_block_collective_unit)
if(MPI_FOUND)
    add_test(NAME PowerLanczosBlockCollective_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_block_collective_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME PowerLanczosBlockCollective_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_block_collective_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
        PowerLanczosBlockCollective_Unit_Rank2_mpi PROPERTIES PROCESSORS 2)
    set_tests_properties(
        PowerLanczosBlockCollective_Unit_Rank4_mpi PROPERTIES PROCESSORS 4)
endif()

add_executable(power_lanczos_primitive_trace_unit
    power_lanczos_primitive_trace_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_primitive_trace.c)
target_include_directories(power_lanczos_primitive_trace_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(power_lanczos_primitive_trace_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_primitive_trace_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(power_lanczos_primitive_trace_unit m)
add_test(NAME PowerLanczosPrimitiveTrace_Unit
    COMMAND power_lanczos_primitive_trace_unit)

add_executable(power_lanczos_numeric_evidence_unit
    power_lanczos_numeric_evidence_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_include_directories(power_lanczos_numeric_evidence_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(power_lanczos_numeric_evidence_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_numeric_evidence_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_numeric_evidence_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(power_lanczos_numeric_evidence_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(power_lanczos_numeric_evidence_unit
        ${MPI_C_LIBRARIES})
endif()
add_test(NAME PowerLanczosNumericEvidence_Unit
    COMMAND power_lanczos_numeric_evidence_unit)

add_executable(power_lanczos_json_writer_unit
    power_lanczos_json_writer_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_json_writer.c)
target_include_directories(power_lanczos_json_writer_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(power_lanczos_json_writer_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_json_writer_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(power_lanczos_json_writer_unit m)
add_test(NAME PowerLanczosJsonWriter_Unit
    COMMAND power_lanczos_json_writer_unit)

add_executable(power_lanczos_output_transaction_unit
    power_lanczos_output_transaction_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_output_transaction.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_sha256.c)
target_include_directories(power_lanczos_output_transaction_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(power_lanczos_output_transaction_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
    MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_output_transaction_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(power_lanczos_output_transaction_unit m)
add_test(NAME PowerLanczosOutputTransaction_Unit
    COMMAND power_lanczos_output_transaction_unit)

add_executable(power_lanczos_runtime_unit
    power_lanczos_selector_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_selector.c)
target_include_directories(power_lanczos_runtime_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_runtime_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
add_test(NAME PowerLanczosRuntime_Unit
    COMMAND power_lanczos_runtime_unit)

add_executable(power_lanczos_observable_correctness
    power_lanczos_observable_correctness.c
    power_lanczos_observable_dense_oracle.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_evaluator.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_action.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_blocks.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_include_directories(power_lanczos_observable_correctness PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(power_lanczos_observable_correctness PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_observable_correctness PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_observable_correctness
    pfapack ltl2inv ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(power_lanczos_observable_correctness blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(power_lanczos_observable_correctness
        ${MPI_C_LIBRARIES})
endif()
add_test(NAME PowerLanczosObservableCorrectness_Unit
    COMMAND power_lanczos_observable_correctness)
add_executable(power_lanczos_observable_scaling_driver
    power_lanczos_observable_scaling_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_census.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_action.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_evaluator.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_blocks.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(power_lanczos_observable_scaling_driver PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(power_lanczos_observable_scaling_driver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include ${CMAKE_CURRENT_SOURCE_DIR})
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_observable_scaling_driver PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_observable_scaling_driver
    pfapack ltl2inv ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(power_lanczos_observable_scaling_driver stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(power_lanczos_observable_scaling_driver
        blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(power_lanczos_observable_scaling_driver
        ${MPI_C_LIBRARIES})
endif()
add_test(NAME PowerLanczosObservableScalingCampaign_Unit
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/power_lanczos_observable_scaling_campaign.py
            self-test)
add_test(NAME PowerLanczosObservableScalingPackage_Lifecycle
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/package_power_lanczos_observable_scaling.py
            self-test-lifecycle)
set(_p6c2_scaling_driver_test_args
    --python ${PYTHON_EXECUTABLE}
    --campaign
        ${CMAKE_CURRENT_SOURCE_DIR}/python/power_lanczos_observable_scaling_campaign.py
    --driver $<TARGET_FILE:power_lanczos_observable_scaling_driver>
    --repo ${CMAKE_SOURCE_DIR}
    --build-dir ${CMAKE_BINARY_DIR})
if(MPI_FOUND)
    list(APPEND _p6c2_scaling_driver_test_args
        --mpiexec ${MPIEXEC_EXECUTABLE})
endif()
add_test(NAME PowerLanczosObservableScalingDriver_Focused
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/runtest_power_lanczos_observable_scaling_driver.py
            ${_p6c2_scaling_driver_test_args})
set_tests_properties(
    PowerLanczosObservableScalingCampaign_Unit
    PowerLanczosObservableScalingPackage_Lifecycle
    PowerLanczosObservableScalingDriver_Focused
    PROPERTIES TIMEOUT 180
    ENVIRONMENT
        "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")

if(MPI_FOUND)
    add_test(NAME PowerLanczosObservableCorrectness_Unit_mpi2
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_observable_correctness>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME PowerLanczosObservableCorrectness_Unit_mpi4
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_observable_correctness>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
        PowerLanczosObservableCorrectness_Unit_mpi2 PROPERTIES PROCESSORS 2)
    set_tests_properties(
        PowerLanczosObservableCorrectness_Unit_mpi4 PROPERTIES PROCESSORS 4)
    add_executable(power_lanczos_observable_census_mpi
        power_lanczos_observable_census_mpi.c
        ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_census.c)
    target_include_directories(power_lanczos_observable_census_mpi PRIVATE
        ${CMAKE_SOURCE_DIR}/src/mVMC/include)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(power_lanczos_observable_census_mpi PRIVATE
            -Wall -Wextra -Wpedantic -Werror)
    endif()
    target_link_libraries(power_lanczos_observable_census_mpi
        ${MPI_C_LIBRARIES})
    add_test(NAME PowerLanczosObservableCensus_Unit_mpi2
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_observable_census_mpi>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME PowerLanczosObservableCensus_Unit_mpi4
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_observable_census_mpi>
                ${MPIEXEC_POSTFLAGS})
endif()

add_executable(power_lanczos_build_metadata_unit
    power_lanczos_fp_metadata_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_fp_metadata.c)
target_include_directories(power_lanczos_build_metadata_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_build_metadata_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
    target_compile_definitions(power_lanczos_build_metadata_unit PRIVATE
        MVMC_POWER_LANCZOS_FP_CONTRACT_OFF_REQUESTED=1)
endif()
add_test(NAME PowerLanczosBuildMetadata_Unit
    COMMAND power_lanczos_build_metadata_unit)

add_executable(classic_krylov_configuration_unit
    classic_krylov_configuration_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c)
target_include_directories(classic_krylov_configuration_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(classic_krylov_configuration_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(classic_krylov_configuration_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
add_test(NAME ClassicKrylovConfiguration_Unit
    COMMAND classic_krylov_configuration_unit)

add_executable(classic_krylov_projection_unit
    classic_krylov_projection_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c)
target_include_directories(classic_krylov_projection_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(classic_krylov_projection_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(classic_krylov_projection_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
add_test(NAME ClassicKrylovProjection_Unit
    COMMAND classic_krylov_projection_unit)

add_executable(update_weight_unit
    update_weight_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/update_weight.c)
target_include_directories(update_weight_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_link_libraries(update_weight_unit m)
add_test(NAME UpdateWeight_Unit COMMAND update_weight_unit)

add_executable(lanczos2_solver_unit
    lanczos2_solver_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/physcal_lanczos2.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/physcal_lanczos.c)
target_include_directories(lanczos2_solver_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_link_libraries(lanczos2_solver_unit ${LAPACK_LIBRARIES} m)
if(MPI_FOUND)
    target_link_libraries(lanczos2_solver_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME Lanczos2Solver_Unit COMMAND lanczos2_solver_unit)

add_executable(krylov_gevp_solver_unit
    krylov_gevp_solver_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_gevp_solver.c)
target_include_directories(krylov_gevp_solver_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(krylov_gevp_solver_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(krylov_gevp_solver_unit ${LAPACK_LIBRARIES} m)
add_test(NAME KrylovGEVPSolver_Unit COMMAND krylov_gevp_solver_unit)

add_executable(power_lanczos_gevp_unit
    power_lanczos_gevp_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_gevp.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_gevp_solver.c)
target_include_directories(power_lanczos_gevp_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_gevp_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
    target_compile_definitions(power_lanczos_gevp_unit PRIVATE
        MVMC_POWER_LANCZOS_FP_CONTRACT_OFF_REQUESTED=1)
endif()
target_link_libraries(power_lanczos_gevp_unit ${LAPACK_LIBRARIES} m)
add_test(NAME PowerLanczosGEVP_Unit COMMAND power_lanczos_gevp_unit)

add_executable(absolute_pfaffian_unit
    absolute_pfaffian_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_include_directories(absolute_pfaffian_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(absolute_pfaffian_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(absolute_pfaffian_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(absolute_pfaffian_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(absolute_pfaffian_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME AbsolutePfaffian_Unit COMMAND absolute_pfaffian_unit)
if(MPI_FOUND)
    add_test(NAME AbsolutePfaffian_Unit_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:absolute_pfaffian_unit>
                ${MPIEXEC_POSTFLAGS})
endif()

add_executable(scaled_pfaffian_unit
    scaled_pfaffian_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(scaled_pfaffian_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(scaled_pfaffian_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(scaled_pfaffian_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(scaled_pfaffian_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(scaled_pfaffian_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(scaled_pfaffian_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME ScaledPfaffian_Unit COMMAND scaled_pfaffian_unit)
if(MPI_FOUND)
    add_test(NAME ScaledPfaffian_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:scaled_pfaffian_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME ScaledPfaffian_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:scaled_pfaffian_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(ScaledPfaffian_Unit_Rank2_mpi PROPERTIES PROCESSORS 2)
    set_tests_properties(ScaledPfaffian_Unit_Rank4_mpi PROPERTIES PROCESSORS 4)
endif()

add_executable(classic_pfaffian_state_unit
    classic_pfaffian_state_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_include_directories(classic_pfaffian_state_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(classic_pfaffian_state_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(classic_pfaffian_state_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(classic_pfaffian_state_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(classic_pfaffian_state_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME ClassicPfaffianState_Unit COMMAND classic_pfaffian_state_unit)

add_executable(classic_pfaffian_collective_unit
    classic_pfaffian_collective_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_include_directories(classic_pfaffian_collective_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(classic_pfaffian_collective_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(classic_pfaffian_collective_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(classic_pfaffian_collective_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(classic_pfaffian_collective_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME ClassicPfaffianCollective_Unit
    COMMAND classic_pfaffian_collective_unit)
if(MPI_FOUND)
    add_test(NAME ClassicPfaffianCollective_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:classic_pfaffian_collective_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME ClassicPfaffianCollective_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:classic_pfaffian_collective_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(ClassicPfaffianCollective_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
    set_tests_properties(ClassicPfaffianCollective_Unit_Rank4_mpi PROPERTIES
        PROCESSORS 4)
endif()

add_executable(classic_pfaffian_sampler_unit
    classic_pfaffian_sampler_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_audit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(classic_pfaffian_sampler_unit PRIVATE
    MVMC_ENABLE_PFAFFIAN_SAMPLER_REFERENCE
    MVMC_ENABLE_PFAFFIAN_STATE_AUDIT)
target_include_directories(classic_pfaffian_sampler_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(classic_pfaffian_sampler_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(classic_pfaffian_sampler_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(classic_pfaffian_sampler_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(classic_pfaffian_sampler_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME ClassicPfaffianSampler_Unit
    COMMAND classic_pfaffian_sampler_unit)
if(MPI_FOUND)
    add_test(NAME ClassicPfaffianSampler_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:classic_pfaffian_sampler_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME ClassicPfaffianSampler_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:classic_pfaffian_sampler_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(ClassicPfaffianSampler_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
    set_tests_properties(ClassicPfaffianSampler_Unit_Rank4_mpi PROPERTIES
        PROCESSORS 4)
endif()

add_executable(classic_pfaffian_legacy_ab_unit
    classic_pfaffian_legacy_ab_unit.c
    classic_pfaffian_legacy_matrix_tu.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_include_directories(classic_pfaffian_legacy_ab_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_link_libraries(classic_pfaffian_legacy_ab_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(classic_pfaffian_legacy_ab_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(classic_pfaffian_legacy_ab_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME ClassicPfaffianLegacyExecutionAB_Unit
    COMMAND classic_pfaffian_legacy_ab_unit)

add_executable(absolute_krylov_reference_unit
    absolute_krylov_reference_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c)
target_compile_definitions(absolute_krylov_reference_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
target_include_directories(absolute_krylov_reference_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(absolute_krylov_reference_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(absolute_krylov_reference_unit m)
add_test(NAME AbsoluteKrylovReference_Unit
    COMMAND absolute_krylov_reference_unit)

add_executable(krylov_fock_proposal_unit
    krylov_fock_proposal_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c)
target_compile_definitions(krylov_fock_proposal_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
target_include_directories(krylov_fock_proposal_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(krylov_fock_proposal_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(krylov_fock_proposal_unit m)
add_test(NAME KrylovFockProposal_Unit
    COMMAND krylov_fock_proposal_unit)

add_executable(bounded_krylov_engine_unit
    bounded_krylov_engine_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(bounded_krylov_engine_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
    MVMC_ENABLE_POWER_LANCZOS_TESTING_HOOKS)
target_include_directories(bounded_krylov_engine_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bounded_krylov_engine_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(bounded_krylov_engine_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(bounded_krylov_engine_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(bounded_krylov_engine_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(bounded_krylov_engine_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME BoundedKrylovEngine_Unit
    COMMAND bounded_krylov_engine_unit)
if(MPI_FOUND)
    add_test(NAME BoundedKrylovEngine_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:bounded_krylov_engine_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME BoundedKrylovEngine_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:bounded_krylov_engine_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(BoundedKrylovEngine_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
    set_tests_properties(BoundedKrylovEngine_Unit_Rank4_mpi PROPERTIES
        PROCESSORS 4)
endif()

add_executable(bounded_krylov_collective_unit
    bounded_krylov_collective_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(bounded_krylov_collective_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(bounded_krylov_collective_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bounded_krylov_collective_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(bounded_krylov_collective_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(bounded_krylov_collective_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(bounded_krylov_collective_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(bounded_krylov_collective_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME BoundedKrylovCollective_Unit
    COMMAND bounded_krylov_collective_unit)
if(MPI_FOUND)
    add_test(NAME BoundedKrylovCollective_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:bounded_krylov_collective_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME BoundedKrylovCollective_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:bounded_krylov_collective_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(BoundedKrylovCollective_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
    set_tests_properties(BoundedKrylovCollective_Unit_Rank4_mpi PROPERTIES
        PROCESSORS 4)
endif()

add_executable(bounded_krylov_integration_unit
    bounded_krylov_integration_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(bounded_krylov_integration_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(bounded_krylov_integration_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bounded_krylov_integration_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(bounded_krylov_integration_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(bounded_krylov_integration_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(bounded_krylov_integration_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(bounded_krylov_integration_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME BoundedKrylovIntegration_Unit
    COMMAND bounded_krylov_integration_unit)
if(MPI_FOUND)
    add_test(NAME BoundedKrylovIntegration_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:bounded_krylov_integration_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME BoundedKrylovIntegration_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:bounded_krylov_integration_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(BoundedKrylovIntegration_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
    set_tests_properties(BoundedKrylovIntegration_Unit_Rank4_mpi PROPERTIES
        PROCESSORS 4)
endif()

add_executable(krylov_positive_guide_unit
    krylov_positive_guide_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(krylov_positive_guide_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(krylov_positive_guide_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(krylov_positive_guide_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(krylov_positive_guide_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(krylov_positive_guide_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(krylov_positive_guide_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(krylov_positive_guide_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME KrylovPositiveGuide_Unit
    COMMAND krylov_positive_guide_unit)

add_executable(krylov_positive_sampler_unit
    krylov_positive_sampler_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(krylov_positive_sampler_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(krylov_positive_sampler_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(krylov_positive_sampler_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(krylov_positive_sampler_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(krylov_positive_sampler_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(krylov_positive_sampler_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(krylov_positive_sampler_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME KrylovPositiveSampler_Unit
    COMMAND krylov_positive_sampler_unit)
if(MPI_FOUND)
    add_test(NAME KrylovPositiveSampler_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                $<TARGET_FILE:krylov_positive_sampler_unit>)
    set_tests_properties(KrylovPositiveSampler_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
    add_test(NAME KrylovPositiveSampler_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                $<TARGET_FILE:krylov_positive_sampler_unit>)
    set_tests_properties(KrylovPositiveSampler_Unit_Rank4_mpi PROPERTIES
        PROCESSORS 4)
endif()

add_executable(krylov_final_state_unit
    krylov_final_state_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(krylov_final_state_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
    MVMC_ENABLE_POWER_LANCZOS_P5_TESTING)
target_include_directories(krylov_final_state_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(krylov_final_state_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(krylov_final_state_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(krylov_final_state_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(krylov_final_state_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(krylov_final_state_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME KrylovFinalState_Unit
    COMMAND krylov_final_state_unit)
if(MPI_FOUND)
    add_test(NAME KrylovFinalState_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                $<TARGET_FILE:krylov_final_state_unit>)
    set_tests_properties(KrylovFinalState_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
endif()

add_executable(krylov_final_state_sensitivity_unit
    krylov_final_state_sensitivity_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state_sensitivity.c)
target_compile_definitions(krylov_final_state_sensitivity_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
    MVMC_ENABLE_POWER_LANCZOS_P5_TESTING)
target_include_directories(krylov_final_state_sensitivity_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(krylov_final_state_sensitivity_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(krylov_final_state_sensitivity_unit m)
if(MPI_FOUND)
    target_link_libraries(krylov_final_state_sensitivity_unit
        ${MPI_C_LIBRARIES})
endif()
add_test(NAME KrylovFinalStateSensitivity_Unit
    COMMAND krylov_final_state_sensitivity_unit)
if(MPI_FOUND)
    add_test(NAME KrylovFinalStateSensitivity_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                $<TARGET_FILE:krylov_final_state_sensitivity_unit>)
    set_tests_properties(
        KrylovFinalStateSensitivity_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
endif()

add_executable(krylov_final_state_chain_unit
    krylov_final_state_chain_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state_chain.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(krylov_final_state_chain_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
    MVMC_ENABLE_POWER_LANCZOS_P5_TESTING)
target_include_directories(krylov_final_state_chain_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(krylov_final_state_chain_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(krylov_final_state_chain_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(krylov_final_state_chain_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(krylov_final_state_chain_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(krylov_final_state_chain_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME KrylovFinalStateChain_Unit
    COMMAND krylov_final_state_chain_unit)
if(MPI_FOUND)
    add_test(NAME KrylovFinalStateChain_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                $<TARGET_FILE:krylov_final_state_chain_unit>)
    set_tests_properties(KrylovFinalStateChain_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
    add_test(NAME KrylovFinalStateChain_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                $<TARGET_FILE:krylov_final_state_chain_unit>)
    set_tests_properties(KrylovFinalStateChain_Unit_Rank4_mpi PROPERTIES
        PROCESSORS 4)
endif()

add_executable(krylov_final_state_chain_driver
    krylov_final_state_chain_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state_chain.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_streaming_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(krylov_final_state_chain_driver PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
    MVMC_ENABLE_POWER_LANCZOS_P5_TESTING)
target_include_directories(krylov_final_state_chain_driver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(krylov_final_state_chain_driver PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(krylov_final_state_chain_driver pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(krylov_final_state_chain_driver stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(krylov_final_state_chain_driver blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(krylov_final_state_chain_driver ${MPI_C_LIBRARIES})
endif()
add_test(NAME KrylovFinalStateChainDriver_Smoke
    COMMAND krylov_final_state_chain_driver
            two_state_support_bridge 256 16 64
            0x5035434453454541 0x5035435354524541 123)
if(MPI_FOUND)
    add_test(NAME KrylovFinalStateChainDriver_Smoke_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                $<TARGET_FILE:krylov_final_state_chain_driver>
                complex_ring4_near_zero 256 16 64
                0x5035434453454541 0x5035435354524541 123)
    set_tests_properties(
        KrylovFinalStateChainDriver_Smoke_Rank2_mpi PROPERTIES
        PROCESSORS 2)
endif()

add_executable(krylov_streaming_statistics_unit
    krylov_streaming_statistics_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_streaming_statistics.c)
target_compile_definitions(krylov_streaming_statistics_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
target_include_directories(krylov_streaming_statistics_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(krylov_streaming_statistics_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(krylov_streaming_statistics_unit m)
add_test(NAME KrylovStreamingStatistics_Unit
    COMMAND krylov_streaming_statistics_unit)

add_executable(krylov_matrix_measurement_unit
    krylov_matrix_measurement_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_matrix_measurement.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_streaming_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(krylov_matrix_measurement_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(krylov_matrix_measurement_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(krylov_matrix_measurement_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(krylov_matrix_measurement_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(krylov_matrix_measurement_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(krylov_matrix_measurement_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(krylov_matrix_measurement_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME KrylovMatrixMeasurement_Unit
    COMMAND krylov_matrix_measurement_unit)

add_executable(classic_krylov_amplitude_unit
    classic_krylov_amplitude_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_amplitude.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(classic_krylov_amplitude_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(classic_krylov_amplitude_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(classic_krylov_amplitude_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(classic_krylov_amplitude_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(classic_krylov_amplitude_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(classic_krylov_amplitude_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME ClassicKrylovAmplitude_Unit
    COMMAND classic_krylov_amplitude_unit)
if(MPI_FOUND)
    add_test(NAME ClassicKrylovAmplitude_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:classic_krylov_amplitude_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME ClassicKrylovAmplitude_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:classic_krylov_amplitude_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(ClassicKrylovAmplitude_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
    set_tests_properties(ClassicKrylovAmplitude_Unit_Rank4_mpi PROPERTIES
        PROCESSORS 4)
endif()

# Production-only anchor: intentionally omit the reference and testing-hook
# macros so this target proves the classic adapters and scaled callback are
# reachable through the bounded engine exactly as they are in vmc.out.
add_executable(classic_krylov_production_anchor_unit
    classic_krylov_production_anchor_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_amplitude.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_classic_bridge.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_sha256.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_rng.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_chain_rng.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_snapshot.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_block_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_census.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_action.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_evaluator.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_blocks.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_observable_confirmation.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_primitive_trace.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_production_session.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_production_session_prepare.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_production_session_stages.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_stabilization_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_stabilization_output.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_corrected_dispatch.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_json_writer.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_output_transaction.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_streaming_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_matrix_measurement.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_gevp_solver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_gevp.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state_chain.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(classic_krylov_production_anchor_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
    MVMC_ENABLE_POWER_LANCZOS_P5_CORE)
target_include_directories(classic_krylov_production_anchor_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(classic_krylov_production_anchor_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(classic_krylov_production_anchor_unit
    pfapack ltl2inv ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(classic_krylov_production_anchor_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(classic_krylov_production_anchor_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(classic_krylov_production_anchor_unit
        ${MPI_C_LIBRARIES})
endif()
add_test(NAME ClassicKrylovProductionAnchor_Unit
    COMMAND classic_krylov_production_anchor_unit)
if(MPI_FOUND)
    add_test(NAME ClassicKrylovProductionAnchor_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:classic_krylov_production_anchor_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME ClassicKrylovProductionAnchor_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:classic_krylov_production_anchor_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(
        ClassicKrylovProductionAnchor_Unit_Rank2_mpi PROPERTIES PROCESSORS 2)
    set_tests_properties(
        ClassicKrylovProductionAnchor_Unit_Rank4_mpi PROPERTIES PROCESSORS 4)
endif()

add_executable(classic_krylov_model_unit
    classic_krylov_model_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c)
target_compile_definitions(classic_krylov_model_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
target_include_directories(classic_krylov_model_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(classic_krylov_model_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(classic_krylov_model_unit m)
if(MPI_FOUND)
    target_link_libraries(classic_krylov_model_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME ClassicKrylovModel_Unit
    COMMAND classic_krylov_model_unit)
if(MPI_FOUND)
    add_test(NAME ClassicKrylovModel_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:classic_krylov_model_unit>
                ${MPIEXEC_POSTFLAGS})
    add_test(NAME ClassicKrylovModel_Unit_Rank4_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:classic_krylov_model_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(ClassicKrylovModel_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
    set_tests_properties(ClassicKrylovModel_Unit_Rank4_mpi PROPERTIES
        PROCESSORS 4)
endif()

add_executable(absolute_krylov_classic_ed_driver
    absolute_krylov_classic_ed_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_amplitude.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(absolute_krylov_classic_ed_driver PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
target_include_directories(absolute_krylov_classic_ed_driver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(absolute_krylov_classic_ed_driver PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(absolute_krylov_classic_ed_driver pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(absolute_krylov_classic_ed_driver blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(absolute_krylov_classic_ed_driver ${MPI_C_LIBRARIES})
endif()

# One-time P6-C1 Cartesian evidence driver.  Deliberately not registered as
# a permanent CTest; the versioned evidence generator invokes it exactly once
# after independent pre-result review.
add_executable(power_lanczos_p6c1_cartesian_driver
    power_lanczos_p6c1_cartesian_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_amplitude.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_gevp_solver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_gevp.c)
target_compile_definitions(power_lanczos_p6c1_cartesian_driver PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(power_lanczos_p6c1_cartesian_driver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_p6c1_cartesian_driver PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_p6c1_cartesian_driver pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(power_lanczos_p6c1_cartesian_driver stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(power_lanczos_p6c1_cartesian_driver blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(power_lanczos_p6c1_cartesian_driver
        ${MPI_C_LIBRARIES})
endif()

# Focused pilot composition test.  This is a distinct executable with an
# exact nine-case grammar; the official driver above keeps its 2016-case gate.
add_executable(power_lanczos_p6c1_preflight
    power_lanczos_p6c1_cartesian_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_amplitude.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_gevp_solver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_gevp.c)
target_compile_definitions(power_lanczos_p6c1_preflight PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE MVMC_P6C1_PREFLIGHT)
target_include_directories(power_lanczos_p6c1_preflight PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_p6c1_preflight PRIVATE
        -Wall -Wextra -Wpedantic -Werror -fno-fast-math -ffp-contract=off)
endif()
target_link_libraries(power_lanczos_p6c1_preflight pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(power_lanczos_p6c1_preflight stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(power_lanczos_p6c1_preflight blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(power_lanczos_p6c1_preflight ${MPI_C_LIBRARIES})
endif()
set(POWER_LANCZOS_P6C1_PREFLIGHT_INPUT
    ${CMAKE_BINARY_DIR}/test/power_lanczos_p6c1_preflight_input.txt)
add_test(NAME PowerLanczosP6C1_PreflightPrepare
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/test/python/power_lanczos_p6c1_cartesian_evidence.py
            prepare-preflight ${POWER_LANCZOS_P6C1_PREFLIGHT_INPUT})
set_tests_properties(PowerLanczosP6C1_PreflightPrepare PROPERTIES
    FIXTURES_SETUP PowerLanczosP6C1_PreflightInput)
add_test(NAME PowerLanczosP6C1_Preflight
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_SOURCE_DIR}/test/python/power_lanczos_p6c1_cartesian_evidence.py
            run-preflight --expected-world-size 1 --
            $<TARGET_FILE:power_lanczos_p6c1_preflight>
            ${POWER_LANCZOS_P6C1_PREFLIGHT_INPUT})
set_tests_properties(PowerLanczosP6C1_Preflight PROPERTIES
    FIXTURES_REQUIRED PowerLanczosP6C1_PreflightInput
    ENVIRONMENT "OMP_NUM_THREADS=1")
if(MPI_FOUND)
    add_test(NAME PowerLanczosP6C1_Preflight_mpi
        COMMAND ${PYTHON_EXECUTABLE}
                ${CMAKE_SOURCE_DIR}/test/python/power_lanczos_p6c1_cartesian_evidence.py
                run-preflight --expected-world-size 2 --
                ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:power_lanczos_p6c1_preflight>
                ${POWER_LANCZOS_P6C1_PREFLIGHT_INPUT}
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(PowerLanczosP6C1_Preflight_mpi PROPERTIES
        FIXTURES_REQUIRED PowerLanczosP6C1_PreflightInput
        ENVIRONMENT "OMP_NUM_THREADS=1")
endif()

add_executable(absolute_krylov_profile_driver
    absolute_krylov_profile_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_amplitude.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(absolute_krylov_profile_driver PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
target_include_directories(absolute_krylov_profile_driver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(absolute_krylov_profile_driver PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(absolute_krylov_profile_driver pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(absolute_krylov_profile_driver blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(absolute_krylov_profile_driver ${MPI_C_LIBRARIES})
endif()
add_test(NAME AbsoluteKrylovProfile_Smoke
    COMMAND absolute_krylov_profile_driver 4 1 3 1)
add_test(NAME AbsoluteKrylovProfile_InvalidOddSite
    COMMAND absolute_krylov_profile_driver 5 1 1 1)
add_test(NAME AbsoluteKrylovProfile_InvalidQP
    COMMAND absolute_krylov_profile_driver 4 2 1 1)
add_test(NAME AbsoluteKrylovProfile_InvalidSampleCount
    COMMAND absolute_krylov_profile_driver 4 1 37 1)
set_tests_properties(
    AbsoluteKrylovProfile_InvalidOddSite
    AbsoluteKrylovProfile_InvalidQP
    AbsoluteKrylovProfile_InvalidSampleCount
    PROPERTIES WILL_FAIL TRUE)
if(MPI_FOUND)
    add_test(NAME AbsoluteKrylovProfile_ResourceLimit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                $<TARGET_FILE:absolute_krylov_profile_driver> 16 4 1 1)
    set_tests_properties(
        AbsoluteKrylovProfile_ResourceLimit_Rank2_mpi PROPERTIES
        WILL_FAIL TRUE TIMEOUT 30 PROCESSORS 2
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
endif()

add_executable(bounded_krylov_profile_driver
    bounded_krylov_profile_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(bounded_krylov_profile_driver PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(bounded_krylov_profile_driver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bounded_krylov_profile_driver PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(bounded_krylov_profile_driver pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(bounded_krylov_profile_driver stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(bounded_krylov_profile_driver blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(bounded_krylov_profile_driver ${MPI_C_LIBRARIES})
endif()
add_test(NAME BoundedKrylovProfile_Smoke
    COMMAND bounded_krylov_profile_driver 4 1 3 0 1)
add_test(NAME BoundedKrylovProfile_CacheSmoke
    COMMAND bounded_krylov_profile_driver 4 4 3 4096 1)
add_test(NAME BoundedKrylovProfile_SessionSmoke
    COMMAND bounded_krylov_profile_driver 4 1 3 4096 1 session-profile)
add_test(NAME BoundedKrylovProfile_InvalidOddSite
    COMMAND bounded_krylov_profile_driver 5 1 1 0 1)
add_test(NAME BoundedKrylovProfile_InvalidQP
    COMMAND bounded_krylov_profile_driver 4 2 1 0 1)
add_test(NAME BoundedKrylovProfile_InvalidCache
    COMMAND bounded_krylov_profile_driver 4 1 1 4294967297 1)
set_tests_properties(
    BoundedKrylovProfile_InvalidOddSite
    BoundedKrylovProfile_InvalidQP
    BoundedKrylovProfile_InvalidCache
    PROPERTIES WILL_FAIL TRUE)
if(MPI_FOUND)
    add_test(NAME BoundedKrylovProfile_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                $<TARGET_FILE:bounded_krylov_profile_driver> 4 4 3 4096 1)
    set_tests_properties(BoundedKrylovProfile_Rank2_mpi PROPERTIES
        PROCESSORS 2
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
endif()

add_executable(bounded_krylov_markov_driver
    bounded_krylov_markov_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_matrix_measurement.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_streaming_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(bounded_krylov_markov_driver PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(bounded_krylov_markov_driver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bounded_krylov_markov_driver PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(bounded_krylov_markov_driver pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(bounded_krylov_markov_driver stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(bounded_krylov_markov_driver blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(bounded_krylov_markov_driver ${MPI_C_LIBRARIES})
endif()
add_test(NAME BoundedKrylovMarkov_Smoke
    COMMAND bounded_krylov_markov_driver
            4 1 128 0 1e-4 0x5034535F52305631)
set_tests_properties(BoundedKrylovMarkov_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME BoundedKrylovMarkovShell_Smoke
    COMMAND bounded_krylov_markov_driver
            4 1 128 0 1e-2 0x50345334534D4B31 shell 1 8 1 2 7)
set_tests_properties(BoundedKrylovMarkovShell_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME BoundedKrylovMarkovP4S9Session_Smoke
    COMMAND bounded_krylov_markov_driver
            4 1 128 4194304 1e-2 0x5034533952305631 0 1 7 session)
set_tests_properties(BoundedKrylovMarkovP4S9Session_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
if(MPI_FOUND)
    add_test(NAME BoundedKrylovMarkovShell_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                $<TARGET_FILE:bounded_krylov_markov_driver>
                4 1 128 0 1e-2 0x50345334534D4B31 shell 1 8 1 2 7)
    add_test(NAME BoundedKrylovMarkovShell_Rank4_EmptyQP_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                $<TARGET_FILE:bounded_krylov_markov_driver>
                4 1 128 0 1e-2 0x50345334534D4B31 shell 1 8 1 2 7)
    add_test(NAME BoundedKrylovMarkovP4S9Session_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                $<TARGET_FILE:bounded_krylov_markov_driver>
                4 1 128 4194304 1e-2 0x5034533952305631 0 1 7 session)
    add_test(NAME BoundedKrylovMarkovP4S9Session_Rank4_EmptyQP_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                $<TARGET_FILE:bounded_krylov_markov_driver>
                4 1 128 4194304 1e-2 0x5034533952305631 0 1 7 session)
    set_tests_properties(
        BoundedKrylovMarkovShell_Rank2_mpi
        BoundedKrylovMarkovShell_Rank4_EmptyQP_mpi
        BoundedKrylovMarkovP4S9Session_Rank2_mpi
        BoundedKrylovMarkovP4S9Session_Rank4_EmptyQP_mpi
        PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    set_tests_properties(
        BoundedKrylovMarkovShell_Rank2_mpi
        BoundedKrylovMarkovP4S9Session_Rank2_mpi PROPERTIES
        PROCESSORS 2)
    set_tests_properties(
        BoundedKrylovMarkovShell_Rank4_EmptyQP_mpi
        BoundedKrylovMarkovP4S9Session_Rank4_EmptyQP_mpi
        PROPERTIES PROCESSORS 4)
endif()

add_executable(p4f_matrix_reparse_driver
    p4f_matrix_reparse_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_matrix_measurement.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_streaming_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_gevp_solver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
add_dependencies(p4f_matrix_reparse_driver bounded_krylov_markov_driver)
target_compile_definitions(p4f_matrix_reparse_driver PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(p4f_matrix_reparse_driver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(p4f_matrix_reparse_driver PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(p4f_matrix_reparse_driver pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(p4f_matrix_reparse_driver stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(p4f_matrix_reparse_driver blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(p4f_matrix_reparse_driver ${MPI_C_LIBRARIES})
endif()
add_test(NAME P4FMatrixReparse_Smoke
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/python/runtest_p4f_matrix_reparse.py
            $<TARGET_FILE:bounded_krylov_markov_driver>
            $<TARGET_FILE:p4f_matrix_reparse_driver>)
set_tests_properties(P4FMatrixReparse_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")

add_executable(bounded_krylov_markov_pilot_driver
    bounded_krylov_markov_pilot_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_matrix_measurement.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_streaming_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(bounded_krylov_markov_pilot_driver PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(bounded_krylov_markov_pilot_driver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bounded_krylov_markov_pilot_driver PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(bounded_krylov_markov_pilot_driver pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(bounded_krylov_markov_pilot_driver stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(bounded_krylov_markov_pilot_driver blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(bounded_krylov_markov_pilot_driver
        ${MPI_C_LIBRARIES})
endif()
add_test(NAME BoundedKrylovMarkovPilot_Smoke
    COMMAND bounded_krylov_markov_pilot_driver 4 1 256 0)
set_tests_properties(BoundedKrylovMarkovPilot_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME BoundedKrylovMarkovShellPilot_Smoke
    COMMAND bounded_krylov_markov_pilot_driver 4 1 256 0 shell)
set_tests_properties(BoundedKrylovMarkovShellPilot_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME BoundedKrylovMarkovFlatGuidePilot_Smoke
    COMMAND bounded_krylov_markov_pilot_driver 4 1 256 0 flat)
set_tests_properties(BoundedKrylovMarkovFlatGuidePilot_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME BoundedKrylovMarkovSurrogatePilot_Smoke
    COMMAND bounded_krylov_markov_pilot_driver 4 1 256 0 surrogate)
set_tests_properties(BoundedKrylovMarkovSurrogatePilot_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME BoundedKrylovMarkovPartialGuidePilot_Smoke
    COMMAND bounded_krylov_markov_pilot_driver 4 1 256 0 partial)
set_tests_properties(BoundedKrylovMarkovPartialGuidePilot_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME BoundedKrylovMarkovPartialNeighborPilot_Smoke
    COMMAND bounded_krylov_markov_pilot_driver
            4 1 256 0 partial-neighbor)
set_tests_properties(BoundedKrylovMarkovPartialNeighborPilot_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME BoundedKrylovMarkovPartialCallbackPilot_Smoke
    COMMAND bounded_krylov_markov_pilot_driver
            4 1 256 0 partial-callback)
set_tests_properties(BoundedKrylovMarkovPartialCallbackPilot_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME BoundedKrylovMarkovLongDirectPilot_Smoke
    COMMAND bounded_krylov_markov_pilot_driver 4 1 256 0 long-direct)
set_tests_properties(BoundedKrylovMarkovLongDirectPilot_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")

add_executable(bounded_krylov_markov_timing_driver
    bounded_krylov_markov_timing_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_matrix_measurement.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_numeric_evidence.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_streaming_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(bounded_krylov_markov_timing_driver PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
target_include_directories(bounded_krylov_markov_timing_driver PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bounded_krylov_markov_timing_driver PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(bounded_krylov_markov_timing_driver pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(bounded_krylov_markov_timing_driver stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(bounded_krylov_markov_timing_driver blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(bounded_krylov_markov_timing_driver
        ${MPI_C_LIBRARIES})
endif()
add_test(NAME BoundedKrylovMarkovPartialCallbackTiming_Smoke
    COMMAND bounded_krylov_markov_timing_driver
            4 1 0 128 partial-callback)
set_tests_properties(BoundedKrylovMarkovPartialCallbackTiming_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME BoundedKrylovMarkovLongDirectSessionTiming_Smoke
    COMMAND bounded_krylov_markov_timing_driver
            4 1 4194304 16 long-direct-session)
set_tests_properties(BoundedKrylovMarkovLongDirectSessionTiming_Smoke PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")

add_subdirectory(python)
