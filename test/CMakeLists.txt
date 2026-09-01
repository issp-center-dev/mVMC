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

add_executable(gc_size_unit gc_size_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/gc_size.c)
target_include_directories(gc_size_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
if(MPI_FOUND)
    target_link_libraries(gc_size_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME GCSize_Unit COMMAND gc_size_unit)
foreach(gc_size_death_mode mul_overflow add_overflow negative int_overflow)
    add_test(NAME GCSize_Unit_${gc_size_death_mode}
        COMMAND gc_size_unit ${gc_size_death_mode})
    set_tests_properties(GCSize_Unit_${gc_size_death_mode}
        PROPERTIES WILL_FAIL TRUE)
endforeach()

add_executable(gc_config_unit gc_config_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/gc_config.c)
target_include_directories(gc_config_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_link_libraries(gc_config_unit m)
add_test(NAME GCConfig_Unit COMMAND gc_config_unit)

add_executable(gc_parameter_unit gc_parameter_unit.c
    ${CMAKE_SOURCE_DIR}/src/sfmt/SFMT.c)
target_include_directories(gc_parameter_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include
    ${CMAKE_SOURCE_DIR}/src/sfmt)
target_link_libraries(gc_parameter_unit m)
if(MPI_FOUND)
    target_link_libraries(gc_parameter_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME GCParameter_Unit COMMAND gc_parameter_unit)

add_executable(gc_matrix_unit gc_matrix_unit.c)
target_include_directories(gc_matrix_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_link_libraries(gc_matrix_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(gc_matrix_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(gc_matrix_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME GCMatrix_Unit COMMAND gc_matrix_unit)

add_executable(gc_pfupdate_unit gc_pfupdate_unit.c)
target_include_directories(gc_pfupdate_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_link_libraries(gc_pfupdate_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(gc_pfupdate_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(gc_pfupdate_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME GCPfUpdate_Unit COMMAND gc_pfupdate_unit)

add_executable(gc_sampler_unit gc_sampler_unit.c
    ${CMAKE_SOURCE_DIR}/src/sfmt/SFMT.c)
target_include_directories(gc_sampler_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include
    ${CMAKE_SOURCE_DIR}/src/sfmt)
target_link_libraries(gc_sampler_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(gc_sampler_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(gc_sampler_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME GCSampler_Unit COMMAND gc_sampler_unit)
add_test(NAME GCSampler_Unit_DebugRebuild COMMAND gc_sampler_unit)
set_tests_properties(GCSampler_Unit_DebugRebuild PROPERTIES
    ENVIRONMENT "MVMC_GC_DEBUG_REBUILD_INTERVAL=1")
if(MPI_FOUND)
    add_test(NAME GCSampler_RebuildStatus_MPI
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:gc_sampler_unit>
                --collective-only ${MPIEXEC_POSTFLAGS})
    set_tests_properties(GCSampler_RebuildStatus_MPI PROPERTIES
        PROCESSORS 2 TIMEOUT 30)
endif()

add_executable(gc_slater_diff_unit gc_slater_diff_unit.c)
target_include_directories(gc_slater_diff_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_link_libraries(gc_slater_diff_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(gc_slater_diff_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(gc_slater_diff_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME GCSlaterDiff_Unit COMMAND gc_slater_diff_unit)

add_executable(gc_green_unit gc_green_unit.c)
target_include_directories(gc_green_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_link_libraries(gc_green_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(Require_BLIS)
    target_link_libraries(gc_green_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(gc_green_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME GCGreen_Unit COMMAND gc_green_unit)
add_test(NAME GCGreen_Unit_OMP4 COMMAND gc_green_unit)
set_tests_properties(GCGreen_Unit_OMP4 PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=4")

add_executable(nbody_operator_unit
    nbody_operator_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/nbody_operator.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/backflow_nbody_scratch.c)
target_include_directories(nbody_operator_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
add_test(NAME NBodyOperator_Unit COMMAND nbody_operator_unit)

add_executable(lanczos2_contract_unit
    lanczos2_contract_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/lanczos2_contract.c)
target_include_directories(lanczos2_contract_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
add_test(NAME Lanczos2Contract_Unit COMMAND lanczos2_contract_unit)

add_executable(update_weight_unit
    update_weight_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/update_weight.c)
target_include_directories(update_weight_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_link_libraries(update_weight_unit m)
add_test(NAME UpdateWeight_Unit COMMAND update_weight_unit)

add_executable(real_rbm_oracle_unit
    real_rbm_oracle_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/rbm.c)
target_include_directories(real_rbm_oracle_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_link_libraries(real_rbm_oracle_unit m)
add_test(NAME RealRBM_PrimitiveOracle COMMAND real_rbm_oracle_unit)
add_test(NAME RealRBM_PrimitiveOracle_omp COMMAND real_rbm_oracle_unit)
set_tests_properties(RealRBM_PrimitiveOracle PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1")
set_tests_properties(RealRBM_PrimitiveOracle_omp PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=4")

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
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
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

add_executable(classic_krylov_amplitude_unit
    classic_krylov_amplitude_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_amplitude.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c)
target_compile_definitions(classic_krylov_amplitude_unit PRIVATE
    MVMC_ENABLE_ABSOLUTE_KRYLOV_REFERENCE)
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

add_executable(power_lanczos_stabilized_unit
    power_lanczos_stabilized_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/absolute_pfaffian.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_matrix.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_pfaffian_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_amplitude.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_engine.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/bounded_krylov_collective.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_proposal.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_guide.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_positive_sampler.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_matrix_measurement.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_streaming_statistics.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_gevp_solver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_final_state_chain.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_classic_bridge.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/power_lanczos_stabilized.c)
target_include_directories(power_lanczos_stabilized_unit PRIVATE
    ${CMAKE_SOURCE_DIR}/src/mVMC/include)
target_compile_definitions(power_lanczos_stabilized_unit PRIVATE
    MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE
    MVMC_ENABLE_POWER_LANCZOS_P5_CORE)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(power_lanczos_stabilized_unit PRIVATE
        -Wall -Wextra -Wpedantic -Werror)
endif()
target_link_libraries(power_lanczos_stabilized_unit pfapack ltl2inv
    ${LAPACK_LIBRARIES} m)
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(power_lanczos_stabilized_unit stdc++)
endif()
if(Require_BLIS)
    target_link_libraries(power_lanczos_stabilized_unit blis pthread)
endif()
if(MPI_FOUND)
    target_link_libraries(power_lanczos_stabilized_unit ${MPI_C_LIBRARIES})
endif()
add_test(NAME PowerLanczosStabilized_Unit
    COMMAND power_lanczos_stabilized_unit)
if(MPI_FOUND)
    add_test(NAME PowerLanczosStabilized_Unit_Rank2_mpi
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:power_lanczos_stabilized_unit>
                ${MPIEXEC_POSTFLAGS})
    set_tests_properties(PowerLanczosStabilized_Unit_Rank2_mpi PROPERTIES
        PROCESSORS 2)
endif()

add_executable(classic_krylov_model_unit
    classic_krylov_model_unit.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_model.c)
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
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_amplitude.c
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

add_executable(absolute_krylov_profile_driver
    absolute_krylov_profile_driver.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/krylov_fock_reference.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_configuration.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_projection.c
    ${CMAKE_SOURCE_DIR}/src/mVMC/classic_krylov_amplitude.c
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
