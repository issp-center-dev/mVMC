file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_mode1.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_mpi.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_expert.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_invalid_expert.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_sr.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_nstore_nsplit.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_UHF.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_UHF_InterAll.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_expert_phys.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_expert_mpi.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_doublon_only.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_reweight.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_rescale.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/backflow_def_helper.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_bf.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_bf_fsz.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_lanczos_expert.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_lanczos_sector_ed.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_lanczos2_step1_compat.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_lanczos2_real_ed.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_lanczos2_complex_ed.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_lanczos2_heisenberg_ed.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_lanczos2_projected_ed.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_lanczos2_projection_smoke.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_absolute_krylov_classic_ed.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_absolute_krylov_production_anchor.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_absolute_krylov_profile.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/generate_absolute_krylov_p3_evidence.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_absolute_krylov_p3_evidence.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_absolute_krylov_p3_pipeline.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/absolute_krylov_p3_manifest.json DESTINATION ${CMAKE_BINARY_DIR}/test/python/data)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/absolute_krylov_p3_measurement_policy.json DESTINATION ${CMAKE_BINARY_DIR}/test/python/data)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_nbodyg.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/nbodyg_exact_oracle.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/nbodyinterall_exact_oracle.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_nbodyinterall_invalid.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_update_weight.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)

function(add_python_vmc_test model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test)

function(add_python_vmc_test_mode1 model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_mode1.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_mode1)

function(add_python_vmc_test_mpi model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_mpi.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_mpi)

function(add_python_vmc_test_expert model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_expert.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_expert)

function(add_python_vmc_test_invalid_expert model expected_error)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py ${model} "${expected_error}")
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_invalid_expert)

function(add_python_vmc_test_invalid_expert_with_updates test_name model expected_error)
    add_test(NAME ${test_name}
        COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py
        ${model} "${expected_error}" ${ARGN})
    set_tests_properties(${test_name} PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_invalid_expert_with_updates)

function(add_python_vmc_test_expert_phys model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_expert_phys.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_expert_phys)

function(add_python_vmc_test_nbodyg model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_nbodyg.py ${model})
    set_tests_properties(${model} PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python"
        RESOURCE_LOCK "nbodyg_${model}_workdir")
endfunction(add_python_vmc_test_nbodyg)

function(add_python_vmc_test_nbodyg_mpi test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_nbodyg.py ${model})
    set_tests_properties(${test_name} PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2"
        RESOURCE_LOCK "nbodyg_${model}_workdir")
endfunction(add_python_vmc_test_nbodyg_mpi)

function(add_python_nbodyg_exact_oracle test_name)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyg_exact_oracle.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1")
endfunction(add_python_nbodyg_exact_oracle)

function(add_python_nbodyg_generated_check test_name omp_threads)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyg_exact_oracle.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_nbodyg_generated_check)

function(add_python_nbodyg_config_oracle test_name)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyg_exact_oracle.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_BF_NBODY_ORACLE_DUMP=backflow_nbody_oracle.dat;MVMC_BF_FSZ_GREEN_REBUILD_CHECK=1;MVMC_BF_FSZ_PF_UPDATE_CHECK=1")
endfunction(add_python_nbodyg_config_oracle)

function(add_python_nbodyg_config_oracle_omp test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyg_exact_oracle.py ${model})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;MVMC_BF_NBODY_ORACLE_DUMP=backflow_nbody_oracle.dat;MVMC_BF_NBODY_TEST_WORK_SUFFIX=_${test_name};MVMC_BF_FSZ_GREEN_REBUILD_CHECK=1;MVMC_BF_FSZ_PF_UPDATE_CHECK=1")
endfunction(add_python_nbodyg_config_oracle_omp)

function(add_python_nbodyg_config_oracle_mpi test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyg_exact_oracle.py ${model})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_BF_NBODY_ORACLE_DUMP=backflow_nbody_oracle.dat;MVMC_BF_NBODY_TEST_WORK_SUFFIX=_${test_name};MVMC_BF_FSZ_GREEN_REBUILD_CHECK=1;MVMC_BF_FSZ_PF_UPDATE_CHECK=1")
endfunction(add_python_nbodyg_config_oracle_mpi)

function(add_python_nbodyg_invalid test_name expected_error)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyg_exact_oracle.py ${test_name} --expect-error "${expected_error}")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1")
endfunction(add_python_nbodyg_invalid)

function(add_python_nbodyg_invalid_mpi test_name expected_error)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyg_exact_oracle.py ${test_name} --expect-error "${expected_error}")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_MPI_PROCS=2")
endfunction(add_python_nbodyg_invalid_mpi)

function(add_python_nbodyinterall_exact_oracle test_name)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyinterall_exact_oracle.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1"
        RESOURCE_LOCK "nbodyinterall_${test_name}_workdir")
endfunction(add_python_nbodyinterall_exact_oracle)

function(add_python_nbodyinterall_generated_check test_name omp_threads)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyinterall_exact_oracle.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_nbodyinterall_generated_check)

function(add_python_nbodyinterall_config_oracle test_name)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyinterall_exact_oracle.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_BF_NBODY_ORACLE_DUMP=backflow_nbody_oracle.dat;MVMC_BF_FSZ_GREEN_REBUILD_CHECK=1;MVMC_BF_FSZ_PF_UPDATE_CHECK=1")
endfunction(add_python_nbodyinterall_config_oracle)

function(add_python_nbodyinterall_config_oracle_omp test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyinterall_exact_oracle.py ${model})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;MVMC_BF_NBODY_ORACLE_DUMP=backflow_nbody_oracle.dat;MVMC_BF_NBODY_TEST_WORK_SUFFIX=_${test_name};MVMC_BF_FSZ_GREEN_REBUILD_CHECK=1;MVMC_BF_FSZ_PF_UPDATE_CHECK=1")
endfunction(add_python_nbodyinterall_config_oracle_omp)

function(add_python_nbodyinterall_config_oracle_mpi test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyinterall_exact_oracle.py ${model})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_BF_NBODY_ORACLE_DUMP=backflow_nbody_oracle.dat;MVMC_BF_NBODY_TEST_WORK_SUFFIX=_${test_name};MVMC_BF_FSZ_GREEN_REBUILD_CHECK=1;MVMC_BF_FSZ_PF_UPDATE_CHECK=1")
endfunction(add_python_nbodyinterall_config_oracle_mpi)

function(add_python_nbodyinterall_mpi test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} nbodyinterall_exact_oracle.py ${model})
    set_tests_properties(${test_name} PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_MPI_PROCS=2"
        RESOURCE_LOCK "nbodyinterall_${model}_workdir")
endfunction(add_python_nbodyinterall_mpi)

function(add_python_nbodyinterall_invalid test_name expected_error)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_nbodyinterall_invalid.py ${test_name} "${expected_error}")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_nbodyinterall_invalid)

function(add_python_vmc_test_expert_phys_mpi test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_expert_phys.py ${model})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2")
endfunction(add_python_vmc_test_expert_phys_mpi)

function(add_python_vmc_test_expert_mpi model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_expert_mpi.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_expert_mpi)

function(add_python_vmc_test_doublon_only model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_doublon_only.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_doublon_only)

function(add_python_vmc_test_lanczos_expert model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos_expert.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_lanczos_expert)

function(add_python_vmc_test_lanczos_sector_ed test_name)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos_sector_ed.py)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_lanczos_sector_ed)

function(add_python_vmc_test_doublon_only_mpi test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_doublon_only.py ${model})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2")
endfunction(add_python_vmc_test_doublon_only_mpi)

function(add_python_vmc_test_reweight test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_reweight.py ${model})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_reweight)


function(add_python_vmc_test_sr model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_sr.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_sr)

function(add_python_vmc_test_rescale test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_rescale.py ${model})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_rescale)

function(add_python_backflow_helper_test test_name)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} backflow_def_helper.py --self-test)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_backflow_helper_test)

function(add_python_vmc_test_bf test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_bf)

function(add_python_vmc_test_bf_lanczos test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --lanczos-mode 1 --compare-no-bf-lanczos)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_BF_LANCZOS_STATE_CHECK=1;MVMC_LANCZOS_ORACLE_DUMP=1;MVMC_BF_TEST_WORK_SUFFIX=_${test_name}")
endfunction(add_python_vmc_test_bf_lanczos)

function(add_python_vmc_test_bf_lanczos_mpi test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --lanczos-mode 1 --compare-no-bf-lanczos)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_BF_LANCZOS_STATE_CHECK=1;MVMC_LANCZOS_ORACLE_DUMP=1;MVMC_BF_TEST_WORK_SUFFIX=_${test_name}")
endfunction(add_python_vmc_test_bf_lanczos_mpi)

function(add_python_vmc_test_bf_lanczos_nonidentity test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --lanczos-mode 1 --use-nonidentity-init --lanczos-samples 128)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_BF_LANCZOS_STATE_CHECK=1;MVMC_BF_TEST_WORK_SUFFIX=_${test_name}")
endfunction(add_python_vmc_test_bf_lanczos_nonidentity)

function(add_python_vmc_test_bf_lanczos_invalid test_name model mode expected_error)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --lanczos-mode ${mode} --expect-error "${expected_error}")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_BF_TEST_WORK_SUFFIX=_${test_name}")
endfunction(add_python_vmc_test_bf_lanczos_invalid)

function(add_python_vmc_test_bf_invalid test_name model expected_error)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-error "${expected_error}")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_bf_invalid)

function(add_python_vmc_test_bf_fsz test_name)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf_fsz.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_fsz)

function(add_python_vmc_test_bf_fsz_mpi test_name)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf_fsz.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2;OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_fsz_mpi)

function(add_python_vmc_test_bf_fsz_mpi_procs test_name mpi_procs)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf_fsz.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=${mpi_procs};OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_fsz_mpi_procs)

function(add_python_vmc_test_bf_fsz_omp test_name omp_threads)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf_fsz.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads};OMP_NESTED=true;OMP_MAX_ACTIVE_LEVELS=2")
endfunction(add_python_vmc_test_bf_fsz_omp)

function(add_python_vmc_test_bf_fsz_omp_nested_off test_name omp_threads)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf_fsz.py ${test_name})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads};OMP_NESTED=false;OMP_MAX_ACTIVE_LEVELS=1")
endfunction(add_python_vmc_test_bf_fsz_omp_nested_off)

function(add_python_vmc_test_bf_compact test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --compact-backflow)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_bf_compact)

function(add_python_vmc_test_bf_nqp test_name model nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_bf_nqp)

function(add_python_vmc_test_bf_mpi_nqp test_name model nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2;OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_mpi_nqp)

function(add_python_vmc_test_bf_omp_nqp test_name model omp_threads nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_omp_nqp)

function(add_python_vmc_test_bf_energy_omp_nqp test_name model omp_threads nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-energy)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_energy_omp_nqp)

function(add_python_vmc_test_bf_energy_mpi_nqp test_name model nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-energy)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2;OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_energy_mpi_nqp)

function(add_python_vmc_test_bf_energy_mpi_nqp_nsplit test_name model nqp_full nsplit_size)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-energy --set-nsplit-size ${nsplit_size})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=${nsplit_size};OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_energy_mpi_nqp_nsplit)

function(add_python_vmc_test_bf_energy_clean_omp_nqp test_name model omp_threads nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-energy --reject-output "StcOpt:" --reject-output "Error:" --reject-output "error:" --reject-output "warning:" --reject-output "waring:")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads};MVMC_BF_TEST_WORK_SUFFIX=_${test_name}")
endfunction(add_python_vmc_test_bf_energy_clean_omp_nqp)

function(add_python_vmc_test_bf_gradient_omp_nqp test_name model omp_threads nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-gradient)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_gradient_omp_nqp)

function(add_python_vmc_test_bf_projbf_fd_omp_nqp test_name model omp_threads nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-proj-bf-finite-diff)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_projbf_fd_omp_nqp)

function(add_python_vmc_test_bf_twobodyg_omp_nqp test_name model omp_threads nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-twobodyg)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_twobodyg_omp_nqp)

function(add_python_vmc_test_bf_twobodyg_mpi test_name model nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-twobodyg)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2")
endfunction(add_python_vmc_test_bf_twobodyg_mpi)

function(add_python_vmc_test_bf_twobodygex_omp_nqp test_name model omp_threads nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-twobodygex)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_twobodygex_omp_nqp)

function(add_python_vmc_test_bf_twobodygex_mpi test_name model nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-twobodygex)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2")
endfunction(add_python_vmc_test_bf_twobodygex_mpi)

function(add_python_vmc_test_bf_green2_bruteforce test_name model omp_threads all_complex_flag)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --check-bf-green2-bruteforce --use-nonidentity-init --expect-all-complex-flag ${all_complex_flag})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_green2_bruteforce)

function(add_python_vmc_test_bf_nbody_components test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --check-bf-nbody-components --use-nonidentity-init)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_nbody_components)

function(add_python_vmc_test_bf_nbody_components_mpi_nsplit test_name model mpi_procs nsplit_size)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --check-bf-nbody-components --use-nonidentity-init --set-nsplit-size ${nsplit_size})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_MPI_PROCS=${mpi_procs}")
endfunction(add_python_vmc_test_bf_nbody_components_mpi_nsplit)

function(add_python_vmc_test_bf_nbody_dispatch test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --check-bf-nbody-dispatch --use-nonidentity-init)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_nbody_dispatch)

function(add_python_vmc_test_bf_nbody_state test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --check-bf-nbody-state --use-nonidentity-init)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_nbody_state)

function(add_python_vmc_test_bf_nbody_failure test_name model mode)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --inject-bf-nbody-failure ${mode} --use-nonidentity-init)
    set_tests_properties(${test_name} PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1"
        RESOURCE_LOCK "bf_nbody_failure_${model}_${mode}_workdir")
endfunction(add_python_vmc_test_bf_nbody_failure)

function(add_python_vmc_test_bf_nbody_failure_mpi test_name model mode)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --inject-bf-nbody-failure ${mode} --use-nonidentity-init)
    set_tests_properties(${test_name} PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2;OMP_NUM_THREADS=1"
        RESOURCE_LOCK "bf_nbody_failure_${model}_${mode}_workdir")
endfunction(add_python_vmc_test_bf_nbody_failure_mpi)

function(add_python_vmc_test_bf_real_complex_nonidentity test_name real_model complex_model omp_threads)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${real_model} --compare-real-complex-nonidentity ${complex_model} --reject-output "Error:" --reject-output "error:")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_real_complex_nonidentity)

function(add_python_vmc_test_bf_opt_output_restart test_name model omp_threads)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --check-opt-output-restart --reject-output "Error:" --reject-output "error:")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_opt_output_restart)

function(add_python_vmc_test_bf_nonidentity_ncond test_name model ncond omp_threads)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --use-nonidentity-init --set-ncond ${ncond} --reject-output "Error:" --reject-output "error:")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_nonidentity_ncond)

function(add_python_vmc_test_bf_nonidentity_ncond_mpi test_name model ncond)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --use-nonidentity-init --set-ncond ${ncond} --reject-output "Error:" --reject-output "error:")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2;OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_nonidentity_ncond_mpi)

function(add_python_vmc_test_bf_nonidentity_mpi test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --use-nonidentity-init --reject-output "Error:" --reject-output "error:")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2;OMP_NUM_THREADS=1")
endfunction(add_python_vmc_test_bf_nonidentity_mpi)

function(add_python_vmc_test_bf_mpi test_name model)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2")
endfunction(add_python_vmc_test_bf_mpi)

function(add_python_vmc_test_invalid_bf model expected_error)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-error "${expected_error}")
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_invalid_bf)


function(add_python_uhf_test model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_UHF.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_uhf_test)

# Suppress warnings about PythonInterp package.
# TODO: replace it with find_package(Python3 REQUIRED
#cmake_policy(SET CMP0148 OLD)
find_package(PythonInterp 3.6 REQUIRED)

add_test(NAME UpdateWeight_E2E
    COMMAND ${PYTHON_EXECUTABLE} runtest_update_weight.py)
set_tests_properties(UpdateWeight_E2E PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1"
    RESOURCE_LOCK "update_weight_e2e_workdir")
if(MPI_FOUND)
    add_test(NAME UpdateWeight_E2E_MPI
        COMMAND ${PYTHON_EXECUTABLE} runtest_update_weight.py)
    set_tests_properties(UpdateWeight_E2E_MPI PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_MPI_PROCS=2"
        RESOURCE_LOCK "update_weight_e2e_workdir")
endif()

add_test(NAME AbsoluteKrylovProfile_AuditParity
    COMMAND ${PYTHON_EXECUTABLE} runtest_absolute_krylov_profile.py
            $<TARGET_FILE:absolute_krylov_profile_driver>)
set_tests_properties(AbsoluteKrylovProfile_AuditParity PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME AbsoluteKrylovP3Evidence_Math
    COMMAND ${PYTHON_EXECUTABLE} runtest_absolute_krylov_p3_evidence.py
            ${CMAKE_BINARY_DIR}/test/python/data/absolute_krylov_p3_measurement_policy.json)
set_tests_properties(AbsoluteKrylovP3Evidence_Math PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
add_test(NAME AbsoluteKrylovP3Evidence_Pipeline
    COMMAND ${PYTHON_EXECUTABLE} runtest_absolute_krylov_p3_pipeline.py
            ${CMAKE_BINARY_DIR}/test/python/data/absolute_krylov_p3_measurement_policy.json
            ${CMAKE_SOURCE_DIR})
set_tests_properties(AbsoluteKrylovP3Evidence_Pipeline PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
if(MPI_FOUND)
    add_test(NAME AbsoluteKrylovProfile_AuditParity_Rank2_mpi
        COMMAND ${PYTHON_EXECUTABLE} runtest_absolute_krylov_profile.py
                $<TARGET_FILE:absolute_krylov_profile_driver>)
    set_tests_properties(AbsoluteKrylovProfile_AuditParity_Rank2_mpi PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_MPIEXEC=${MPIEXEC_EXECUTABLE};MVMC_MPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        PROCESSORS 2)
endif()

add_test(NAME AbsoluteKrylovProductionAnchor
    COMMAND ${PYTHON_EXECUTABLE}
            runtest_absolute_krylov_production_anchor.py
            $<TARGET_FILE:absolute_krylov_classic_ed_driver>)
set_tests_properties(AbsoluteKrylovProductionAnchor PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1"
    RESOURCE_LOCK "absolute_krylov_production_anchor")
if(MPI_FOUND)
    add_test(NAME AbsoluteKrylovProductionAnchor_Rank2_mpi
        COMMAND ${PYTHON_EXECUTABLE}
                runtest_absolute_krylov_production_anchor.py
                $<TARGET_FILE:absolute_krylov_classic_ed_driver>)
    set_tests_properties(AbsoluteKrylovProductionAnchor_Rank2_mpi PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_MPIEXEC=${MPIEXEC_EXECUTABLE};MVMC_MPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        PROCESSORS 2
        RESOURCE_LOCK "absolute_krylov_production_anchor_mpi")
endif()

set(python_test_vmc_model
  HeisenbergChain
  HeisenbergChain_cmp
  HeisenbergChain_fsz
  HubbardChain
  HubbardChain_cmp
  HubbardChain_fsz
  HubbardTetragonal
  HubbardTetragonal_MomentumProjection
  KondoChain
  KondoChain_cmp
  KondoChain_fsz
  KondoChain_Stot1_cmp
)

set(python_test_vmc_model_mpi
  HeisenbergChain_mpi
  HeisenbergChain_cmp_mpi
  HeisenbergChain_fsz_mpi
  HubbardChain_mpi
  HubbardChain_cmp_mpi
  HubbardChain_fsz_mpi
  HubbardTetragonal_mpi
  HubbardTetragonal_MomentumProjection_mpi
  KondoChain_mpi
  KondoChain_cmp_mpi
  KondoChain_fsz_mpi
  KondoChain_Stot1_cmp_mpi
)

set(python_test_vmc_model_expert_mpi
  HubbardChain_SpinJastrow_mpi
)

set(python_test_vmc_model_mode1
  SpinChainLanczos
  HubbardChainLanczos
)

set(python_test_vmc_model_expert
  GeneralRBM_cmp
  TJChainL4_Ncond2_J_path4
  TJChainL4_Ncond2_J_path5
  TJChainL4_Ncond4_J_path5
  HeisenbergChain_SpinJastrow
  HeisenbergChain_fsz_SpinJastrow
  HubbardChain_SpinJastrow
  HubbardChain_fsz_SpinJastrow
)

set(python_test_vmc_model_expert_phys
  Twist_SSH_cmp
  Twist_SSH_fsz
)


set(python_test_vmc_model_sr
  HubbardChain_DiagCG
  HubbardChain_cmp_DiagCG
  HubbardChain_fsz_DiagCG
)


set(python_test_uhf_model
  UHF_HubbardSquare
  UHF_HubbardTriangular
  UHF_InterAll_N2
)


foreach(model ${python_test_vmc_model})
    add_python_vmc_test(${model})
endforeach(model)
set_tests_properties(HubbardTetragonal_MomentumProjection PROPERTIES
    FIXTURES_SETUP lanczos2_projection_input)

foreach(model ${python_test_vmc_model_mpi})
    add_python_vmc_test(${model})
endforeach(model)

foreach(model ${python_test_vmc_model_mode1})
    add_python_vmc_test_mode1(${model})
endforeach(model)

foreach(model ${python_test_vmc_model_expert})
    add_python_vmc_test_expert(${model})
endforeach(model)

add_python_vmc_test_reweight(GeneralRBM_cmp_reweight GeneralRBM_cmp)
add_python_vmc_test_reweight(HubbardChain_fsz_SpinJastrow_reweight HubbardChain_fsz_SpinJastrow)

foreach(model ${python_test_vmc_model_expert_phys})
    add_python_vmc_test_expert_phys(${model})
endforeach(model)

add_python_vmc_test_expert_phys_mpi(Twist_SSH_cmp_mpi Twist_SSH_cmp)
add_python_vmc_test_expert_phys_mpi(Twist_SSH_fsz_mpi Twist_SSH_fsz)

foreach(model ${python_test_vmc_model_expert_mpi})
    add_python_vmc_test_expert_mpi(${model})
endforeach(model)

add_python_vmc_test_doublon_only(DoublonOnly_GreenFunc2)
add_test(NAME DoublonOnly_GreenFunc2_Complex
         COMMAND ${PYTHON_EXECUTABLE} runtest_doublon_only.py DoublonOnly_GreenFunc2 --complex-orbital)
set_tests_properties(DoublonOnly_GreenFunc2_Complex PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
add_python_vmc_test_doublon_only_mpi(DoublonOnly_GreenFunc2_mpi DoublonOnly_GreenFunc2)

add_python_vmc_test_lanczos_expert(TJChainL4_Ncond2_J_path4_Lanczos)
add_python_vmc_test_lanczos_expert(TJChainL4_Ncond2_J_path5_Lanczos)
add_python_vmc_test_lanczos_expert(TJChainL4_Ncond2_J_path4_Lanczos2)
add_python_vmc_test_lanczos_expert(DoublonOnly_Lanczos)
add_python_vmc_test_lanczos_sector_ed(LanczosSectorED)
add_test(NAME Lanczos2_Step1_OmittedExplicit_ByteEquivalent
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_step1_compat.py
    HubbardChainLanczos)
set_tests_properties(Lanczos2_Step1_OmittedExplicit_ByteEquivalent PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1"
    RESOURCE_LOCK "lanczos2_step1_compat_workdir")

add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidStepZero HubbardChain_SpinJastrow
    "NLanczosStep must be 1 or 2" NLanczosStep=0)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidStepNegative HubbardChain_SpinJastrow
    "NLanczosStep must be 1 or 2" NLanczosStep=-1)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidStepThree HubbardChain_SpinJastrow
    "NLanczosStep must be 1 or 2" NLanczosStep=3)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidStepFractional HubbardChain_SpinJastrow
    "NLanczosStep must be 1 or 2" NLanczosStep=1.5)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidLanczosMode HubbardChain_SpinJastrow
    "NLanczosStep=2 requires NLanczosMode=1"
    NLanczosStep=2 NLanczosMode=0 NVMCCalMode=1)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidVmcCalMode HubbardChain_SpinJastrow
    "NLanczosStep=2 requires NVMCCalMode=1"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=0)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidUpdatePath HubbardChain_SpinJastrow
    "NLanczosStep=2 supports NExUpdatePath=0"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1 NExUpdatePath=1)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidFSZ HubbardChain_fsz_SpinJastrow
    "does not support orbital-general (FSZ) mode"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1
    NExUpdatePath=0 NMPTrans=1)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidRBM GeneralRBM_cmp
    "does not support RBM"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidReweight HubbardChain_SpinJastrow
    "NLanczosStep=2 does not support reweight=1"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1 reweight=1)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidReweight_mpi HubbardChain_SpinJastrow
    "NLanczosStep=2 does not support reweight=1"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1 reweight=1)
set_tests_properties(Lanczos2_InvalidReweight_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_MPI_PROCS=2")
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos_InvalidSupportMode HubbardChain_SpinJastrow
    "NLanczosSupportMode must be 0 or 1"
    NLanczosSupportMode=2)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos_InvalidSupportMode_mpi HubbardChain_SpinJastrow
    "NLanczosSupportMode must be 0 or 1"
    NLanczosSupportMode=2)
set_tests_properties(Lanczos_InvalidSupportMode_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_MPI_PROCS=2")
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidPairHopping DoublonOnly_Lanczos
    "does not support PairHopping"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1 NExUpdatePath=0)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidExchange TJChainL4_Ncond2_J_path4_Lanczos
    "does not support Exchange"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1 NExUpdatePath=0)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidInterAll UHF_InterAll_N2
    "does not support InterAll"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidNBodyG NBodyG_Compare
    "does not support NBodyG"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidSpinFlipTransfer HubbardChain_SpinJastrow
    "supports only spin-conserving Transfer terms"
    spin_flip_transfer NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1)
add_test(NAME Lanczos2_InvalidBackFlow
    COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Real
    --expect-error "NLanczosStep=2 does not support BackFlow"
    --lanczos-mode 1 --lanczos-step 2 --vmc-cal-mode 1)
set_tests_properties(Lanczos2_InvalidBackFlow PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1")
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_InvalidLanczosMode_mpi HubbardChain_SpinJastrow
    "NLanczosStep=2 requires NLanczosMode=1"
    NLanczosStep=2 NLanczosMode=0 NVMCCalMode=1)
set_tests_properties(Lanczos2_InvalidLanczosMode_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_MPI_PROCS=2")
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_Heisenberg_InvalidLanczosMode HeisenbergChain_SpinJastrow
    "NLanczosStep=2 requires NLanczosMode=1"
    NLanczosStep=2 NLanczosMode=2 NVMCCalMode=1 NMPTrans=1)
add_python_vmc_test_invalid_expert_with_updates(
    Lanczos2_Heisenberg_InvalidMultiQP HeisenbergChain_SpinJastrow
    "pure-spin mode currently requires NQPFull=1"
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1 NMPTrans=2)
add_test(NAME Lanczos2_ContractPositive_TwoSz0_NonFSZ
    COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py
    HubbardChain_SpinJastrow "Finish calculation." allow_success
    NLanczosStep=2 NLanczosMode=1 NVMCCalMode=1)
set_tests_properties(Lanczos2_ContractPositive_TwoSz0_NonFSZ PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME Lanczos2_Real_ED
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_real_ed.py)
set_tests_properties(Lanczos2_Real_ED PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_real_ed_workdir")

add_test(NAME AbsoluteKrylovClassic_ED
    COMMAND ${PYTHON_EXECUTABLE} runtest_absolute_krylov_classic_ed.py
            $<TARGET_FILE:absolute_krylov_classic_ed_driver>)
set_tests_properties(AbsoluteKrylovClassic_ED PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_test(NAME AbsoluteKrylovClassic_ED_omp
    COMMAND ${PYTHON_EXECUTABLE} runtest_absolute_krylov_classic_ed.py
            $<TARGET_FILE:absolute_krylov_classic_ed_driver>)
set_tests_properties(AbsoluteKrylovClassic_ED_omp PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
if(MPI_FOUND)
    add_test(NAME AbsoluteKrylovClassic_ED_Rank2_mpi
        COMMAND ${PYTHON_EXECUTABLE} runtest_absolute_krylov_classic_ed.py
                $<TARGET_FILE:absolute_krylov_classic_ed_driver>)
    set_tests_properties(AbsoluteKrylovClassic_ED_Rank2_mpi PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_MPIEXEC=${MPIEXEC_EXECUTABLE};MVMC_MPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        PROCESSORS 2)
    add_test(NAME AbsoluteKrylovClassic_ED_Rank4_mpi
        COMMAND ${PYTHON_EXECUTABLE} runtest_absolute_krylov_classic_ed.py
                $<TARGET_FILE:absolute_krylov_classic_ed_driver>)
    set_tests_properties(AbsoluteKrylovClassic_ED_Rank4_mpi PROPERTIES
        ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=4;MVMC_MPIEXEC=${MPIEXEC_EXECUTABLE};MVMC_MPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
        PROCESSORS 4)
endif()
add_test(NAME Lanczos2_Real_ED_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_real_ed.py)
set_tests_properties(Lanczos2_Real_ED_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_real_ed_mpi_workdir")
add_test(NAME Lanczos2_Real_ED_Split2_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_real_ed.py)
set_tests_properties(Lanczos2_Real_ED_Split2_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_NSPLIT_SIZE=2;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_real_ed_split2_mpi_workdir")
add_test(NAME Lanczos2_Real_ED_omp
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_real_ed.py)
set_tests_properties(Lanczos2_Real_ED_omp PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_real_ed_workdir")
add_test(NAME Lanczos2_NonfiniteSample_Fails
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_real_ed.py
    "Lanczos2 moment or transformed matrix is not finite")
set_tests_properties(Lanczos2_NonfiniteSample_Fails PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_TEST_NONFINITE_SAMPLE=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_real_ed_workdir")
add_test(NAME Lanczos2_NonfiniteSample_Fails_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_real_ed.py
    "Lanczos2 moment or transformed matrix is not finite")
set_tests_properties(Lanczos2_NonfiniteSample_Fails_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_LANCZOS2_TEST_NONFINITE_SAMPLE=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_real_ed_mpi_workdir")
add_test(NAME Lanczos2_Complex_ED
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_complex_ed.py)
set_tests_properties(Lanczos2_Complex_ED PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_complex_ed_workdir")
add_test(NAME Lanczos2_Complex_ED_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_complex_ed.py)
set_tests_properties(Lanczos2_Complex_ED_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_complex_ed_mpi_workdir")
add_test(NAME Lanczos2_Complex_ED_omp
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_complex_ed.py)
set_tests_properties(Lanczos2_Complex_ED_omp PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_complex_ed_workdir")
add_test(NAME Lanczos2_Heisenberg_Real_ED
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py)
set_tests_properties(Lanczos2_Heisenberg_Real_ED PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_heisenberg_real_ed_workdir")
add_test(NAME Lanczos2_Heisenberg_Complex_ED
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py --complex)
set_tests_properties(Lanczos2_Heisenberg_Complex_ED PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_heisenberg_complex_ed_workdir")
add_test(NAME Lanczos2_Heisenberg_Ising_ED
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py --ising)
set_tests_properties(Lanczos2_Heisenberg_Ising_ED PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_heisenberg_ising_ed_workdir")
add_test(NAME Lanczos2_Heisenberg_Standard_Smoke
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py --standard-smoke)
set_tests_properties(Lanczos2_Heisenberg_Standard_Smoke PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_heisenberg_standard_smoke_workdir")
add_test(NAME Lanczos2_Heisenberg_InvalidLocSpin
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py --invalid-locspin)
set_tests_properties(Lanczos2_Heisenberg_InvalidLocSpin PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1"
    RESOURCE_LOCK "lanczos2_heisenberg_invalid_locspin_workdir")
add_test(NAME Lanczos2_Heisenberg_InvalidLocSpin_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py --invalid-locspin)
set_tests_properties(Lanczos2_Heisenberg_InvalidLocSpin_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2"
    RESOURCE_LOCK "lanczos2_heisenberg_invalid_locspin_mpi_workdir")
add_test(NAME Lanczos2_Heisenberg_Real_ED_omp
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py)
set_tests_properties(Lanczos2_Heisenberg_Real_ED_omp PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_heisenberg_real_ed_workdir")
add_test(NAME Lanczos2_Heisenberg_Complex_ED_omp
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py --complex)
set_tests_properties(Lanczos2_Heisenberg_Complex_ED_omp PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_heisenberg_complex_ed_workdir")
add_test(NAME Lanczos2_Heisenberg_Real_ED_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py)
set_tests_properties(Lanczos2_Heisenberg_Real_ED_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_heisenberg_real_ed_mpi_workdir")
add_test(NAME Lanczos2_Heisenberg_Complex_ED_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py --complex)
set_tests_properties(Lanczos2_Heisenberg_Complex_ED_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_heisenberg_complex_ed_mpi_workdir")
add_test(NAME Lanczos2_Heisenberg_Real_ED_Split2_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py)
set_tests_properties(Lanczos2_Heisenberg_Real_ED_Split2_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_NSPLIT_SIZE=2;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_heisenberg_real_ed_split2_mpi_workdir")
add_test(NAME Lanczos2_Heisenberg_Real_ED_Rank4_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_heisenberg_ed.py)
set_tests_properties(Lanczos2_Heisenberg_Real_ED_Rank4_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=4;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_heisenberg_real_ed_mpi_workdir")
add_test(NAME Lanczos2_Projected_ED
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_projected_ed.py)
set_tests_properties(Lanczos2_Projected_ED PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1;MVMC_LANCZOS2_TEST_CALHCA_GUARD_AUDIT=1"
    RESOURCE_LOCK "lanczos2_projected_ed_workdir")
add_test(NAME Lanczos2_Projected_Complex_ED_Guard
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_projected_ed.py --complex)
set_tests_properties(Lanczos2_Projected_Complex_ED_Guard PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1;MVMC_LANCZOS2_TEST_CALHCA_GUARD_AUDIT=1"
    RESOURCE_LOCK "lanczos2_projected_complex_ed_workdir")
add_test(NAME Lanczos2_Projected_ED_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_projected_ed.py)
set_tests_properties(Lanczos2_Projected_ED_mpi PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_projected_ed_mpi_workdir")
add_test(NAME Lanczos2_Projected_ED_omp
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_projected_ed.py)
set_tests_properties(Lanczos2_Projected_ED_omp PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    RESOURCE_LOCK "lanczos2_projected_ed_workdir")
add_test(NAME Lanczos2_Projection_Smoke
    COMMAND ${PYTHON_EXECUTABLE} runtest_lanczos2_projection_smoke.py)
set_tests_properties(Lanczos2_Projection_Smoke PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;MVMC_LANCZOS2_STATE_CHECK=1"
    FIXTURES_REQUIRED lanczos2_projection_input
    RESOURCE_LOCK "lanczos2_projection_smoke_workdir")
foreach(model ${python_test_vmc_model_sr})
    add_python_vmc_test_sr(${model})
endforeach(model)

add_python_vmc_test_rescale(RescaleSmat_DiagCG HubbardChain_DiagCG)
add_python_vmc_test_rescale(RescaleSmat_GeneralRBM_smoke GeneralRBM_cmp)

# reweight=1 must give consistent results between the direct (dposv) SR
# solver and the SR-CG Store path (w-weighted <O> accumulation).
add_test(NAME HubbardChain_DiagCG_reweight
         COMMAND ${PYTHON_EXECUTABLE} runtest_sr.py HubbardChain_DiagCG --reweight)
set_tests_properties(HubbardChain_DiagCG_reweight PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
add_test(NAME HubbardChain_cmp_DiagCG_reweight
         COMMAND ${PYTHON_EXECUTABLE} runtest_sr.py HubbardChain_cmp_DiagCG --reweight)
set_tests_properties(HubbardChain_cmp_DiagCG_reweight PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")

add_python_vmc_test_invalid_expert(SpinJastrow_InvalidClassRange "class index")
add_python_vmc_test_invalid_expert(SpinJastrow_InvalidSymmetry "inconsistent class for pair")
add_python_vmc_test_invalid_expert(SpinJastrow_NegativeNIdx "non-negative")
add_python_vmc_test_invalid_expert(DoublonOnly_InvalidNeZero "requires Ne > 0")
add_python_vmc_test_invalid_expert(DoublonOnly_InvalidNeFull "requires Ne < Nsite")
add_python_vmc_test_invalid_expert(DoublonOnly_InvalidLocSpin "NLocalSpin > 0")
add_python_vmc_test_invalid_expert(DoublonOnly_InvalidOrbitalGeneral "OrbitalGeneral=1")

add_python_backflow_helper_test(BackFlow_Helper_ChainNN)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_Identity_Complex)
add_python_vmc_test_bf_fsz_mpi(FSZ_Lanczos_Identity_Complex_mpi)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_Identity_Real)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_SpinChanging_Complex)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_CrossFamilyProbe)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_NonfiniteThreshold)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_NonfiniteWarning)
add_python_vmc_test_bf_fsz_mpi(FSZ_Lanczos_NonfiniteAggregateWarning_mpi)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_RebuildFailureWarning)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_InvalidNegativeMode)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_InvalidMode2)
add_python_vmc_test_bf_fsz_mpi(FSZ_Lanczos_InvalidMode2_mpi)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_InvalidOptimization)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_InvalidSPGaussLeg)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_InvalidPairHop)
add_python_vmc_test_bf_fsz(FSZ_Lanczos_InvalidReweight)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_Identity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_Lanczos_Identity_Complex_mpi)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_Lanczos_Identity_Complex_omp 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_Identity_Real)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_Identity_Real_Experimental)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_SpinChanging_Identity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_Lanczos_SpinChanging_Identity_Complex_mpi)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_Lanczos_SpinChanging_Identity_Complex_omp 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_SpinChanging_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_Lanczos_SpinChanging_NonIdentity_Complex_omp 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_MomentumProjection_Identity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_MomentumProjection_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_NonfiniteThreshold)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_Lanczos_NonfiniteThreshold_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_NonfiniteWarning)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_Lanczos_NonfiniteAggregateWarning_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_RebuildFailureWarning)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_ProfileIsolation)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_InvalidMode2)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_InvalidOptimization)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_InvalidSPGaussLeg)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_InvalidPairHop)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_InvalidExchange)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Lanczos_InvalidInterAll)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_IdentityEnergy_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_IdentityEnergy_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_NBody_Dispatch_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SpinChanging_C1_Identity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_SpinChanging_C1_Identity_Complex_mpi)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_SpinChanging_C1_Identity_Complex_omp 4)
add_python_vmc_test_bf_fsz_omp_nested_off(BackFlow_FSZ_SpinChanging_C1_Identity_Complex_omp_nested_off 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SpinChanging_C1_MomentumProjection_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SpinChanging_C1_MomentumProjection_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_SpinChanging_C1_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SpinChanging_C1_Reject_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_SpinChanging_C1_RankFallback_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SpinChanging_C2_Identity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_SpinChanging_C2_Identity_Complex_mpi)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_SpinChanging_C2_Identity_Complex_omp 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SpinChanging_C2_MomentumProjection_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SpinChanging_C2_MomentumProjection_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SpinChanging_C2_Order_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_C2_Detail_Profile_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_C2_Detail_Profile_Complex_mpi)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_C2_Detail_Profile_Complex_omp 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_C2_MatrixFree_Paths_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_C2_MatrixFree_Paths_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_C2_MatrixFree_Paths_NonIdentity_Complex_omp 4)
add_python_vmc_test_bf_fsz(BackFlow_InvalidTwoBodyHamiltonian_NonFSZ)
add_python_vmc_test_bf_fsz(BackFlow_InvalidReweight)
add_python_vmc_test_bf_fsz_mpi(BackFlow_InvalidReweight_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvalidSpinChangingTransfer)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvalidSpinChangingOneBodyG)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvalidSpinChangingTwoBodyG)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvalidSpinChangingTwoBodyGEx)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvalidReweight)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_InvalidReweight_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_TwoBodyG_StaleBase_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_AffectedRows_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_AffectedRows_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_PfUpdate_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_PfUpdate_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz_mpi_procs(BackFlow_FSZ_PfUpdate_MoreRanks_NonIdentity_Complex_mpi 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_PfUpdate_Fallback_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_PfUpdate_ExactZeroFallback_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_PfUpdate_LapackFallback_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_PfUpdate_RankFallback_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_PfUpdate_DirectFull_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_PfUpdate_PermutedLabels_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_PfUpdate_ExplicitState_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_PfUpdate_InvalidArguments_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvUpdate_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_InvUpdate_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvDetailProfile_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_InvDetailProfile_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvResidualGemm_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_InvResidualGemm_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvWoodburyGemm_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_InvWoodburyGemm_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvScanFusion_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_InvScanFusion_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz_mpi_procs(BackFlow_FSZ_InvUpdate_MoreRanks_NonIdentity_Complex_mpi 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvUpdate_Fallback_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvUpdate_GetrfFallback_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvUpdate_GetriFallback_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_InvUpdate_RankFallback_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvUpdate_DirectFull_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvUpdate_ExplicitState_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_InvUpdate_InvalidArguments_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_MatrixFreeRows_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_MatrixFreeRows_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_MatrixFreeRows_NonIdentity_Complex_omp 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_MatrixFreeRows_InvalidArguments_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_MultiMoveRows_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_MultiMoveRows_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_MultiMoveRows_NonIdentity_Complex_omp 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_MultiMoveRows_InvalidArguments_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_GreenMatrixFree_Lazy_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_GreenMatrixFree_Lazy_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_GreenMatrixFree_Lazy_NonIdentity_Complex_omp 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SamplingMatrixFree_Lazy_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_SamplingMatrixFree_Lazy_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz_omp(BackFlow_FSZ_SamplingMatrixFree_Lazy_NonIdentity_Complex_omp 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SamplingMatrixFree_Reject_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_SamplingMatrixFree_Reject_NonIdentity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_MomentumProjection_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_MomentumProjection_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_MomentumProjection_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SRDiff_Identity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_SRDiff_Identity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SRDiff_MomentumProjection_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SRDiff_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SRDiff_MomentumProjection_NonIdentity_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SRDiff_SpinChanging_NonIdentity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_InvalidSpinChangingOneBodyG_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Optimization_Identity_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_Optimization_Identity_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Optimization_Identity_Real)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Optimization_MomentumProjection_Complex)
add_python_vmc_test_bf_fsz_mpi(BackFlow_FSZ_Optimization_MomentumProjection_Complex_mpi)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Optimization_BackFlowParams_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Optimization_StoredO_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Optimization_SRCG_Complex)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_Optimization_AllFixed_Invalid)
add_python_vmc_test_bf_fsz_mpi_procs(BackFlow_FSZ_Optimization_NSplitSize_mpi 4)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_SpinChanging_C_Optimization_NonIdentity_Complex)
add_python_vmc_test_bf(BackFlow_Identity_Complex BackFlow_Identity_Complex)
add_python_vmc_test_bf_lanczos(BackFlow_Lanczos_Identity_Complex BackFlow_Identity_Complex)
add_test(NAME BackFlow_Lanczos_Identity_Real
         COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Real --lanczos-mode 1 --expect-error "power-Lanczos support mismatch")
set_tests_properties(BackFlow_Lanczos_Identity_Real PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_BF_LANCZOS_STATE_CHECK=1;MVMC_BF_TEST_WORK_SUFFIX=_BackFlow_Lanczos_Identity_Real")
add_test(NAME BackFlow_Lanczos_Identity_Real_Experimental
         COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Real --lanczos-mode 1 --lanczos-support-mode 1 --compare-no-bf-lanczos)
set_tests_properties(BackFlow_Lanczos_Identity_Real_Experimental PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_BF_LANCZOS_STATE_CHECK=1;MVMC_LANCZOS_ORACLE_DUMP=1;MVMC_BF_TEST_WORK_SUFFIX=_BackFlow_Lanczos_Identity_Real_Experimental")
add_python_vmc_test_bf_lanczos_mpi(BackFlow_Lanczos_Identity_Complex_mpi BackFlow_Identity_Complex)
add_test(NAME BackFlow_Lanczos_Identity_Complex_omp
         COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Complex --lanczos-mode 1 --compare-no-bf-lanczos)
set_tests_properties(BackFlow_Lanczos_Identity_Complex_omp PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;MVMC_BF_LANCZOS_STATE_CHECK=1;MVMC_LANCZOS_ORACLE_DUMP=1;MVMC_BF_TEST_WORK_SUFFIX=_BackFlow_Lanczos_Identity_Complex_omp")
add_python_vmc_test_bf_lanczos_nonidentity(BackFlow_Lanczos_NonIdentity_Complex BackFlow_Identity_Complex)
add_test(NAME BackFlow_Lanczos_NonIdentity_Complex_omp
         COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Complex --lanczos-mode 1 --use-nonidentity-init --lanczos-samples 128)
set_tests_properties(BackFlow_Lanczos_NonIdentity_Complex_omp PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;MVMC_BF_LANCZOS_STATE_CHECK=1;MVMC_BF_TEST_WORK_SUFFIX=_BackFlow_Lanczos_NonIdentity_Complex_omp")
add_python_vmc_test_bf_lanczos_nonidentity(BackFlow_Lanczos_NonIdentity_Real BackFlow_Identity_Real)
add_python_vmc_test_bf_lanczos_invalid(BackFlow_Lanczos_InvalidMode2 BackFlow_Identity_Complex 2 "supports only NLanczosMode==1")
add_test(NAME BackFlow_Lanczos_InvalidOptimization
         COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Complex --lanczos-mode 1 --vmc-cal-mode 0 --expect-error "requires NVMCCalMode==1")
set_tests_properties(BackFlow_Lanczos_InvalidOptimization PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_BF_TEST_WORK_SUFFIX=_BackFlow_Lanczos_InvalidOptimization")
add_test(NAME BackFlow_Lanczos_NonfiniteThreshold
         COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Complex --lanczos-mode 1 --expect-lanczos-nonfinite)
set_tests_properties(BackFlow_Lanczos_NonfiniteThreshold PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_LANCZOS_TEST_NONFINITE_SAMPLE=0;MVMC_BF_TEST_WORK_SUFFIX=_BackFlow_Lanczos_NonfiniteThreshold")
add_test(NAME BackFlow_Lanczos_NonfiniteWarning
         COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Complex --lanczos-mode 1 --expect-lanczos-warning)
set_tests_properties(BackFlow_Lanczos_NonfiniteWarning PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_LANCZOS_TEST_NONFINITE_SAMPLE=0;MVMC_BF_TEST_WORK_SUFFIX=_BackFlow_Lanczos_NonfiniteWarning")
add_test(NAME BackFlow_Lanczos_NonfiniteAggregateWarning_mpi
         COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Complex --lanczos-mode 1 --lanczos-samples 75 --expect-lanczos-warning-count 1/150)
set_tests_properties(BackFlow_Lanczos_NonfiniteAggregateWarning_mpi PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_MPI_PROCS=2;MVMC_LANCZOS_TEST_NONFINITE_SAMPLE=0;MVMC_LANCZOS_TEST_NONFINITE_PARENT_RANK=0;MVMC_BF_LANCZOS_STATE_CHECK=1;MVMC_BF_TEST_WORK_SUFFIX=_BackFlow_Lanczos_NonfiniteAggregateWarning_mpi")
add_test(NAME BackFlow_Lanczos_RebuildFailureWarning
         COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Complex --lanczos-mode 1 --expect-lanczos-warning)
set_tests_properties(BackFlow_Lanczos_RebuildFailureWarning PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1;MVMC_LANCZOS_TEST_REBUILD_FAILURES=1;MVMC_BF_LANCZOS_STATE_CHECK=1;MVMC_BF_TEST_WORK_SUFFIX=_BackFlow_Lanczos_RebuildFailureWarning")
add_python_vmc_test_bf_compact(BackFlow_Compact_Identity_Complex BackFlow_Identity_Complex)
add_python_vmc_test_bf_invalid(BackFlow_Identity_Complex_MultiQP BackFlow_Identity_Complex_MultiQP "momentum projection in non-FSZ mode yet")
add_python_vmc_test_bf_invalid(BackFlow_Identity_Real_MultiQP BackFlow_Identity_Real_MultiQP "momentum projection in non-FSZ mode yet")
add_python_vmc_test_bf_mpi(BackFlow_Identity_Complex_mpi BackFlow_Identity_Complex)
add_python_vmc_test_bf_energy_omp_nqp(BackFlow_IdentityEnergy_Real BackFlow_Identity_Real 1 1)
add_python_vmc_test_bf_invalid(BackFlow_Identity_Real_MultiQP_L8s1 BackFlow_Identity_Real_MultiQP_L8s1 "momentum projection in non-FSZ mode yet")
add_python_vmc_test_bf_invalid(BackFlow_Identity_Real_MultiQP_L8s4 BackFlow_Identity_Real_MultiQP_L8s4 "momentum projection in non-FSZ mode yet")
add_python_vmc_test_bf(BackFlow_Optimization_Complex BackFlow_Optimization_Complex)
add_python_vmc_test_bf_compact(BackFlow_Compact_Optimization_Complex BackFlow_Optimization_Complex)
add_python_vmc_test_bf_compact(BackFlow_Compact_Optimization_Real BackFlow_Optimization_Real)
add_python_vmc_test_bf_energy_clean_omp_nqp(BackFlow_Optimization_Real BackFlow_Optimization_Real 1 1)
add_python_vmc_test_bf_invalid(BackFlow_Optimization_Real_MultiQP BackFlow_Optimization_Real_MultiQP "momentum projection in non-FSZ mode yet")
add_python_vmc_test_bf_invalid(BackFlow_Optimization_Complex_MultiQP BackFlow_Optimization_Complex_MultiQP "momentum projection in non-FSZ mode yet")
add_python_vmc_test_bf_twobodyg_omp_nqp(BackFlow_Identity_TwoBodyG_Complex BackFlow_Identity_Complex 1 1)
add_python_vmc_test_bf_twobodyg_mpi(BackFlow_Identity_TwoBodyG_Complex_mpi BackFlow_Identity_Complex 1)
add_python_vmc_test_bf_twobodygex_omp_nqp(BackFlow_Identity_TwoBodyGEx_Complex BackFlow_Identity_Complex 1 1)
add_python_vmc_test_bf_twobodygex_mpi(BackFlow_Identity_TwoBodyGEx_Complex_mpi BackFlow_Identity_Complex 1)
add_python_vmc_test_bf_green2_bruteforce(BackFlow_NonIdentity_TwoBodyG_BruteForce BackFlow_Identity_Complex 1 1)
add_python_vmc_test_bf_green2_bruteforce(BackFlow_NonIdentity_TwoBodyG_BruteForce_Real BackFlow_Identity_Real 1 0)
add_python_vmc_test_bf_nbody_components(BackFlow_NBody_Components_NonFSZ BackFlow_Identity_Complex)
add_python_vmc_test_bf_nbody_components_mpi_nsplit(BackFlow_NBody_Components_NonFSZ_MPI4Split2_mpi BackFlow_Identity_Complex 4 2)
add_python_vmc_test_bf_nbody_dispatch(BackFlow_NBody_Dispatch_NonFSZ BackFlow_Identity_Complex)
add_python_vmc_test_bf_nbody_state(BackFlow_NBody_State_NonFSZ BackFlow_Identity_Complex)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_Workspace BackFlow_Identity_Complex workspace)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_Candidate BackFlow_Identity_Complex candidate)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_Pfaffian BackFlow_Identity_Complex pfaffian)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_InvalidStage BackFlow_Identity_Complex invalid-stage)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_InvalidTerm BackFlow_Identity_Complex invalid-term)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_InvalidStateCheck BackFlow_Identity_Complex invalid-state-check)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_MissingStage BackFlow_Identity_Complex missing-stage)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_MissingTerm BackFlow_Identity_Complex missing-term)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBodyG_Failure_Workspace BackFlow_Identity_Complex nbodyg-workspace)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBodyG_Failure_Candidate BackFlow_Identity_Complex nbodyg-candidate)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBodyG_Failure_Pfaffian BackFlow_Identity_Complex nbodyg-pfaffian)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_InvalidTermNonNumeric BackFlow_Identity_Complex invalid-term-nonnumeric)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_InvalidTermOverflow BackFlow_Identity_Complex invalid-term-overflow)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_InvalidTermPlus BackFlow_Identity_Complex invalid-term-plus)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_InvalidTermSpace BackFlow_Identity_Complex invalid-term-space)
add_python_vmc_test_bf_nbody_failure(BackFlow_NBody_Failure_InvalidTermOutOfRange BackFlow_Identity_Complex invalid-term-out-of-range)
add_python_vmc_test_bf_nbody_failure_mpi(BackFlow_NBody_Failure_InvalidTermNonNumeric_mpi BackFlow_Identity_Complex invalid-term-nonnumeric)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_NBody_State)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_NBody_Failure_Workspace)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_NBody_Failure_Candidate)
add_python_vmc_test_bf_fsz(BackFlow_FSZ_NBody_Failure_Pfaffian)
add_test(
    NAME BackFlow_Complex_CopyIndex_Private_Source
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_backflow_omp_private.py
            ${CMAKE_SOURCE_DIR}/src/mVMC/backflow_measure.c)
add_test(
    NAME BackFlow_Complex_CopyIndex_OMP4
    COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Identity_Complex
            --expect-nqp-full 1
            --compare-no-bf-twobodyg
            --reject-output "Error:"
            --reject-output "error:")
set_tests_properties(
    BackFlow_Complex_CopyIndex_OMP4
    PROPERTIES ENVIRONMENT
               "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=4;MVMC_BF_TEST_WORK_SUFFIX=_BackFlow_Complex_CopyIndex_OMP4")
add_python_vmc_test_bf_real_complex_nonidentity(BackFlow_NonIdentity_RealVsComplex BackFlow_Identity_Real BackFlow_Identity_Complex 1)
add_python_vmc_test_bf_opt_output_restart(BackFlow_Optimization_OptOutputRestart BackFlow_Optimization_Complex 1)
add_test(
    NAME BackFlow_Optimization_SingleSampleOptOutputRestart
    COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Optimization_Complex
            --check-opt-output-restart
            --opt-output-samples 1
            --reject-output "Error:"
            --reject-output "error:")
set_tests_properties(
    BackFlow_Optimization_SingleSampleOptOutputRestart
    PROPERTIES ENVIRONMENT
               "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1")
add_test(
    NAME BackFlow_Optimization_SingleSampleOptOutputRestart_mpi
    COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Optimization_Complex
            --check-opt-output-restart
            --opt-output-samples 1
            --reject-output "Error:"
            --reject-output "error:")
set_tests_properties(
    BackFlow_Optimization_SingleSampleOptOutputRestart_mpi
    PROPERTIES ENVIRONMENT
               "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS=2;OMP_NUM_THREADS=1")
add_test(
    NAME Optimization_SingleSampleOptOutputRestart_NoBackFlow
    COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py BackFlow_Optimization_Complex
            --check-opt-output-restart
            --opt-output-samples 1
            --opt-output-no-backflow
            --reject-output "Error:"
            --reject-output "error:")
set_tests_properties(
    Optimization_SingleSampleOptOutputRestart_NoBackFlow
    PROPERTIES ENVIRONMENT
               "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=1")
add_python_vmc_test_bf_nonidentity_ncond(BackFlow_NonIdentity_HighFilling_Complex BackFlow_Identity_Complex 6 1)
add_python_vmc_test_bf_nonidentity_ncond(BackFlow_NonIdentity_HighFilling_Real BackFlow_Identity_Real 6 1)
add_python_vmc_test_bf_nonidentity_ncond_mpi(BackFlow_NonIdentity_HighFilling_Complex_mpi BackFlow_Identity_Complex 6)
add_python_vmc_test_bf_nonidentity_ncond_mpi(BackFlow_NonIdentity_HighFilling_Real_mpi BackFlow_Identity_Real 6)
add_python_vmc_test_invalid_expert(DoublonOnly_InvalidRBM "not compatible with RBM")
add_python_vmc_test_invalid_expert(BackFlow_InvalidSpinProjection "NSPGaussLeg==1")
add_python_vmc_test_invalid_expert(TJChainL4_Ncond2_J_path4_InvalidRBM "does not support RBM")
add_python_vmc_test_invalid_expert(TJChainL4_Ncond2_J_path5_InvalidRBM "does not support RBM")
add_python_vmc_test_invalid_expert(TJChainL4_Ncond4_J_path4_InvalidFullFilling "requires 2*Ne < Nsite")
add_test(NAME NSplitSize_Zero_Invalid
         COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py HubbardChain_SpinJastrow
                 "NSplitSize (= 0) must be a positive integer" NSplitSize=0 NStore=0)
set_tests_properties(NSplitSize_Zero_Invalid PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
add_test(NAME NSplitSize_Negative_Invalid
         COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py HubbardChain_SpinJastrow
                 "NSplitSize (= -1) must be a positive integer" NSplitSize=-1 NStore=0)
set_tests_properties(NSplitSize_Negative_Invalid PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
add_test(NAME NSplitSize_NStore_Compare_mpi
         COMMAND ${PYTHON_EXECUTABLE} runtest_nstore_nsplit.py HubbardChain_SpinJastrow)
set_tests_properties(NSplitSize_NStore_Compare_mpi PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS_REF=2;MVMC_MPI_PROCS_SPLIT=4;OMP_NUM_THREADS=1")
add_test(NAME NSplitSize_NStore_fsz_Compare_mpi
         COMMAND ${PYTHON_EXECUTABLE} runtest_nstore_nsplit.py HubbardChain_fsz_SpinJastrow)
set_tests_properties(NSplitSize_NStore_fsz_Compare_mpi PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;MVMC_MPI_PROCS_REF=2;MVMC_MPI_PROCS_SPLIT=4;OMP_NUM_THREADS=1")
add_test(NAME NSplitSize_SRCG_Invalid
         COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py HubbardChain_SpinJastrow
                 "NSplitSize > 1 cannot be combined" NSplitSize=2 NStore=0 NSRCG=1)
set_tests_properties(NSplitSize_SRCG_Invalid PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
add_test(NAME NSROptItrStep_Zero_Invalid
         COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py HubbardChain_SpinJastrow
                 "NSROptItrStep (= 0) must be a positive integer" NSROptItrStep=0 NSROptItrSmp=1)
set_tests_properties(NSROptItrStep_Zero_Invalid PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
add_test(NAME NSROptItrSmp_Zero_Invalid
         COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py HubbardChain_SpinJastrow
                 "NSROptItrSmp (= 0) must satisfy" NSROptItrStep=2 NSROptItrSmp=0)
set_tests_properties(NSROptItrSmp_Zero_Invalid PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
add_test(NAME NSROptItrSmp_TooLarge_Invalid
         COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py HubbardChain_SpinJastrow
                 "NSROptItrSmp (= 3) must satisfy" NSROptItrStep=2 NSROptItrSmp=3)
set_tests_properties(NSROptItrSmp_TooLarge_Invalid PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
add_python_vmc_test_nbodyg(NBodyG_Compare)
add_python_vmc_test_nbodyg_mpi(NBodyG_Compare_mpi NBodyG_Compare)
add_python_vmc_test_nbodyg(NBodyG_fsz_Compare)
add_python_vmc_test_nbodyg(NBodyG_fsz_AP_Compare)
add_python_vmc_test_nbodyg(NBodyG_ZeroTerms)
add_python_nbodyg_exact_oracle(NBodyG_Complex_N3_Oracle)
add_python_nbodyg_exact_oracle(NBodyG_fsz_N3_Oracle)
add_python_nbodyg_generated_check(NBodyG_Projector_N3_MultiTerm 2)
add_python_nbodyg_generated_check(NBodyG_MultiBin 1)
add_python_nbodyg_generated_check(BackFlow_NBodyG_NonIdentity_N3_Output 1)
add_python_nbodyg_generated_check(BackFlow_FSZ_NBodyG_SpinChanging_N3_Output 1)
add_python_nbodyg_generated_check(BackFlow_FSZ_NBodyG_Mixed_N4_Output 1)
add_python_nbodyg_config_oracle(BackFlow_NBodyG_NonIdentity_N3_ConfigOracle)
add_python_nbodyg_config_oracle(BackFlow_NBodyG_EffectiveN2_IndependentExactOracle)
add_python_nbodyg_config_oracle(BackFlow_NBodyG_NonIdentity_N4_ConfigOracle)
add_python_nbodyg_config_oracle(BackFlow_FSZ_NBodyG_SpinChanging_N3_ConfigOracle)
add_python_nbodyg_config_oracle(BackFlow_FSZ_NBodyG_EffectiveN2_IndependentExactOracle)
add_python_nbodyg_config_oracle(BackFlow_FSZ_NBodyG_SpinChanging_N4_ConfigOracle)
add_python_nbodyg_config_oracle(BackFlow_FSZ_NBodyG_SpinChanging_N5_ConfigOracle)
add_python_nbodyg_config_oracle(BackFlow_FSZ_NBody_MultiQP_ConfigOracle)
add_python_nbodyg_config_oracle_omp(BackFlow_NBodyG_NonIdentity_N3_ConfigOracle_omp BackFlow_NBodyG_NonIdentity_N3_ConfigOracle)
add_python_nbodyg_config_oracle_omp(BackFlow_FSZ_NBodyG_SpinChanging_N3_ConfigOracle_omp BackFlow_FSZ_NBodyG_SpinChanging_N3_ConfigOracle)
add_python_nbodyg_config_oracle_mpi(BackFlow_NBodyG_NonIdentity_N3_ConfigOracle_mpi BackFlow_NBodyG_NonIdentity_N3_ConfigOracle)
add_python_nbodyg_config_oracle_mpi(BackFlow_FSZ_NBodyG_SpinChanging_N3_ConfigOracle_mpi BackFlow_FSZ_NBodyG_SpinChanging_N3_ConfigOracle)
add_python_vmc_test_invalid_expert(NBodyG_InvalidNegativeCount "NNBodyG")
add_python_vmc_test_invalid_expert(NBodyG_InvalidZeroOrder "term order")
add_python_vmc_test_invalid_expert(NBodyG_InvalidFieldCount "malformed NBodyG term")
add_python_vmc_test_invalid_expert(NBodyG_InvalidExtraToken "extra token")
add_python_vmc_test_invalid_expert(NBodyG_InvalidSite "site index")
add_python_vmc_test_invalid_expert(NBodyG_InvalidSpin "spin index")
add_python_vmc_test_invalid_expert(NBodyG_InvalidSpinChangeNonFsz "spin-changing")
add_python_nbodyg_invalid(NBodyG_InvalidBackFlow "BackFlow NBodyG requires complex variational parameters")
add_python_nbodyg_invalid_mpi(NBodyG_InvalidBackFlow_mpi "BackFlow NBodyG requires complex variational parameters")
add_python_nbodyg_invalid(NBodyG_InvalidBackFlowLanczos "BackFlow NBodyG is not supported with Lanczos mode")
add_python_nbodyg_invalid(NBodyG_InvalidBackFlowSpinChangeNonFSZ "spin-changing NBodyG factor requires orbital-general mode")
add_python_nbodyg_invalid(NBodyG_InvalidBackFlowFszSpinChangeFixedSz "BackFlow FSZ spin-changing NBodyG factors require TwoSz==-1")
add_python_nbodyinterall_exact_oracle(NBodyInterAll_N1_TransferLike)
add_python_nbodyinterall_exact_oracle(NBodyInterAll_N2_InterAllLike)
add_python_nbodyinterall_exact_oracle(NBodyInterAll_N3_Oracle)
add_python_nbodyinterall_exact_oracle(NBodyInterAll_fsz_N3_Oracle)
add_python_nbodyinterall_exact_oracle(NBodyInterAll_ZeroTerms)
add_python_nbodyinterall_generated_check(NBodyInterAll_Projector_N3_MultiTerm 2)
add_python_nbodyinterall_generated_check(BackFlow_NBodyInterAll_NonIdentity_N1_Energy 1)
add_python_nbodyinterall_generated_check(BackFlow_NBodyInterAll_NonIdentity_N2_Energy 1)
add_python_nbodyinterall_generated_check(BackFlow_NBodyInterAll_NonIdentity_Mixed_Energy 1)
add_python_nbodyinterall_generated_check(BackFlow_NBody_ZeroTerms_NonFSZ 1)
add_python_nbodyinterall_generated_check(BackFlow_FSZ_NBodyInterAll_SpinChanging_Mixed_Energy 1)
add_python_nbodyinterall_generated_check(BackFlow_FSZ_NBody_ZeroTerms 1)
add_python_nbodyinterall_config_oracle(BackFlow_NBodyInterAll_NonIdentity_Mixed_ConfigOracle)
add_python_nbodyinterall_config_oracle(BackFlow_FSZ_NBodyInterAll_Mixed_ConfigOracle)
add_python_nbodyinterall_config_oracle_omp(BackFlow_NBodyInterAll_NonIdentity_Mixed_ConfigOracle_omp BackFlow_NBodyInterAll_NonIdentity_Mixed_ConfigOracle)
add_python_nbodyinterall_config_oracle_omp(BackFlow_FSZ_NBodyInterAll_Mixed_ConfigOracle_omp BackFlow_FSZ_NBodyInterAll_Mixed_ConfigOracle)
add_python_nbodyinterall_config_oracle_mpi(BackFlow_NBodyInterAll_NonIdentity_Mixed_ConfigOracle_mpi BackFlow_NBodyInterAll_NonIdentity_Mixed_ConfigOracle)
add_python_nbodyinterall_config_oracle_mpi(BackFlow_FSZ_NBodyInterAll_Mixed_ConfigOracle_mpi BackFlow_FSZ_NBodyInterAll_Mixed_ConfigOracle)
add_python_nbodyinterall_mpi(NBodyInterAll_N3_Oracle_mpi NBodyInterAll_N3_Oracle)
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidNegativeCount "NNBodyInterAll")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidZeroOrder "term order")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidFieldCount "invalid NBodyInterAll coefficient")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidExtraToken "extra token")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidSite "site index")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidSpin "spin index")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidCoeff "invalid NBodyInterAll coefficient")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidSpinChangeNonFsz "spin-changing")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidRealPath "requires complex variational parameters")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidBackFlow "BackFlow NBodyInterAll is not supported with Lanczos mode")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidBackFlowReal "BackFlow NBodyInterAll requires complex variational parameters")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidBackFlowSpinChangeNonFSZ "spin-changing NBodyInterAll factor requires orbital-general mode")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidBackFlowFszSpinChangeFixedSz "BackFlow FSZ spin-changing NBodyInterAll factors require TwoSz==-1")
add_python_nbodyinterall_invalid(NBodyInterAll_InvalidLanczos "not implemented for Lanczos mode")
add_test(NAME SpinJastrow_InvalidInIndex
         COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py SpinJastrow_InvalidInIndex "duplicated index")
set_tests_properties(SpinJastrow_InvalidInIndex PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")

foreach(model ${python_test_uhf_model})
    add_python_uhf_test(${model})
endforeach(model)

add_test(NAME UHF_InterAll COMMAND ${PYTHON_EXECUTABLE} test_UHF_InterAll.py)
set_tests_properties(UHF_InterAll PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")

add_subdirectory(data)
