file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_mode1.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_mpi.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_expert.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_invalid_expert.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_sr.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_UHF.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_UHF_InterAll.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_expert_phys.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_expert_mpi.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_doublon_only.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_reweight.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_rescale.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/backflow_def_helper.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_bf.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_lanczos_expert.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_lanczos_sector_ed.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)

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

function(add_python_vmc_test_expert_phys model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_expert_phys.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_expert_phys)

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

function(add_python_vmc_test_bf_nqp test_name model nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_bf_nqp)

function(add_python_vmc_test_bf_omp_nqp test_name model omp_threads nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full})
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_omp_nqp)

function(add_python_vmc_test_bf_energy_omp_nqp test_name model omp_threads nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-energy)
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_energy_omp_nqp)

function(add_python_vmc_test_bf_energy_clean_omp_nqp test_name model omp_threads nqp_full)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --expect-nqp-full ${nqp_full} --compare-no-bf-energy --reject-output "StcOpt:" --reject-output "Error:" --reject-output "error:" --reject-output "warning:" --reject-output "waring:")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
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

function(add_python_vmc_test_bf_real_complex_nonidentity test_name real_model complex_model omp_threads)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${real_model} --compare-real-complex-nonidentity ${complex_model} --reject-output "Error:" --reject-output "error:")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_real_complex_nonidentity)

function(add_python_vmc_test_bf_opt_output_restart test_name model omp_threads)
    add_test(NAME ${test_name} COMMAND ${PYTHON_EXECUTABLE} runtest_bf.py ${model} --check-opt-output-restart --reject-output "Error:" --reject-output "error:")
    set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python;OMP_NUM_THREADS=${omp_threads}")
endfunction(add_python_vmc_test_bf_opt_output_restart)

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
add_python_vmc_test_bf(BackFlow_Identity_Complex BackFlow_Identity_Complex)
add_python_vmc_test_bf_nqp(BackFlow_Identity_Complex_MultiQP BackFlow_Identity_Complex_MultiQP 2)
add_python_vmc_test_bf_mpi(BackFlow_Identity_Complex_mpi BackFlow_Identity_Complex)
add_python_vmc_test_bf_energy_omp_nqp(BackFlow_IdentityEnergy_Real BackFlow_Identity_Real 1 1)
add_python_vmc_test_bf(BackFlow_Optimization_Complex BackFlow_Optimization_Complex)
add_python_vmc_test_bf_energy_clean_omp_nqp(BackFlow_Optimization_Real BackFlow_Optimization_Real 1 1)
add_python_vmc_test_bf_omp_nqp(BackFlow_Optimization_Complex_MultiQP_omp4 BackFlow_Optimization_Complex_MultiQP 4 2)
add_python_vmc_test_bf_energy_omp_nqp(BackFlow_Optimization_IdentityEnergy_MultiQP_omp4 BackFlow_Optimization_Complex_MultiQP 4 2)
add_python_vmc_test_bf_gradient_omp_nqp(BackFlow_Optimization_IdentityGradient_MultiQP_omp4 BackFlow_Optimization_Complex_MultiQP 4 2)
add_python_vmc_test_bf_projbf_fd_omp_nqp(BackFlow_Optimization_ProjBFFiniteDiff_MultiQP_omp4 BackFlow_Optimization_Complex_MultiQP 4 2)
add_python_vmc_test_bf_twobodyg_omp_nqp(BackFlow_Identity_TwoBodyG_Complex BackFlow_Identity_Complex 1 1)
add_python_vmc_test_bf_twobodyg_mpi(BackFlow_Identity_TwoBodyG_Complex_mpi BackFlow_Identity_Complex 1)
add_python_vmc_test_bf_twobodygex_omp_nqp(BackFlow_Identity_TwoBodyGEx_Complex BackFlow_Identity_Complex 1 1)
add_python_vmc_test_bf_twobodygex_mpi(BackFlow_Identity_TwoBodyGEx_Complex_mpi BackFlow_Identity_Complex 1)
add_python_vmc_test_bf_real_complex_nonidentity(BackFlow_NonIdentity_RealVsComplex BackFlow_Identity_Real BackFlow_Identity_Complex 1)
add_python_vmc_test_bf_opt_output_restart(BackFlow_Optimization_OptOutputRestart BackFlow_Optimization_Complex 1)
add_python_vmc_test_invalid_expert(DoublonOnly_InvalidRBM "not compatible with RBM")
add_python_vmc_test_invalid_expert(BackFlow_InvalidSpinProjection "NSPGaussLeg==1")
add_python_vmc_test_invalid_expert(TJChainL4_Ncond2_J_path4_InvalidRBM "does not support RBM")
add_python_vmc_test_invalid_expert(TJChainL4_Ncond2_J_path5_InvalidRBM "does not support RBM")
add_python_vmc_test_invalid_expert(TJChainL4_Ncond4_J_path4_InvalidFullFilling "requires 2*Ne < Nsite")
add_test(NAME SpinJastrow_InvalidInIndex
         COMMAND ${PYTHON_EXECUTABLE} runtest_invalid_expert.py SpinJastrow_InvalidInIndex "duplicated index")
set_tests_properties(SpinJastrow_InvalidInIndex PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")

foreach(model ${python_test_uhf_model})
    add_python_uhf_test(${model})
endforeach(model)

add_test(NAME UHF_InterAll COMMAND ${PYTHON_EXECUTABLE} test_UHF_InterAll.py)
set_tests_properties(UHF_InterAll PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")

add_subdirectory(data)
