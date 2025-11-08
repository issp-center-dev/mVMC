file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_mode1.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_mpi.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_expert.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_sr.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_UHF.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_UHF_InterAll.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest_expert_phys.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)

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

function(add_python_vmc_test_expert_phys model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_expert_phys.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_expert_phys)


function(add_python_vmc_test_sr model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_sr.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_vmc_test_sr)


function(add_python_uhf_test model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest_UHF.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_uhf_test)

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

set(python_test_vmc_model_mode1
  SpinChainLanczos
  HubbardChainLanczos
)

set(python_test_vmc_model_expert
  GeneralRBM_cmp
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

foreach(model ${python_test_vmc_model_expert_phys})
    add_python_vmc_test_expert_phys(${model})
endforeach(model)


foreach(model ${python_test_vmc_model_sr})
    add_python_vmc_test_sr(${model})
endforeach(model)

foreach(model ${python_test_uhf_model})
    add_python_uhf_test(${model})
endforeach(model)

add_test(NAME UHF_InterAll COMMAND ${PYTHON_EXECUTABLE} test_UHF_InterAll.py)
set_tests_properties(UHF_InterAll PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")

add_subdirectory(data)
