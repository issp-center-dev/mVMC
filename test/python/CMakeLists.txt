file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtest.py DESTINATION ${CMAKE_BINARY_DIR}/test/python)
function(add_python_test model)
    add_test(NAME ${model} COMMAND ${PYTHON_EXECUTABLE} runtest.py ${model})
    set_tests_properties(${model} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_test)

find_package(PythonInterp 2.7 REQUIRED)

set(python_test_model
  HeisenbergChain
  HubbardTetragonal
  KondoChain
  HeisenbergChain_cmp
  HubbardTetragonal_cmp
  KondoChain_cmp
  KondoChain_Stot1_cmp
)

foreach(model ${python_test_model})
    add_python_test(${model})
endforeach(model)

add_subdirectory(data)
