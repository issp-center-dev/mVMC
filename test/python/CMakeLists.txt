file(COPY ${CMAKE_SOURCE_DIR}/python/ DESTINATION ${CMAKE_BINARY_DIR}/test/python FILES_MATCHING PATTERN "*.py")
function(add_python_test test)
    add_test(NAME python_${test} COMMAND ${PYTHON_EXECUTABLE} ${test}.py)
    set_tests_properties(python_${test} PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/test/python")
endfunction(add_python_test)

find_package(PythonInterp REQUIRED)

file(GLOB_RECURSE TEST_FILES "${CMAKE_CURRENT_SOURCE_DIR}/*.py")
foreach(testfile ${TEST_FILES})
    configure_file(${testfile} ${CMAKE_BINARY_DIR}/test/python COPYONLY)
endforeach()

set(python_test_src
  test_HeisenbergChain
  test_UHF_HubbardSquare
  test_UHF_HubbardTriangular
  test_UHF_InterAll
)

foreach(test ${python_test_src})
    add_python_test(${test})
endforeach(test)

add_subdirectory(data)
