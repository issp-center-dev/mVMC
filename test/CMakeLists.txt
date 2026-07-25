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

add_subdirectory(python)
