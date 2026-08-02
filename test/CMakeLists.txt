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

add_subdirectory(python)
