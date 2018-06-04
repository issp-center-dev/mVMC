foreach(datafile ${DATA_FILES})
    configure_file(${datafile} ${CMAKE_BINARY_DIR}/test/ COPYONLY)
endforeach()
add_subdirectory(python)
