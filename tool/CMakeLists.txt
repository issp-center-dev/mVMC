# include guard
cmake_minimum_required(VERSION 2.8.0 )

add_definitions(-D_mVMC)
if(${CMAKE_PROJECT_NAME} STREQUAL "Project")
  message(FATAL_ERROR "cmake should be executed not for 'src' subdirectory, but for the top directory of mVMC.")
endif(${CMAKE_PROJECT_NAME} STREQUAL "Project")

add_executable(greenr2k greenr2k.F90)
target_link_libraries(greenr2k ${LAPACK_LIBRARIES})

install(TARGETS greenr2k RUNTIME DESTINATION bin)
#
# Scripts
#
configure_file(wout2geom.sh wout2geom.sh COPYONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/wout2geom.sh DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
configure_file(respack2wan90.py respack2wan90.py COPYONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/respack2wan90.py DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
