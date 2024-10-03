if(NOT DEFINED BLIS_ARTIFACT_CONFIG)
  execute_process(
    COMMAND uname -m
    COMMAND tr -d '\n'
    OUTPUT_VARIABLE ARCHITECTURE)
  if(${ARCHITECTURE} STREQUAL "x86_64")
      set(BLIS_ARTIFACT_CONFIG "x86_64")
  elseif(${ARCHITECTURE} STREQUAL "aarch64" OR ${ARCHITECTURE} STREQUAL "arm64")
      set(BLIS_ARTIFACT_CONFIG "arm64")
  else()
    message(FATAL_ERROR "Failed to recognize architecture ${ARCHITECTURE}.")
  endif()
endif(NOT DEFINED BLIS_ARTIFACT_CONFIG)

message(STATUS "Downloading BLIS artifact...")

set(BLIS_ARCHIVE ${CMAKE_CURRENT_BINARY_DIR}/libblis_artifact.tar.gz)
if(${BLIS_ARTIFACT_CONFIG} STREQUAL "x86_64")
  set(BLIS_ARTIFACT_URL https://github.com/JuliaBinaryWrappers/blis_jll.jl/releases/download/blis-v0.9.0%2B4/blis.v0.9.0.x86_64-linux-gnu.tar.gz)
elseif(${BLIS_ARTIFACT_CONFIG} STREQUAL "arm64" OR ${BLIS_ARTIFACT_CONFIG} STREQUAL "armsve")
  set(BLIS_ARTIFACT_URL https://github.com/JuliaBinaryWrappers/blis_jll.jl/releases/download/blis-v0.9.0%2B4/blis.v0.9.0.aarch64-linux-gnu.tar.gz)
elseif(${BLIS_ARTIFACT_CONFIG} STREQUAL "a64fx")
  set(BLIS_ARTIFACT_URL https://github.com/xrq-phys/blis/releases/download/sv0.8.1+arm/libblis_a64fx_aarch64-linux-gnu-gcc-10.tar.gz)
else()
  message(FATAL_ERROR "Not BLIS artifact available.")
endif()

message("${BLIS_ARTIFACT_URL} ${BLIS_ARCHIVE}")
file(DOWNLOAD "${BLIS_ARTIFACT_URL}" "${BLIS_ARCHIVE}")

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib/libblis.a
    ${CMAKE_CURRENT_BINARY_DIR}/include/blis/blis.h
  COMMAND ${CMAKE_COMMAND} -E tar -zxf ${BLIS_ARCHIVE}
  DEPENDS ${BLIS_ARCHIVE})
add_custom_target(blis_target
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/lib/libblis.a)
add_custom_target(blis_include
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/include/blis/blis.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/blis)

add_library(blis STATIC IMPORTED GLOBAL)
set_target_properties(blis PROPERTIES
  IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lib/libblis.a)
add_dependencies(blis blis_target)

