if(NOT DEFINED BLIS_ARTIFACT_CONFIG)
  execute_process(
    COMMAND uname -m
    COMMAND tr -d '\n'
    OUTPUT_VARIABLE ARCHITECTURE)
  if(${ARCHITECTURE} STREQUAL "x86_64")
    execute_process(
      COMMAND cat /proc/cpuinfo
      COMMAND grep Intel
      COMMAND head -1
      COMMAND tr -d '\n'
      OUTPUT_VARIABLE UARCH_INTEL)
    if(NOT "${UARCH_INTEL} " STREQUAL " ")
      set(BLIS_ARTIFACT_CONFIG "intel64")
    else()
      set(BLIS_ARTIFACT_CONFIG "amd64")
    endif()
  elseif(${ARCHITECTURE} STREQUAL "aarch64" OR
      ${ARCHITECTURE} STREQUAL "arm64")
      # Currently SVE requires being set manually.
      set(BLIS_ARTIFACT_CONFIG "arm64")
  else()
    message(FATAL_ERROR "Failed to recognize architecture ${ARCHITECTURE}.")
  endif()
endif(NOT DEFINED BLIS_ARTIFACT_CONFIG)

message(STATUS "Downloading BLIS artifact...")

set(BLIS_ARCHIVE ${CMAKE_CURRENT_BINARY_DIR}/libblis_artifact.tar.gz)
if(${BLIS_ARTIFACT_CONFIG} STREQUAL "intel64")
  set(BLIS_ARTIFACT_URL https://github.com/xrq-phys/blis/releases/download/sv0.8.1+arm/libblis_intel64_gcc.tar.gz)
elseif(${BLIS_ARTIFACT_CONFIG} STREQUAL "amd64")
  set(BLIS_ARTIFACT_URL https://github.com/xrq-phys/blis/releases/download/sv0.8.1+arm/libblis_amd64_gcc.tar.gz)
elseif(${BLIS_ARTIFACT_CONFIG} STREQUAL "arm64")
  set(BLIS_ARTIFACT_URL https://github.com/xrq-phys/blis/releases/download/sv0.8.1+arm/libblis_cortexa57_aarch64-linux-gnu-gcc.tar.gz)
elseif(${BLIS_ARTIFACT_CONFIG} STREQUAL "armsve")
  set(BLIS_ARTIFACT_URL https://github.com/xrq-phys/blis/releases/download/sv0.8.1+arm/libblis_armsve_aarch64-linux-gnu-gcc-10.tar.gz)
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

