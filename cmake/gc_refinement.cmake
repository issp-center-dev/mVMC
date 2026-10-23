# Optional exact-input GC fallback. ON requires dependencies; AUTO preserves
# builds on systems without the development packages and retains fail-closed
# behavior there. Production validation must record whether this is enabled.
set(MVMC_GC_REFINEMENT AUTO CACHE STRING "GC refinement: AUTO, ON, or OFF")
set_property(CACHE MVMC_GC_REFINEMENT PROPERTY STRINGS AUTO ON OFF)
if(NOT MVMC_GC_REFINEMENT MATCHES "^(AUTO|ON|OFF)$")
  message(FATAL_ERROR "MVMC_GC_REFINEMENT must be AUTO, ON, or OFF")
endif()
if(NOT MVMC_GC_REFINEMENT STREQUAL "OFF")
  find_path(MVMC_MPFR_INCLUDE_DIR mpfr.h HINTS
    "$ENV{HOMEBREW_PREFIX}/opt/mpfr/include")
  find_path(MVMC_GMP_INCLUDE_DIR gmpxx.h HINTS
    "$ENV{HOMEBREW_PREFIX}/opt/gmp/include")
  find_library(MVMC_MPFR_LIBRARY mpfr HINTS
    "$ENV{HOMEBREW_PREFIX}/opt/mpfr/lib")
  find_library(MVMC_GMP_LIBRARY gmp HINTS
    "$ENV{HOMEBREW_PREFIX}/opt/gmp/lib")
  if(MVMC_MPFR_INCLUDE_DIR AND MVMC_GMP_INCLUDE_DIR AND
     MVMC_MPFR_LIBRARY AND MVMC_GMP_LIBRARY)
    add_library(mvmc_gc_refinement STATIC src/mVMC/krylov_refinement.cpp)
    # In GNU dialect mode libstdc++'s <complex.h> undefines the C ABI
    # compatibility macro on every include. Strict C++ leaves it intact.
    set_target_properties(mvmc_gc_refinement PROPERTIES CXX_EXTENSIONS OFF)
    target_include_directories(mvmc_gc_refinement PRIVATE
      src/mVMC/include ${MVMC_MPFR_INCLUDE_DIR} ${MVMC_GMP_INCLUDE_DIR})
    target_compile_definitions(mvmc_gc_refinement PRIVATE
      MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE MVMC_ENABLE_POWER_LANCZOS_P5_CORE)
    # Finite-input guards and interval boundaries must survive optimization.
    # Do not rely on compiler defaults (Intel changed fast in 2024.2).
    if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang|FujitsuClang)$")
      target_compile_options(mvmc_gc_refinement PRIVATE
        -Wall -Wextra -fno-fast-math -ffp-contract=off)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
      target_compile_options(mvmc_gc_refinement PRIVATE
        -fp-model=strict -ffp-contract=off)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      target_compile_options(mvmc_gc_refinement PRIVATE -fp-model strict)
    else()
      message(FATAL_ERROR
        "GC refinement has no floating-point contract for C++ compiler "
        "${CMAKE_CXX_COMPILER_ID}. Use a supported compiler or set "
        "MVMC_GC_REFINEMENT=OFF to retain fail-closed behavior.")
    endif()
    target_link_libraries(mvmc_gc_refinement
      ${MVMC_MPFR_LIBRARY} ${MVMC_GMP_LIBRARY})
    message(STATUS "GC exact-input refinement enabled (MPFR/GMP)")
  elseif(MVMC_GC_REFINEMENT STREQUAL "ON")
    message(FATAL_ERROR "GC refinement requires MPFR and GMP development packages")
  else()
    message(STATUS "GC refinement unavailable: unresolved amplitudes remain errors")
  endif()
endif()
function(mvmc_enable_gc_refinement target)
  if(TARGET mvmc_gc_refinement)
    target_compile_definitions(${target} PRIVATE MVMC_HAVE_MPFR_REFINEMENT)
    # The C++ archive uses C symbols supplied by the executable.
    target_link_libraries(${target} mvmc_gc_refinement)
  endif()
endfunction()
