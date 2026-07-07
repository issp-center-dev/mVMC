# for GCC Compiler (Homebrew)
# Disable USE_GEMMT to build without BLIS (uses reference dskr2k/zskr2k instead)
set(USE_GEMMT OFF CACHE BOOL "Use GEMMT (requires BLIS)" FORCE)
# Prefer Accelerate for the mac_gcc verification build. Homebrew BLIS 2.0 has
# shown incorrect dgemm_ results for TN/N=4,16 shapes used by BackFlow batching.
set(BLA_VENDOR Apple CACHE STRING "BLAS/LAPACK vendor" FORCE)

# Auto-detect Homebrew GCC (avoid Apple Clang's /usr/bin/gcc)
execute_process(
    COMMAND brew --prefix gcc
    OUTPUT_VARIABLE GCC_PREFIX
    RESULT_VARIABLE GCC_PREFIX_RESULT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT GCC_PREFIX_RESULT EQUAL 0 OR GCC_PREFIX STREQUAL "")
    message(FATAL_ERROR "Homebrew GCC not found. Install it with: brew install gcc")
endif()

file(GLOB _gcc_bins "${GCC_PREFIX}/bin/gcc-[0-9]*")
file(GLOB _gxx_bins "${GCC_PREFIX}/bin/g++-[0-9]*")
file(GLOB _gfc_bins "${GCC_PREFIX}/bin/gfortran-[0-9]*")
if(NOT _gcc_bins OR NOT _gxx_bins OR NOT _gfc_bins)
    message(FATAL_ERROR "Homebrew GCC binaries not found in ${GCC_PREFIX}/bin/")
endif()

# Select the highest GCC version numerically (avoids lexicographic pitfall where
# gcc-9 sorts after gcc-14 as a string).
set(GCC_C_COMPILER "")
set(_gcc_best "0")
foreach(_bin IN LISTS _gcc_bins)
    get_filename_component(_name "${_bin}" NAME)
    if(_name MATCHES "^gcc-([0-9]+)$")
        if("${CMAKE_MATCH_1}" VERSION_GREATER "${_gcc_best}")
            set(_gcc_best "${CMAKE_MATCH_1}")
            set(GCC_C_COMPILER "${_bin}")
        endif()
    endif()
endforeach()

set(GCC_CXX_COMPILER "")
set(_gxx_best "0")
foreach(_bin IN LISTS _gxx_bins)
    get_filename_component(_name "${_bin}" NAME)
    if(_name MATCHES "^g\\+\\+-([0-9]+)$")
        if("${CMAKE_MATCH_1}" VERSION_GREATER "${_gxx_best}")
            set(_gxx_best "${CMAKE_MATCH_1}")
            set(GCC_CXX_COMPILER "${_bin}")
        endif()
    endif()
endforeach()

set(GCC_Fortran_COMPILER "")
set(_gfc_best "0")
foreach(_bin IN LISTS _gfc_bins)
    get_filename_component(_name "${_bin}" NAME)
    if(_name MATCHES "^gfortran-([0-9]+)$")
        if("${CMAKE_MATCH_1}" VERSION_GREATER "${_gfc_best}")
            set(_gfc_best "${CMAKE_MATCH_1}")
            set(GCC_Fortran_COMPILER "${_bin}")
        endif()
    endif()
endforeach()

if(GCC_C_COMPILER STREQUAL "" OR GCC_CXX_COMPILER STREQUAL "" OR GCC_Fortran_COMPILER STREQUAL "")
    message(FATAL_ERROR "No versioned GCC binaries matching gcc-N/g++-N/gfortran-N found in ${GCC_PREFIX}/bin/")
endif()
set(CMAKE_C_COMPILER "${GCC_C_COMPILER}" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "${GCC_CXX_COMPILER}" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "${GCC_Fortran_COMPILER}" CACHE STRING "" FORCE)
message(STATUS "Using Homebrew GCC: ${GCC_C_COMPILER}")

set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall  -Wformat -Werror=format-security")
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wno-unknown-pragmas ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unknown-pragmas ")
