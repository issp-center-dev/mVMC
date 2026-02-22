# for GCC Compiler (Homebrew)
# Disable USE_GEMMT to build without BLIS (uses reference dskr2k/zskr2k instead)
set(USE_GEMMT OFF CACHE BOOL "Use GEMMT (requires BLIS)" FORCE)

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

list(SORT _gcc_bins)
list(SORT _gxx_bins)
list(SORT _gfc_bins)
list(GET _gcc_bins -1 GCC_C_COMPILER)
list(GET _gxx_bins -1 GCC_CXX_COMPILER)
list(GET _gfc_bins -1 GCC_Fortran_COMPILER)
set(CMAKE_C_COMPILER "${GCC_C_COMPILER}" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "${GCC_CXX_COMPILER}" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "${GCC_Fortran_COMPILER}" CACHE STRING "" FORCE)
message(STATUS "Using Homebrew GCC: ${GCC_C_COMPILER}")

set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall  -Wformat -Werror=format-security")
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wno-unknown-pragmas ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unknown-pragmas ")
