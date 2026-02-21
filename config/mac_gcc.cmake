# for GCC Compiler (Homebrew)
# BLIS を入れずにビルドするため USE_GEMMT=OFF（GEMMT は参照実装の dskr2k/zskr2k を使用）
set(USE_GEMMT OFF CACHE BOOL "Use GEMMT (requires BLIS)" FORCE)

set(CMAKE_C_COMPILER "/opt/homebrew/bin/gcc-15" CACHE STRING "" FORCE)
set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-15" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "/opt/homebrew/bin/gfortran-15" CACHE STRING "" FORCE)

set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall  -Wformat -Werror=format-security")
set(CMAKE_C_FLAGS_RELEASE "-O3 -Wno-unknown-pragmas ")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wno-unknown-pragmas ")
