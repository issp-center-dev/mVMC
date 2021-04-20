#!/bin/bash
if [ -z ${1} ] || [ ${1} = "help" ]; then
    echo ""
    echo "Usage:"
    echo "./mVMCconfig.sh system_name"
    echo " system_name should be chosen from below:"
    echo "  gcc-fujitsu : GCC/FCC Mixed Compile"
    echo "               (Customized config for Supercomputer Fugaku.)"
    echo "   intel-impi : Intel Compiler + IntelMPI"
    echo "               (ISSP System-C users should use this config.)"
    echo "    intel-mpi : Intel Compiler + OpenMPI/MPICH2"
    echo "    aocc-aocl : AMD AOCC + AOCL + OpenMPI/MPICH2"
    echo "     gcc-aocl : GCC + AOCL + OpenMPI/MPICH2"
    echo "               (ISSP System-B users should use this config.)"
    echo "  gcc-mkl-mpi : GCC + MKL + OpenMPI/MPICH2"
    echo "  gcc-x86-mpi : GCC + OpenMPI/MPICH2 on generic x86_64"
    echo "               (Not recommended. Install MKL, AOCL or libFLAME.)"
    echo "  gcc-arm-mpi : GCC + OpenMPI/MPICH2 on Armv8-A"
    echo ""
else
    if [ ! -e ${PWD}/mVMCconfig.sh ]; then
	echo "Error: Out-out-place configuration not supported."
	echo "       Please run config in mVMC source directory."
	exit
    fi

    # Check config entries one by one.
    if [ ${1} = "gcc-fujitsu" ]; then #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cat > src/make.sys <<EOF
CC = mpifcc -Nclang
CXX = mpiFCC -Nclang
F90 = frt
FFLAGS = -DNDEBUG -DFUJITSU -Kfast
CFLAGS = -DNDEBUG -Ofast -fopenmp
CXXFLAGS = -DNDEBUG -Ofast -fopenmp
CFLAGS += -D_lapack -D_pf_block_update # -D_pfaffine
CXXFLAGS += -DBLAS_EXTERNAL

BLIS_ROOT = # Specify your installation of BLIS > v0.8.
LIBS = -SSL2 -lm -lpthread
SFMTFLAGS =
EOF
	echo "Notice: mVMC for Fugaku cannot be compiled on the login node."
	echo "        Please submit a single-node job to compile:"
	echo "        $ pjsub compile.psh"
	cat > compile.psh <<EOF
#!/bin/bash
#PJM -L "node=1"
#PJM -L "rscunit=rscunit_ft01"
#PJM -L "rscgrp=small"
#PJM -L "elapse=10:00"
#PJM -S
make -j48 mvmc
EOF
    elif [ ${1} = "intel-impi" ]; then #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cat > src/make.sys <<EOF
CC = mpiicc
CXX = mpiicpc
F90 = ifort
FFLAGS = -O3 -implicitnone
CFLAGS = -O3 -qopenmp -Wno-unknown-pragmas
CXXFLAGS = -O3 -std=gnu++14 -fpic -qopenmp
CFLAGS += -D_lapack -D_pf_block_update -D_pfaffine
CXXFLAGS += -DMKL -DBLAS_EXTERNAL

BLIS_ROOT = # Specify your installation of BLIS > v0.8.
LIBS = -mkl=sequential -lm -lpthread
SFMTFLAGS =
EOF
	echo "Notice: Recommended modules for System-C Enaga:"
	echo "        intel/20.0.1 intel-mkl/20.0.1 intel-mpi/2019.7 gcc/7.2.0"
    elif [ ${1} = "intel-mpi" ]; then #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cat > src/make.sys <<EOF
export MPICH_CC = icc
export MPICH_CXX = icpc
export OMPI_CC = icc
export OMPI_CXX = icpc
CC = mpicc
CXX = mpicxx
F90 = ifort
FFLAGS = -O3 -implicitnone
CFLAGS = -O3 -qopenmp -Wno-unknown-pragmas
CXXFLAGS = -O3 -std=gnu++14 -fpic -qopenmp
CFLAGS += -D_lapack -D_pf_block_update -D_pfaffine
CXXFLAGS += -DMKL -DBLAS_EXTERNAL

BLIS_ROOT = # Specify your installation of BLIS > v0.8.
LIBS = -mkl=sequential -lm -lpthread
SFMTFLAGS =
EOF
    elif [ ${1} = "aocc-aocl" ]; then #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cat > src/make.sys <<EOF
export MPICH_CC = clang
export MPICH_CXX = clang++
export OMPI_CC = clang
export OMPI_CXX = clang++
CC = mpicc
CXX = mpicxx
F90 = flang
FFLAGS = -O3
CFLAGS = -O3 -fopenmp
CXXFLAGS = -O3 -fPIC
CFLAGS += -D_lapack -D_pf_block_update -D_pfaffine
CXXFLAGS +=

BLIS_ROOT = # Specify your installation of BLIS > v0.8.
LIBS = -lflame -lm -lpthread
SFMTFLAGS =
EOF
    elif [ ${1} = "gcc-aocl" ]; then #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cat > src/make.sys <<EOF
export MPICH_CC = gcc
export MPICH_CXX = g++
export OMPI_CC = gcc
export OMPI_CXX = g++
CC = mpicc
CXX = mpicxx
F90 = gfortran
FFLAGS = -O3 -fimplicit-none
CFLAGS = -O3 -fopenmp
CXXFLAGS = -O3 -fPIC
CFLAGS += -D_lapack -D_pf_block_update -D_pfaffine
CXXFLAGS +=

BLIS_ROOT = # Specify your installation of BLIS > v0.8.
LIBS = -lflame -lm -lpthread
SFMTFLAGS =
EOF
    elif [ ${1} = "gcc-mkl-mpi" ]; then #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cat > src/make.sys <<EOF
export MPICH_CC = gcc
export MPICH_CXX = g++
export OMPI_CC = gcc
export OMPI_CXX = g++
CC = mpicc
CXX = mpicxx
F90 = gfortran
FFLAGS = -O3 -fimplicit-none
CFLAGS = -O3 -fopenmp
CXXFLAGS = -O3
CFLAGS += -D_lapack -D_pf_block_update -D_pfaffine
CXXFLAGS += -DBLAS_EXTERNAL

BLIS_ROOT = # Specify your installation of BLIS > v0.8.
LIBS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -lpthread
SFMTFLAGS =
EOF
    elif [ ${1} = "gcc-x86-mpi" ]; then #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cat > src/make.sys <<EOF
export MPICH_CC = gcc
export MPICH_CXX = g++
export OMPI_CC = gcc
export OMPI_CXX = g++
CC = mpicc
CXX = mpicxx
F90 = gfortran
FFLAGS = -O3 -fimplicit-none
CFLAGS = -O3 -fopenmp
CXXFLAGS = -O3
CFLAGS += -D_lapack -D_pf_block_update -D_pfaffine
CXXFLAGS +=

BLIS_ROOT = # Specify your installation of BLIS > v0.8.
LIBS = -llapack -lm -lpthread
SFMTFLAGS =
EOF
    elif [ ${1} = "gcc-arm-mpi" ]; then #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cat > src/make.sys <<EOF
export MPICH_CC = gcc
export MPICH_CXX = g++
export OMPI_CC = gcc
export OMPI_CXX = g++
CC = mpicc
CXX = mpicxx
F90 = gfortran
FFLAGS = -O3 -fimplicit-none
CFLAGS = -O3 -fopenmp
CXXFLAGS = -O3
CFLAGS += -D_lapack -D_pf_block_update -D_pfaffine
CXXFLAGS += # -DBLAS_EXTERNAL

BLIS_ROOT = # Specify your installation of BLIS > v0.8.
LIBS = -llapack -lm -lpthread
SFMTFLAGS =
EOF
    else
        echo ""
        echo "Unsupported system. Please type"
        echo "./config.sh help"
        echo ""
        exit
    fi
    cat src/make.sys
    echo 
    echo "Checking out submodules for fast Pfaffian computation."
    # Clone Pfaffine if submodule not loaded.
    if [ ! -e src/pfaffine/src ]; then
        git submodule update --init --recursive
    fi
    # Link make.sys to make.inc for Pfaffine
    ln -s $PWD/src/make.sys src/pfaffine/make.inc;

    echo
    echo "Config is not done yet."
    echo "Please edit make.sys and specify BLIS_ROOT, etc."
    echo

    cat > makefile <<EOF
help:
	@echo ""
	@echo "Usage :"
	@echo "make <entry>"
	@echo ""
	@echo "<entry> is chosen from below"
	@echo "      mvmc : Build simulator mVMC in src/ and tool/"
	@echo "     clean : Remove all generated files excepting makefile"
	@echo " veryclean : Remove all generated files including makefile"
	@echo ""

mvmc:
	make -C src/mVMC -f makefile_src
	make -C tool     -f makefile_tool

clean:
	rm -rf  src/blis-install
	make -C src/blis-build            clean
	make -C src/mVMC -f makefile_src  clean
	make -C tool     -f makefile_tool clean

veryclean:
	make clean
	rm -f src/make.sys makefile
EOF
fi
