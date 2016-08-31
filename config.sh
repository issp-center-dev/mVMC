#!/bin/bash
if [ -z ${1} ] || [ ${1} = "help" ]; then
    echo ""
    echo "Usage:"
    echo "./config.sh system_name"
    echo " system_name should be chosen from below:"
    echo "        sekirei : ISSP system-B"
    echo "            kei : Fujitsu K computer & FX10"
    echo "  openmpi-intel : OpenMPI + Intel Compiler + MKL"
    echo "    mpich-intel : MPICH + Intel Compiler + MKL"
    echo "            gnu : GNU"
    echo "        jupiter : "
    echo "        kashiwa : "
    echo "            sol : "
    echo "          reims : "
    echo "         manual : Configure manualy"
    echo ""
else
    if [ ${1} = "sekirei" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
LIB = -L\$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_sgimpt_lp64 -lpthread -lm
CFLAGS = -O3 -no-prec-div -xHost -qopenmp -Wno-unknown-pragmas
REPORT = -qopt-report-phase=openmp -qopt-report-phase=par
OPTION = -D_mpi_use
CP = cp -f -v
AR = ar rv
FORT = ifort
FFLAGS = -O3 -implicitnone -xHost
SMFTFLAGS = -O3 -no-ansi-alias -xHost -DMEXP=19937 -DHAVE_SSE2
EOF
    elif [ ${1} = "openmpi-intel" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
LIB = -L \${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_openmpi_lp64 -lpthread -lm -ldl
CFLAGS = -O3 -no-prec-div -xHost -qopenmp -Wno-unknown-pragmas -I ${MKLROOT}/include
REPORT = -qopt-report-phase=openmp -qopt-report-phase=par
OPTION = -D_mpi_use
CP = cp -f -v
AR = ar rv
FORT = ifort
FFLAGS = -O3 -implicitnone -xHost
SMFTFLAGS = -O3 -no-ansi-alias -xHost -DMEXP=19937 -DHAVE_SSE2
EOF
    elif [ ${1} = "mpich-intel" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
LIB = -L \${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl
CFLAGS = -O3 -no-prec-div -xHost -qopenmp -Wno-unknown-pragmas -I ${MKLROOT}/include
REPORT = -qopt-report-phase=openmp -qopt-report-phase=par
OPTION = -D_mpi_use
CP = cp -f -v
AR = ar rv
FORT = ifort
FFLAGS = -O3 -implicitnone -xHost
SMFTFLAGS = -O3 -no-ansi-alias -xHost -DMEXP=19937 -DHAVE_SSE2
EOF
    elif [ ${1} = "jupiter" ]; then
        cat > src/make.sys <<EOF
CC      = mpicc
LIB = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread -lm -openmp
CFLAGS = -O3 -no-prec-div -xP -Wno-unknown-pragmas
REPORT = -vec-report1
OPTION = -D_mpi_use -D_lapack
CP = cp -f -v
AR = ar rv
FORT = ifort
FFLAGS = -O3 -implicitnone
SMFTFLAGS = -O3 -no-ansi-alias -DMEXP=19937
EOF
    elif [ ${1} = "kashiwa" ]; then
        cat > src/make.sys <<EOF
CC = icc -lmpi
LIB = -L \$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64\ -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core\ -lmkl_blacs_sgimpt_lp64 -openmp -lpthread -lm
CFLAGS = -O3 -no-prec-div -xHost -openmp -Wno-unknown-pragmas
REPORT = -openmp-report1 -vec-report=1
OPTION = -D_mpi_use
CP = cp -f -v
AR = ar rv
FORT = ifort
FFLAGS = -O3 -implicitnone -xSSE2
SMFTFLAGS = -O3 -no-ansi-alias -xSSE2 -DMEXP=19937 -DHAVE_SSE2
EOF
    elif [ ${1} = "kei" ]; then
        cat > src/make.sys <<EOF
FC = mpifrtpx
CC = mpifccpx
LIB = -SCALAPACK -SSL2BLAMP
CFLAGS = -Kfast,parallel,ocl,openmp
REPORT = -Koptmsg=2
OPTION = -D_mpi_use
CP = cp -f -v
AR = ar rv
FORT = frtpx
FFLAGS = -Kfast,ocl,auto,optmsg=2 -AT
SMFTFLAGS = -Kfast,ocl,nomemalias -DMEXP=19937
EOF
    elif [ ${1} = "reims" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
LIB = -L \$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -openmp -lpthread -lm
CFLAGS = -O3 -no-prec-div -xSSE2 -openmp -Wno-unknown-pragmas
REPORT = -openmp-report1 -vec-report=1
OPTION = -D_mpi_use
CP = cp -f -v
AR = ar rv
FORT = ifort
FFLAGS = -O3 -implicitnone -xSSE2
SMFTFLAGS = -O3 -no-ansi-alias -xSSE2 -DMEXP=19937 -DHAVE_SSE2
EOF
    elif [ ${1} = "sol" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
LIB = -L \$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -openmp -lpthread -lm
CFLAGS = -O3 -no-prec-div -xHost -openmp -Wno-unknown-pragmas
REPORT = -openmp-report1 -vec-report=1
OPTION = -D_mpi_use
CP = cp -f -v
AR = ar rv
FORT = ifort
FFLAGS = -O3 -implicitnone -xSSE2
SMFTFLAGS = -O3 -no-ansi-alias -xSSE2 -DMEXP=19937 -DHAVE_SSE2
EOF
    elif [ ${1} = "gnu" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
LIB = -fopenmp  -framework Accelerate -llapack -lblas -lm
CFLAGS = -O3 -fopenmp
REPORT = 
OPTION = -D_mpi_use -D_lapack
CP = cp -f -v
AR = ar rv
FORT = gfortran
FFLAGS = -O3 -fimplicit-none
SMFTFLAGS = -DMEXP=19937
EOF
    elif [ ${1} == "manual" ]; then
        cat > src/make.sys <<EOF
CC = 
LIB = 
CFLAGS = 
REPORT = 
OPTION = 
CP = cp -f -v
AR = ar rv
FORT = 
FFLAGS = 
SMFTFLAGS = 
EOF
        echo ""
        echo "Generate src/make.sys without any setting."
        echo "Please fill src/make.sys manualy."
        echo ""
    else
        echo ""
        echo "Unsupported system. Please type"
        echo "./config.sh help"
        echo ""
        exit
    fi

    echo "cat src/make.sys"
    cat src/make.sys

    echo
    echo "config DONE"
    echo

    cat > makefile <<EOF
help:
	@echo ""
	@echo "Usage :"
	@echo "make <entry>"
	@echo ""
	@echo "<entry> is chosen from below"
	@echo "      mvmc : Build simulator mVMC in src/"
	@echo " userguide : Generate userguid_jp.pdf & userguide_en.pdf in doc/"
	@echo "     clean : Remove all generated files excepting makefile"
	@echo " veryclean : Remove all generated files including makefile"
	@echo ""

mvmc:
	cd src;make -f makefile_src

userguide:
	cd doc/jp/;make -f makefile_doc_jp;mv userguide_jp.pdf ../
	cd doc/en/;make -f makefile_doc_en;mv userguide_en.pdf ../

clean:
	cd src; make -f makefile_src clean
	cd doc/jp; make -f makefile_doc_jp clean
	cd doc/en; make -f makefile_doc_en clean
	rm -f doc/userguide_??.pdf

veryclean:
	make clean
	rm -f src/make.sys makefile
EOF
fi
