#!/bin/bash
if [ -z ${1} ] || [ ${1} = "help" ]; then
    echo ""
    echo "Usage:"
    echo "./config.sh system_name"
    echo " system_name should be chosen from below:"
    echo "        sekirei : ISSP system-B"
    echo "            kei : Fujitsu K computer & FX10"
    echo "  intel-openmpi : Intel Compiler + OpenMPI"
    echo "    intel-mpich : Intel Compiler + MPICH2 (or IntelMPI)"
    echo "    gcc-openmpi : GCC + OpenMPI"
    echo "  gcc-mpich-mkl : GCC + MPICH + MKL"
    echo "        kashiwa : Remain for compatibility"
    echo "        jupiter : Remain for compatibility"
    echo "            sol : Remain for compatibility"
    echo "          reims : Remain for compatibility"
    echo ""
else
    if [ ${1} = "sekirei" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -O3 -xHost -qopenmp -no-prec-div -Wno-unknown-pragmas
FFLAGS = -O3 -xHost -qopenmp -implicitnone
LIBS = -qopenmp -L \$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_sgimpt_lp64 -lpthread -lm
SFMTFLAGS = -no-ansi-alias -DHAVE_SSE2
EOF
    elif [ ${1} = "kashiwa" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -O3 -xHost -openmp -no-prec-div -Wno-unknown-pragmas
FFLAGS = -O3 -xHost -openmp -implicitnone
LIBS = -openmp -L \$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_sgimpt_lp64 -lpthread -lm -openmp
SFMTFLAGS = -no-ansi-alias -DHAVE_SSE2
EOF
    elif [ ${1} = "kei" ]; then
        cat > src/make.sys <<EOF
CC = mpifccpx
F90 = mpifrtpx
CFLAGS = -Kfast,parallel,ocl,openmp
FFLAGS = -Kfast,parallel,ocl,openmp,ocl,auto -AT -Cpp -D FUJITSU
LIBS = -SCALAPACK -SSL2BLAMP
SFMTFLAGS = -Kfast,ocl,nomemalias
EOF
    elif [ ${1} = "intel-openmpi" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -O3 -xHost -qopenmp -no-prec-div -Wno-unknown-pragmas -I \${MKLROOT}/include
FFLAGS = -O3 -xHost -qopenmp -implicitnone
LIBS = -qopenmp -L \${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_openmpi_lp64 -lpthread -lm -ldl
SFMTFLAGS = -no-ansi-alias -DHAVE_SSE2
EOF
    elif [ ${1} = "sol" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -O3 -xHost -openmp -no-prec-div -Wno-unknown-pragmas
FFLAGS = -O3 -xHost -openmp -implicitnone
LIBS = -openmp -L \$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lpthread -lm -openmp
SFMTFLAGS = -no-ansi-alias -DHAVE_SSE2
EOF
    elif [ ${1} = "intel-mpich" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -O3 -xHost -qopenmp -no-prec-div -Wno-unknown-pragmas -I \${MKLROOT}/include
FFLAGS = -O3 -xHost -qopenmp -implicitnone
LIBS = -qopenmp -L \${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl
SFMTFLAGS = -no-ansi-alias -DHAVE_SSE2
EOF
    elif [ ${1} = "reims" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -O3 -xHost -openmp -no-prec-div -Wno-unknown-pragmas
FFLAGS = -O3 -xHost -openmp -implicitnone
LIBS = -openmp -L \$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -openmp
SFMTFLAGS = -no-ansi-alias -DHAVE_SSE2
EOF
    elif [ ${1} = "gcc-openmpi" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = gfortran
CFLAGS = -O3 -fopenmp
FFLAGS = -O3 -fopenmp -fimplicit-none
LIBS = -fopenmp -lblacs-openmpi -lscalapack-openmpi -llapack -lblas -lm
SFMTFLAGS =
EOF
    elif [ ${1} = "gcc-mpich-mkl" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = gfortran
CFLAGS = -O3 -fopenmp -I\${MKLROOT}/include
FFLAGS = -O3 -fopenmp -fimplicit-none
LIBS = -fopenmp -L\${MKLROOT}/lib -Wl,-rpath,\${MKLROOT}/lib -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_mpich_lp64 -lpthread -lm -ldl
SFMTFLAGS =
EOF
    elif [ ${1} = "jupiter" ]; then
        cat > src/make.sys <<EOF
CC = mpicc
F90 = mpif90
CFLAGS = -O3 -no-prec-div -xP -Wno-unknown-pragmas -D_lapack
FFLAGS = -O3 -implicitnone
LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread -lm -openmp
SFMTFLAGS = -no-ansi-alias
EOF
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
	@echo "      mvmc : Build simulator mVMC in src/ and tool/"
	@echo " userguide : Generate userguid_jp.pdf & userguide_en.pdf in doc/"
	@echo "     clean : Remove all generated files excepting makefile and doc/"
	@echo " veryclean : Remove all generated files including makefile and doc/"
	@echo ""

mvmc:
	cd src;make -f makefile_src
	cd tool;make -f makefile_tool

userguide:
	cd doc/jp/;make -f makefile_doc_jp;mv userguide_jp.pdf ../
	cd doc/en/;make -f makefile_doc_en;mv userguide_en.pdf ../
	cd doc/fourier/ja; make html latexpdfja
	cd doc/fourier/en; make html latexpdfja

clean:
	cd src; make -f makefile_src clean
	cd tool; make -f makefile_tool clean

veryclean:
	make clean
	cd doc/jp; make -f makefile_doc_jp clean
	cd doc/en; make -f makefile_doc_en clean
	cd doc/fourier/ja; make clean
	cd doc/fourier/en; make clean
	rm -f doc/userguide_??.pdf
	rm -f src/make.sys makefile
EOF
fi
