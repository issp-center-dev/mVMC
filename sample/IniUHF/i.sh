#!/bin/sh
#QSUB -queue i18cpu
#QSUB -node 1
#QSUB -mpi 24
#QSUB -omp 1
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=00:30:00  
#PBS -N hybrid
cd ${PBS_O_WORKDIR}

MPI="mpijob"
VMC="./vmc.out -F 5"
DEF1="namelist.def"

date
echo "opt"
${MPI} ${VMC} ${DEF1}
date

