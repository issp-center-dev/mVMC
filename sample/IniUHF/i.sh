#!/bin/sh
#QSUB -queue i18cpu
#QSUB -node 1
#QSUB -mpi 1
#QSUB -omp 24
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=00:10:00  
#PBS -N hybrid
cd ${PBS_O_WORKDIR}

MPI="mpijob"
VMC="./vmc.out"
DEF1="namelist.def"

date
echo "opt"
${MPI} ${VMC} ${DEF1} > std.out
date

