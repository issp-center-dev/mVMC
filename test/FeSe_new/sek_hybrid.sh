#!/bin/sh
#QSUB -queue i9acc
#QSUB -node 8
#QSUB -mpi  16
#QSUB -omp  12
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=00:30:00  
#PBS -N hybrid
cd ${PBS_O_WORKDIR}

MPI="mpijob"
VMC="./vmc.out -F 5"
#DEF="xnamelist_pro.def zinitial.def"
DEF="xnamelist_pro.def"

date
echo "opt_pro"
${MPI} ${VMC} ${DEF}
date

