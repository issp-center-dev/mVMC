#!/bin/sh
#QSUB -queue F144cpu
#QSUB -node 72
#QSUB -mpi  72
#QSUB -omp  24
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=12:00:00  
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

