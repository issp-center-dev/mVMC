#!/bin/sh
#QSUB -queue F18acc
#QSUB -node 8
#QSUB -mpi 96
#QSUB -omp 2
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=00:30:00  
#PBS -N hybrid
cd ${PBS_O_WORKDIR}
. /etc/profile.d/modules.sh
module load intel/15.0.3.187 mpt/2.12 

MPI="mpijob"
VMC="./vmc.out -m 4 DirList.def"
DEF2="xnamelist_pro.def  zinitial.def"
DEFaft="xnamelist_aft.def  zqp_opt.dat"
DEFSC2="xnamelist_aft_SC2.def  zqp_opt.dat"



date
echo "opt_pro"
${MPI} ${VMC} ${DEF2}
date

date
echo "aft"
${MPI} ${VMC} ${DEFaft}
date

date
echo "SC2"
${MPI} ${VMC} ${DEFSC2}
date
