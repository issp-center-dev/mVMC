#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q gr1.q
#$ -pe openmpi 48
export OMP_NUM_THREADS=1

MPI="mpirun -np 12"
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

