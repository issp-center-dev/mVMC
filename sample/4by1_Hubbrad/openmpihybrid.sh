#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q gr1.q
#$ -pe openmpi12 12

export OMP_NUM_THREADS=1

mpirun  -npernode 12 -np 12 ./vmc.out xnamelist.def zinitial.def

mkdir sol_1
cp * sol_1

cp zqp_opt.dat zqp_ipt.dat
mpirun  -npernode 12 -np 12 ./vmc.out xnamelist_pro.def zqp_ipt.dat

